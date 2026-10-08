// Offline, full-data lookahead. Requires the isolated instrumentation patch.
#define SPLINESKETCH_TESTING
#include <splinesketch/paper_splinesketch.hpp>
#include "../accuracy_workload.hpp"
#include "../../tests/data/paper_width_counterexample.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>

using Sketch = splinesketch::CertifiedPaperSplineSketch;
static std::uint64_t bits(double value) { std::uint64_t r; std::memcpy(&r, &value, sizeof r); return r; }
static double from_bits(std::uint64_t value) { double r; std::memcpy(&r, &value, sizeof r); return r; }
struct Error {
  double high = 0, low = 0;
  bool operator<(const Error& other) const { return std::tie(high, low) < std::tie(other.high, other.low); }
};
// Error-free TwoSum: retain the exact difference of binary64 estimate and
// exactly representable integer truth, including the subtraction's rounding.
static Error rank_error(double estimate, std::uint64_t truth) {
  const double b = -static_cast<double>(truth), high = estimate + b;
  const double v = high - estimate;
  const double low = (estimate - (high - v)) + (b - v);
  return high < 0 || (high == 0 && low < 0) ? Error{-high, -low} : Error{high, low};
}
struct Hash {
  std::uint64_t value = 0xcbf29ce484222325ULL;
  void add(std::uint64_t x) { value = (value ^ x) * 0x100000001b3ULL; }
  void floating(long double x) {
    const double high = static_cast<double>(x);
    add(bits(high)); add(bits(std::isfinite(high) ? static_cast<double>(x - high) : 0));
  }
};
struct Alternative { std::uint64_t x; };
struct Event {
  std::size_t index; std::uint64_t state, count, original, selected, split;
  unsigned kind; std::vector<Alternative> alternatives;
};
struct Decision { std::uint64_t state, x; };
using Plan = std::map<std::size_t, Decision>;
struct Context { const Plan* plan; std::vector<Event> events; std::size_t forced = 0; };
static Context* active = nullptr;
namespace splinesketch {
struct PaperCertificateInspector {
  static std::pair<double, double> mass_interval(const Sketch& sketch) {
    std::pair<double, double> result{-1, 1}; long double largest = -1;
    for (std::size_t i = 1; i < sketch.nodes_.size(); ++i)
      if (std::nextafter(sketch.nodes_[i - 1].x, INFINITY) < sketch.nodes_[i].x && sketch.nodes_[i].mass > largest) {
        largest = sketch.nodes_[i].mass; result = {sketch.nodes_[i - 1].x, sketch.nodes_[i].x};
      }
    return result;
  }
  static void hook(bool enabled) { Sketch::oracle_join_hook_ = enabled ? select : nullptr; }
  static std::size_t select(const Sketch& sketch, std::size_t original, std::size_t split, unsigned kind) {
    if (!active || sketch.use_heap_backend()) throw std::runtime_error("oracle requires a scanner context");
    Hash h; h.add(sketch.count_); h.add(sketch.epoch_end_); h.add(sketch.capacity_); h.floating(sketch.factor_);
    for (const auto& n : sketch.nodes_) {
      h.add(bits(n.x)); h.floating(n.mass); h.floating(n.prefix); h.floating(n.slope);
      h.add(n.lower); h.add(n.upper); h.add(n.atom); h.add(n.protected_threshold);
    }
    for (double x : sketch.buffer_) h.add(bits(x));
    for (const auto& item : sketch.heavy_) { h.add(bits(item.first)); h.add(item.second.exact); h.add(item.second.residual); }
    Event event{active->events.size(), h.value, sketch.count_, bits(sketch.nodes_[original].x),
                bits(sketch.nodes_[original].x), split == Sketch::npos ? 0 : bits(sketch.nodes_[split].x), kind, {}};
    std::size_t selected = original;
    const auto override = active->plan->find(event.index);
    if (override != active->plan->end() && override->second.state != event.state)
      throw std::runtime_error("oracle override reached a different prefix state");
    bool found_original = false, found_override = override == active->plan->end();
    for (std::size_t i = 1; i + 1 < sketch.nodes_.size(); ++i) {
      if (i == split || i + 1 == split || sketch.nodes_[i].protected_threshold ||
          !(sketch.nodes_[i].mass + sketch.nodes_[i + 1].mass <= 0.75L * sketch.bound())) continue;
      const auto x = bits(sketch.nodes_[i].x); event.alternatives.push_back({x});
      found_original |= i == original;
      if (override != active->plan->end() && override->second.x == x) { selected = i; found_override = true; }
    }
    if (!found_original || !found_override) throw std::runtime_error("oracle selected an ineligible or overlapping join");
    event.selected = bits(sketch.nodes_[selected].x);
    active->forced += override != active->plan->end(); active->events.push_back(std::move(event));
    return selected;
  }
};
}
struct HookSession {
  explicit HookSession(Context& context) {
    if (active) throw std::runtime_error("nested oracle replay");
    active = &context; splinesketch::PaperCertificateInspector::hook(true);
  }
  ~HookSession() { splinesketch::PaperCertificateInspector::hook(false); active = nullptr; }
};
struct Query { double x; std::uint64_t truth; };
struct Metrics {
  Error error; double worst_x = 0, estimate = 0, uniform = 0;
  std::uint64_t truth = 0, width = 0, hash = 0;
};
static std::vector<Query> query_set(const std::vector<double>& data) {
  auto sorted = data; std::sort(sorted.begin(), sorted.end());
  auto xs = accuracy_workload::make_queries(sorted);
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    xs.push_back(sorted[i]); xs.push_back(std::nextafter(sorted[i], -INFINITY)); xs.push_back(std::nextafter(sorted[i], INFINITY));
    if (i && sorted[i] > sorted[i - 1]) xs.push_back(sorted[i - 1] + (sorted[i] - sorted[i - 1]) / 2);
  }
  const auto ordered = [](double x) { const auto b = bits(x); return b >> 63 ? ~b : b ^ (std::uint64_t{1} << 63); };
  const auto decode = [](std::uint64_t x) { return from_bits(x >> 63 ? x ^ (std::uint64_t{1} << 63) : ~x); };
  const auto lo = ordered(sorted.front()), span = ordered(sorted.back()) - lo;
  for (unsigned i = 0; i <= 1024; ++i) {
    xs.push_back(sorted.front() + (sorted.back() - sorted.front()) * i / 1024);
    xs.push_back(decode(lo + span / 1024 * i + span % 1024 * i / 1024));
  }
  xs.push_back(-INFINITY); xs.push_back(INFINITY);
  std::sort(xs.begin(), xs.end()); xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  std::vector<Query> result;
  for (double x : xs) {
    if (std::isnan(x)) throw std::logic_error("invalid oracle query");
    result.push_back({x, static_cast<std::uint64_t>(std::upper_bound(sorted.begin(), sorted.end(), x) - sorted.begin())});
  }
  return result;
}
static std::uint64_t certificate_checks = 0;
static Metrics measure(const Sketch& sketch, const std::vector<Query>& queries, bool exact_output = false) {
  Metrics result; Hash hash;
  result.uniform = sketch.max_rank_error(); result.width = sketch.max_rank_uncertainty();
  for (const auto& query : queries) {
    const auto answer = sketch.rank_with_error(query.x);
    const auto error = rank_error(answer.estimate, query.truth);
    if (answer.lower_rank > query.truth || answer.upper_rank < query.truth ||
        Error{answer.max_error, 0} < error || Error{result.uniform, 0} < error)
      throw std::runtime_error("invalid oracle certificate");
    if (result.error < error) { result.error = error; result.worst_x = query.x; result.estimate = answer.estimate; result.truth = query.truth; }
    hash.add(bits(query.x)); hash.add(bits(answer.estimate)); hash.add(answer.lower_rank); hash.add(answer.upper_rank); hash.add(bits(answer.max_error));
    if (exact_output)
      std::cout << bits(answer.estimate) << ',' << query.truth << ',' << answer.lower_rank << ',' << answer.upper_rank
                << ',' << bits(answer.max_error) << ',' << bits(result.uniform) << ',' << bits(error.high) << ',' << bits(error.low) << '\n';
    ++certificate_checks;
  }
  result.hash = hash.value; return result;
}
struct Replay { bool success = false; Metrics metrics; std::vector<Event> events; std::string failure; };
static Replay replay(const std::vector<double>& data, const std::vector<Query>& queries, unsigned k,
                     const Plan& plan, bool exact_output = false, bool hooks = true) {
  Context context{&plan, {}, 0}; Replay result;
  try {
    Sketch sketch(k);
    if (hooks) {
      HookSession session(context);
      for (double x : data) sketch.add(x);
      sketch.consolidate();
    } else {
      for (double x : data) sketch.add(x);
      sketch.consolidate();
    }
    if (context.forced != plan.size()) throw std::runtime_error("unreached oracle override");
    result.metrics = measure(sketch, queries, exact_output);
    result.events = std::move(context.events); result.success = true;
  } catch (const std::logic_error& error) { result.failure = error.what(); }
  return result;
}
static void emit(const std::string& phase, std::size_t event, std::uint64_t choice, bool original, const Replay& result,
                 std::size_t queries, std::size_t n) {
  const auto& m = result.metrics;
  std::cout << phase << ',' << event << ',' << choice << ',' << original << ',' << result.success << ',' << queries << ','
            << m.error.high << ',' << m.error.low << ',' << 100 * m.error.high / static_cast<double>(n) << ','
            << bits(m.worst_x) << ',' << bits(m.estimate) << ',' << m.truth << ',' << m.width << ',' << m.uniform
            << ',' << m.hash << ',' << result.events.size() << ',' << result.failure << '\n';
}
static std::vector<double> input(const std::string& path) {
  std::vector<double> data;
  if (path.empty()) for (auto x : paper_width_counterexample) data.push_back(from_bits(x));
  else { std::ifstream stream(path); std::uint64_t x; while (stream >> x) data.push_back(from_bits(x)); }
  if (data.empty()) throw std::invalid_argument("empty input");
  for (double x : data) if (!std::isfinite(x)) throw std::invalid_argument("nonfinite input");
  return data;
}
static void write_plan(const std::string& path, const Plan& plan) {
  std::ofstream stream(path);
  for (const auto& item : plan) stream << item.first << ' ' << item.second.state << ' ' << item.second.x << '\n';
  if (!stream) throw std::runtime_error("cannot save oracle plan");
}
int main(int argc, char** argv) {
  std::string data_path, output_prefix = "/tmp/paper-join-oracle", plan_path;
  unsigned k = 32, generate_seed = 0, observations = 640;
  bool exact_output = false, evaluate_only = false, generate = false;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--exact-queries") exact_output = true;
    else if (option == "--evaluate-only") evaluate_only = true;
    else if (i + 1 < argc && option == "--data") data_path = argv[++i];
    else if (i + 1 < argc && option == "--capacity") k = static_cast<unsigned>(std::stoul(argv[++i]));
    else if (i + 1 < argc && option == "--output-prefix") output_prefix = argv[++i];
    else if (i + 1 < argc && option == "--plan") plan_path = argv[++i];
    else if (i + 1 < argc && option == "--generate-mass") { generate = true; generate_seed = static_cast<unsigned>(std::stoul(argv[++i])); }
    else if (i + 1 < argc && option == "--observations") observations = static_cast<unsigned>(std::stoul(argv[++i]));
    else throw std::invalid_argument("invalid oracle arguments");
  }
  if (k >= 512) throw std::invalid_argument("oracle instrumentation supports scanner capacities only");
  if (generate) {
    Sketch sketch(k); std::mt19937_64 random(generate_seed); std::pair<double, double> interval{-1, 1};
    for (unsigned i = 0; i < observations; ++i) {
      if (i % (5 * k) == 0) interval = splinesketch::PaperCertificateInspector::mass_interval(sketch);
      const double u = static_cast<double>(random() % 1000001) / 1000001;
      double x = interval.first + (interval.second - interval.first) * (0.0001 + u * 0.0001);
      if (x == interval.first) x = std::nextafter(x, interval.second);
      sketch.add(x); std::cout << bits(x) << '\n';
    }
    return 0;
  }
  const auto data = input(data_path); const auto queries = query_set(data);
  Plan plan;
  if (!plan_path.empty()) {
    std::ifstream stream(plan_path); std::size_t event; Decision d;
    while (stream >> event >> d.state >> d.x) plan.emplace(event, d);
    if (!stream.eof()) throw std::runtime_error("cannot read oracle plan");
  }
  if (exact_output) {
    const auto result = replay(data, queries, k, plan, true);
    if (!result.success) throw std::runtime_error(result.failure);
    std::cout << "Oracle exact: " << queries.size() << " queries\n"; return 0;
  }
  std::cout << std::setprecision(17)
            << "phase,event,choice_bits,is_original,success,queries,error_high,error_low,error_percent,worst_query_bits,estimate_bits,truth,width,uniform,query_hash,joins,failure\n";
  const auto baseline = replay(data, queries, k, {});
  const auto unhooked = replay(data, queries, k, {}, false, false);
  if (!baseline.success || !unhooked.success || baseline.metrics.hash != unhooked.metrics.hash ||
      baseline.metrics.uniform != unhooked.metrics.uniform) throw std::logic_error("oracle logging changed the baseline");
  emit("baseline", 0, 0, true, baseline, queries.size(), data.size());
  if (evaluate_only) return 0;
  std::ofstream events(output_prefix + "-events.csv");
  events << "event,count,kind,state_hash,original_bits,split_bits,legal_choices\n";
  for (const auto& event : baseline.events) {
    events << event.index << ',' << event.count << ',' << event.kind << ',' << event.state << ',' << event.original << ','
           << event.split << ',' << event.alternatives.size() << '\n';
    for (const auto& alternative : event.alternatives) {
      const Plan single{{event.index, {event.state, alternative.x}}};
      const auto arm = replay(data, queries, k, single);
      emit("single", event.index, alternative.x, alternative.x == event.original, arm, queries.size(), data.size());
      if (alternative.x == event.original && (!arm.success || arm.metrics.hash != baseline.metrics.hash))
        throw std::logic_error("original oracle arm did not reproduce baseline");
    }
  }
  if (!events) throw std::runtime_error("cannot save oracle events");
  // Sequential lookahead: each arm completes the remaining stream using the
  // default rule; lock the best current decision, then inspect the new path.
  auto current = baseline; Plan greedy;
  for (std::size_t at = 0; at < current.events.size(); ++at) {
    if (at > 10000) throw std::logic_error("oracle decision limit exceeded");
    const auto event = current.events[at]; auto best = current; Plan best_plan = greedy;
    for (const auto& alternative : event.alternatives) {
      auto candidate = greedy; candidate[at] = {event.state, alternative.x};
      const auto arm = replay(data, queries, k, candidate);
      emit("lookahead", at, alternative.x, alternative.x == event.original, arm, queries.size(), data.size());
      if (arm.success && arm.metrics.error < best.metrics.error) { best = arm; best_plan = std::move(candidate); }
    }
    if (current.metrics.error < best.metrics.error) throw std::logic_error("lookahead worsened terminal objective");
    greedy = std::move(best_plan); current = std::move(best);
  }
  emit("oracle", greedy.size(), 0, false, current, queries.size(), data.size());
  write_plan(output_prefix + "-plan.txt", greedy);
  std::cerr << "Oracle: " << baseline.events.size() << " baseline joins; " << greedy.size()
            << " changed choices; " << certificate_checks << " exact-error certificate checks\n";
}
