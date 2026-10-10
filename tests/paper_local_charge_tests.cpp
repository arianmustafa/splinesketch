// Compiled against a temporary, observer-instrumented copy of the headers by
// paper_local_charge_checks.py. The production headers are never modified.
#include <cstddef>
#include <cstdint>
#include <type_traits>
#define SPLINESKETCH_TESTING
namespace splinesketch {
struct PaperCertificateInspector {
  template<class S, class T> static void begin(const S&, const T&);
  template<class S, class T> static void born(const S&, const T&, double, double, double,
                                             long double, long double, long double, long double);
  template<class S> static void joined(const S&, long double, long double);
  template<class S, class T, class N> static void end(const S&, const T&, const N&);
  template<class S> static void merged(const S&, const S&, const S&);
  template<class S> static void completed(const S&);
  static void operation_cases();
  static void initialization_cases();
  static void selected(std::size_t, std::uint64_t, std::uint64_t, std::size_t);
};
}
#include <splinesketch/paper_splinesketch.hpp>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include "paper_resize_history_fixture.hpp"

namespace {
struct Birth {
  double x, a, b, live_a, live_b;
  std::uint64_t epoch, count;
  bool root;
};
std::vector<Birth> births;
std::map<double, std::size_t> live, frozen;
// Actual split ancestry is separate from frozen certificate ancestry.
std::map<double, std::size_t> split_live, split_frozen;
std::size_t split_next = 0;
bool quiet = false, geometry_only = false, initialization_mode = false;
bool operation_mode = false;
std::uint64_t rounding_q = 0;
std::uint64_t violations = 0, passes = 0;
std::uint64_t pass_splits = 0, pass_joins = 0;
std::vector<std::pair<double, long double>> pass_births;
std::uint64_t active_state = 0, next_state = 0;
std::map<double, std::uint64_t> absorbed_truth, full_truth;
std::vector<double> extra_queries;

std::uint64_t exact_rank(const std::map<double, std::uint64_t>& truth, double x) {
  std::uint64_t rank = 0;
  for (auto it = truth.begin(); it != truth.end() && it->first <= x; ++it) rank += it->second;
  return rank;
}

std::uint64_t bits(double x) {
  std::uint64_t result; std::memcpy(&result, &x, sizeof result); return result;
}
void dyadic(long double value) {
  assert(std::isfinite(value) && value >= 0 && value < 0x1p64L);
  const auto whole = static_cast<std::uint64_t>(value);
  auto remainder = value - static_cast<long double>(whole);
  std::string digits;
  while (remainder) {
    assert(digits.size() < 2048);
    remainder *= 2;
    const bool one = remainder >= 1;
    digits.push_back(one ? '1' : '0');
    if (one) remainder -= 1;
  }
  std::cout << '"' << whole << ':' << digits << '"';
}
template<class N>
void grid(const N& nodes, bool with_truth = false, bool with_bounds = false) {
  std::cout << '[';
  bool comma = false;
  for (const auto& node : nodes) {
    if (comma) std::cout << ',';
    comma = true;
    std::cout << '[' << bits(node.x) << ','; dyadic(node.prefix);
    if (with_truth) std::cout << ',' << exact_rank(absorbed_truth, node.x);
    if (with_bounds) std::cout << ',' << node.lower << ',' << node.upper << ',' << node.atom;
    std::cout << ']';
  }
  std::cout << ']';
}
template<class N>
void ledger_state(const N& nodes, std::uint64_t count, const char* stage = "complete") {
  if (!operation_mode) return;
  std::cout << "{\"kind\":\"ledger_state\",\"q\":" << rounding_q
            << ",\"id\":" << active_state << ",\"stage\":\"" << stage << '"'
            << ",\"count\":" << count << ",\"absorbed\":"
            << (nodes.empty() ? 0 : nodes.back().upper) << ",\"top\":";
  dyadic(nodes.empty() ? 0 : nodes.back().prefix);
  std::cout << ",\"nodes\":"; grid(nodes, true, true);
  assert(exact_rank(absorbed_truth, INFINITY) == (nodes.empty() ? 0 : nodes.back().upper));
  std::cout << "}\n";
}
std::size_t root(double x, std::uint64_t epoch, std::uint64_t count) {
  const auto id = births.size();
  births.push_back({x, 0, 0, 0, 0, epoch, count, true}); live[x] = id;
  if (!quiet)
    std::cout << "{\"kind\":\"root\",\"id\":" << id << ",\"x\":" << bits(x) << "}\n";
  return id;
}
std::size_t split_root(double x) {
  const auto id = split_next++;
  split_live[x] = id;
  if (!quiet)
    std::cout << "{\"kind\":\"live_root\",\"id\":" << id << ",\"x\":" << bits(x) << "}\n";
  return id;
}
std::size_t split_identity(double x) {
  const auto found = split_live.find(x);
  return found == split_live.end() ? split_root(x) : found->second;
}
}

template<class S, class T>
void splinesketch::PaperCertificateInspector::begin(const S& sketch, const T& before) {
  if constexpr (std::is_same_v<S, CertifiedPaperSplineSketch>) {
    pass_splits = pass_joins = 0;
    if (operation_mode) {
      pass_births.clear();
      for (const auto& item : before.items) absorbed_truth[item.first] += item.second;
      return;
    }
    frozen.clear();
    split_frozen.clear();
    for (const auto& node : before.nodes) {
      auto it = live.find(node.x);
      const auto id = it == live.end() ? root(node.x, sketch.epoch_end_, sketch.count_) : it->second;
      frozen.emplace(node.x, id);
      split_frozen.emplace(node.x, split_identity(node.x));
    }
  }
}

template<class S, class T>
void splinesketch::PaperCertificateInspector::born(const S& sketch, const T& before, double x,
                                                double live_a, double live_b, long double chosen,
                                                long double low, long double high, long double mass) {
  if constexpr (std::is_same_v<S, CertifiedPaperSplineSketch>) {
    assert(live_a < x && x < live_b);
    ++pass_splits;
    if (operation_mode) {
      pass_births.emplace_back(x, S::spline_rank(before.nodes, x));
      return;
    }
    if (!quiet) {
      std::cout << "{\"kind\":\"mass_split\",\"n\":" << sketch.count_
                << ",\"initial\":" << (frozen.empty() ? "true" : "false")
                << ",\"factor\":";
      dyadic(sketch.factor_);
      std::cout << ",\"low\":"; dyadic(low);
      std::cout << ",\"high\":"; dyadic(high);
      std::cout << ",\"mass\":"; dyadic(mass);
      std::cout << "}\n";
      if (frozen.empty()) split_root(x); // Initialization has exact rank roots.
      else {
        const auto left_id = split_identity(live_a), right_id = split_identity(live_b);
        const auto id = split_next++;
        split_live[x] = id;
        std::cout << "{\"kind\":\"live_birth\",\"id\":" << id << ",\"x\":" << bits(x)
                  << ",\"a\":" << bits(live_a) << ",\"b\":" << bits(live_b)
                  << ",\"lower\":" << left_id << ",\"upper\":" << right_id
                  << ",\"epoch\":" << sketch.epoch_end_ << ",\"chosen\":";
        dyadic(chosen);
        std::cout << "}\n";
      }
    }
    const auto right = frozen.lower_bound(x);
    if (right != frozen.end() && right->first == x) { live[x] = right->second; return; }
    if (right == frozen.begin() || right == frozen.end()) {
      root(x, sketch.epoch_end_, sketch.count_); return;
    }
    const auto left = std::prev(right);
    const auto id = births.size();
    const Birth child{x, left->first, right->first, live_a, live_b,
                      sketch.epoch_end_, sketch.count_, false};
    for (auto parent : {left->second, right->second}) {
      const auto& old = births[parent];
      if (old.root || old.epoch != child.epoch) continue;
      if (8 * (static_cast<long double>(child.b) - child.a) >
          7 * (static_cast<long double>(old.b) - old.a)) ++violations;
    }
    births.push_back(child); live[x] = id;
    if (!quiet)
      std::cout << "{\"kind\":\"birth\",\"id\":" << id << ",\"x\":" << bits(x)
                << ",\"a\":" << bits(child.a) << ",\"b\":" << bits(child.b)
                << ",\"lower\":" << left->second << ",\"upper\":" << right->second
                << ",\"live_a\":" << bits(live_a) << ",\"live_b\":" << bits(live_b)
                << ",\"epoch\":" << child.epoch << ",\"n\":" << child.count << "}\n";
  }
}

template<class S>
void splinesketch::PaperCertificateInspector::joined(const S& sketch, long double sum, long double limit) {
  if constexpr (std::is_same_v<S, CertifiedPaperSplineSketch>) {
    ++pass_joins;
    if (operation_mode) return;
    if (!quiet) {
      std::cout << "{\"kind\":\"mass_join\",\"n\":" << sketch.count_ << ",\"factor\":";
      dyadic(sketch.factor_);
      std::cout << ",\"sum\":"; dyadic(sum);
      std::cout << ",\"limit\":"; dyadic(limit);
      std::cout << "}\n";
    }
  }
}

template<class S, class T, class N>
void splinesketch::PaperCertificateInspector::end(const S& sketch, const T& before, const N& nodes) {
  if constexpr (std::is_same_v<S, CertifiedPaperSplineSketch>) {
    ++passes;
    if (operation_mode) {
      const auto dimension = std::max({sketch.capacity_, before.nodes.size(), nodes.size()}) + 3;
      const auto cost = 10 * dimension * (1 + pass_splits + pass_joins);
      rounding_q += cost;
      std::cout << "{\"kind\":\"rank_pass\",\"cost\":" << cost
                << ",\"splits\":" << pass_splits << ",\"released\":"
                << (before.sums.empty() ? 0 : before.sums.back()) << ",\"source\":";
      grid(before.nodes, false, true);
      std::cout << ",\"births\":[";
      for (std::size_t i = 0; i < pass_births.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << '[' << bits(pass_births[i].first) << ',';
        dyadic(pass_births[i].second);
        std::cout << ']';
      }
      std::cout << ']';
      // Exact checker reconstructs the transfer defect from finite stored
      // values and exact incoming counts, without consulting absorbed truth.
      std::cout << ",\"transfers\":[";
      bool comma = false;
      for (const auto& node : nodes) {
        if (comma) std::cout << ',';
        comma = true;
        std::cout << '[' << bits(node.x) << ','; dyadic(node.prefix);
        std::cout << ','; dyadic(S::spline_rank(before.nodes, node.x));
        std::cout << ',' << before.integer_rank(node.x) << ']';
      }
      std::cout << ']';
      std::cout << "}\n";
      ledger_state(nodes, sketch.count_, "pass");
      return;
    }
    if (quiet) return;
    for (const auto& node : nodes)
      if (live.find(node.x) == live.end()) root(node.x, sketch.epoch_end_, sketch.count_);
    for (const auto& node : nodes) split_identity(node.x);
    if (geometry_only) return;
    std::cout << "{\"kind\":\"pass\",\"n\":" << sketch.count_ << ",\"epoch\":"
              << sketch.epoch_end_ << ",\"splits\":" << pass_splits
              << ",\"joins\":" << pass_joins << ",\"factor\":";
    dyadic(sketch.factor_);
    std::cout << ",\"source\":[";
    const auto write = [](const auto& sequence, const auto& ids) {
      bool comma = false;
      for (const auto& node : sequence) {
        if (comma) std::cout << ',';
        comma = true;
        std::cout << '[' << bits(node.x) << ','; dyadic(node.prefix);
        std::cout << ',' << node.lower << ',' << node.upper << ',' << node.atom << ','
                  << ids.at(node.x) << ']';
      }
    };
    write(before.nodes, frozen);
    std::cout << "],\"out\":["; write(nodes, live);
    std::cout << "],\"masses\":[";
    bool mass_comma = false;
    for (const auto& node : nodes) {
      if (mass_comma) std::cout << ',';
      mass_comma = true; dyadic(node.mass);
    }
    const auto write_ids = [](const auto& sequence, const auto& ids) {
      bool comma = false;
      for (const auto& node : sequence) {
        if (comma) std::cout << ',';
        comma = true; std::cout << ids.at(node.x);
      }
    };
    std::cout << "],\"source_live_ids\":["; write_ids(before.nodes, split_frozen);
    std::cout << "],\"out_live_ids\":["; write_ids(nodes, split_live);
    std::cout << "],\"items\":[";
    bool item_comma = false;
    for (const auto& item : before.items) {
      if (item_comma) std::cout << ',';
      item_comma = true;
      std::cout << '[' << bits(item.first) << ',' << item.second << ']';
    }
    std::cout << "],\"protected\":[";
    bool protected_comma = false;
    for (const auto& node : nodes) {
      if (!node.protected_threshold) continue;
      if (protected_comma) std::cout << ',';
      protected_comma = true;
      std::cout << bits(node.x);
    }
    std::cout << "],\"evaluations\":[";
    bool comma = false;
    for (const auto& node : nodes) {
      if (comma) std::cout << ',';
      comma = true;
      std::cout << '['; dyadic(S::spline_rank(before.nodes, node.x));
      std::cout << ',' << before.integer_rank(node.x) << ',' << before.atom(node.x) << ']';
    }
    std::cout << "]}\n";
  }
}

template<class S>
void splinesketch::PaperCertificateInspector::merged(const S& result, const S& left, const S& right) {
  if constexpr (std::is_same_v<S, CertifiedPaperSplineSketch>) {
    if (!operation_mode) return;
    const auto cost = 10 * (result.nodes_.size() + 3);
    rounding_q += cost;
    std::cout << "{\"kind\":\"merge_samples\",\"cost\":" << cost << ",\"left\":";
    grid(left.nodes_, false, true);
    std::cout << ",\"right\":"; grid(right.nodes_, false, true);
    std::cout << ",\"cuts\":[";
    for (std::size_t i = 0; i < result.nodes_.size(); ++i) {
      if (i) std::cout << ',';
      std::cout << bits(result.nodes_[i].x);
    }
    std::cout << "],\"samples\":[";
    bool comma = false;
    for (const auto& node : result.nodes_) {
      if (comma) std::cout << ',';
      comma = true;
      std::cout << '[';
      const auto left_sample = S::spline_rank(left.nodes_, node.x);
      const auto right_sample = S::spline_rank(right.nodes_, node.x);
      dyadic(left_sample + right_sample);
      std::cout << ','; dyadic(node.mass);
      std::cout << ','; dyadic(node.prefix);
      std::cout << ','; dyadic(left_sample);
      std::cout << ','; dyadic(right_sample);
      std::cout << ']';
    }
    std::cout << "]}\n";
    ledger_state(result.nodes_, result.count_, "merge");
  }
}

template<class S>
void splinesketch::PaperCertificateInspector::completed(const S& sketch) {
  auto reconstructed = absorbed_truth;
  for (double x : sketch.buffer_) ++reconstructed[x];
  for (const auto& item : sketch.heavy_) reconstructed[item.first] += item.second.exact;
  assert(reconstructed == full_truth);
  assert(exact_rank(full_truth, INFINITY) == sketch.count_);
  std::vector<double> probes{-INFINITY, -10001, -1, 0,
      std::numeric_limits<double>::denorm_min(), 0.5, 1, 10001, INFINITY};
  probes.insert(probes.end(), extra_queries.begin(), extra_queries.end());
  for (std::size_t i = 0; i < sketch.nodes_.size(); ++i) {
    const auto x = sketch.nodes_[i].x;
    probes.push_back(x);
    probes.push_back(std::nextafter(x, -INFINITY));
    probes.push_back(std::nextafter(x, INFINITY));
    if (i) probes.push_back(S::midpoint(sketch.nodes_[i - 1].x, x));
  }
  std::cout << "{\"kind\":\"rank_queries\",\"additions\":"
            << sketch.buffer_.size() + sketch.heavy_.size() << ",\"queries\":[";
  bool comma = false;
  for (double x : probes) {
    if (comma) std::cout << ',';
    comma = true;
    std::cout << '['; dyadic(sketch.raw_rank(x));
    std::cout << ','; dyadic(sketch.rank(x));
    std::cout << ',' << exact_rank(full_truth, x) << ',' << bits(x) << ',';
    dyadic(S::spline_rank(sketch.nodes_, x));
    std::uint64_t pending = 0;
    for (double value : sketch.buffer_) if (value <= x) ++pending;
    for (const auto& item : sketch.heavy_) if (item.first <= x) pending += item.second.exact;
    const auto interval = sketch.rank_bounds(x);
    std::cout << ',' << pending << ',' << interval.lower << ',' << interval.upper;
    std::cout << ']';
  }
  std::cout << "]}\n";
}

void splinesketch::PaperCertificateInspector::operation_cases() {
  using S = CertifiedPaperSplineSketch;
  operation_mode = true;
  std::cout << "{\"kind\":\"format\",\"precision\":"
            << std::numeric_limits<long double>::digits << ",\"min_quantum\":"
            << std::numeric_limits<long double>::min_exponent -
                std::numeric_limits<long double>::digits << "}\n";
  struct State {
    S sketch;
    std::uint64_t q = 0, id = next_state++;
    std::map<double, std::uint64_t> absorbed = {}, truth = {};
  };
  const auto mutate = [](State& state, const auto& action) {
    rounding_q = state.q;
    active_state = state.id; absorbed_truth = state.absorbed; full_truth = state.truth;
    std::cout << "{\"kind\":\"rank_context\",\"parents\":[" << state.id << "]}\n";
    action(state.sketch);
    state.q = rounding_q;
    state.absorbed = absorbed_truth;
    ledger_state(state.sketch.nodes_, state.sketch.count_);
    completed(state.sketch);
  };
  const auto add = [&](State& state, double value) {
    ++state.truth[value];
    mutate(state, [&](S& s) { s.add(value); });
  };
  const auto merge = [](State& left, const State& right) {
    rounding_q = std::max(left.q, right.q);
    active_state = left.id; absorbed_truth = left.absorbed; full_truth = left.truth;
    for (const auto& item : right.absorbed) absorbed_truth[item.first] += item.second;
    for (const auto& item : right.truth) full_truth[item.first] += item.second;
    std::cout << "{\"kind\":\"rank_context\",\"parents\":[" << left.id << ',' << right.id << "]}\n";
    left.sketch.merge(right.sketch);
    left.q = rounding_q;
    left.absorbed = absorbed_truth; left.truth = full_truth;
    ledger_state(left.sketch.nodes_, left.sketch.count_);
    completed(left.sketch);
  };
  State large{S(6)}, nearby{S(6)};
  for (int i = 0; i < 6; ++i) add(large, i);
  mutate(large, [](S& s) { s.consolidate(); });
  for (unsigned i = 0; i < 52; ++i) merge(large, large);
  double x = 0x1.968087cfc2a08p+0;
  for (unsigned i = 0; i < 6; ++i) {
    add(nearby, x);
    x = std::nextafter(x, INFINITY);
  }
  mutate(nearby, [](S& s) { s.consolidate(); });
  merge(large, nearby); // Public binary64 decreasing-sample counterexample.
  extra_queries.push_back(paper_test::resize_history_query);
  for (auto policy : {S::BoundPolicy::practical, S::BoundPolicy::theoretical}) {
    State oscillating{S(16, policy)};
    for (double value : paper_test::resize_history_input) add(oscillating, value);
    mutate(oscillating, [](S& s) { s.consolidate(); });
    std::cout << "{\"kind\":\"resize_history\",\"begin\":true,\"query\":"
              << bits(paper_test::resize_history_query) << ",\"source_bound\":";
    dyadic(oscillating.sketch.max_rank_error());
    std::cout << ",\"source_rank\":"; dyadic(oscillating.sketch.rank(paper_test::resize_history_query));
    std::cout << ",\"truth\":" << exact_rank(full_truth, paper_test::resize_history_query);
    std::cout << "}\n";
    for (unsigned i = 0; i < 24; ++i)
      mutate(oscillating, [&](S& s) { s.resize(i % 2 ? 16 : 32); });
    std::cout << "{\"kind\":\"resize_history\",\"begin\":false}\n";
  }
  extra_queries.clear();
  for (unsigned seed : {1U, 8U, 19U}) {
    State state{S(7)};
    std::mt19937_64 random(seed);
    for (unsigned i = 0; i < 512; ++i) {
      const double value = i % 13 ? std::ldexp(static_cast<double>(random() % 31), -20) :
          (i % 2 ? -1.0 : 1.0) * (i + 1);
      add(state, value);
      if (i % 7 == 0) mutate(state, [](S& s) { s.consolidate(); });
      if (i % 11 == 0) {
        const auto at = (i / 11) % 20;
        mutate(state, [&](S& s) { s.resize(6 + (at <= 10 ? at : 20 - at)); });
      }
      if (i % 47 == 0) {
        State other{S(6 + (i % 17))};
        for (double value : {static_cast<double>(i + 10000), -static_cast<double>(i + 10000), 0.0})
          add(other, value);
        merge(state, other);
      }
    }
    mutate(state, [](S& s) { s.consolidate(); });
    // These operations add arithmetic work without adding observations.
    for (unsigned i = 0; i < 32; ++i)
      mutate(state, [&](S& s) { s.resize(i % 2 ? 7 : 8); });
    mutate(state, [](S& s) { s.finalize(); });
  }
}

void splinesketch::PaperCertificateInspector::initialization_cases() {
  using S = CertifiedPaperSplineSketch;
  quiet = true; // These are isolated states satisfying the local premises.
  initialization_mode = true;
  std::cout << "{\"kind\":\"format\",\"precision\":"
            << std::numeric_limits<long double>::digits << "}\n";
  const auto emit = [](std::size_t capacity, std::vector<S::Item> items, std::uint64_t count) {
    S sketch(capacity);
    sketch.count_ = count;
    S::Snapshot initial{{}, items, {}};
    initial.prepare(); sketch.initialize(initial);
    std::cout << "{\"kind\":\"initialization\",\"capacity\":" << capacity
              << ",\"n\":" << count << ",\"items\":[";
    for (std::size_t i = 0; i < items.size(); ++i) {
      if (i) std::cout << ',';
      std::cout << '[' << bits(items[i].first) << ',' << items[i].second << ']';
    }
    std::cout << "],\"nodes\":[";
    for (std::size_t i = 0; i < sketch.nodes_.size(); ++i) {
      if (i) std::cout << ',';
      const auto& node = sketch.nodes_[i];
      std::cout << '[' << bits(node.x) << ','; dyadic(node.prefix);
      std::cout << ',' << node.lower << ',' << node.upper << ',' << node.atom << ']';
    }
    std::cout << "]}\n";
  };
  std::mt19937_64 random(0x494e49544d415353ULL);
  for (std::size_t capacity : {6, 7, 32, 128}) {
    for (auto size : {capacity + 2, 2 * capacity + 1}) {
      for (auto limit : {std::uint64_t(17 * size), (std::uint64_t{1} << 53) - 1,
                         (std::uint64_t{1} << 53) + 1, (std::uint64_t{1} << 63) - 1}) {
        const auto weight = (limit - 1) / (size - 1);
        std::vector<S::Item> items;
        for (std::size_t i = 0; i + 1 < size; ++i) items.emplace_back(i, weight);
        items.emplace_back(size - 1, limit - (size - 1) * weight);
        emit(capacity, items, limit);
      }
      for (unsigned trial = 0; trial < 32; ++trial) {
        const auto weight = 1 + random() % ((std::uint64_t{1} << 62) / (2 * size));
        std::vector<S::Item> items;
        std::uint64_t count = 0;
        for (std::size_t i = 0; i < size; ++i) {
          const auto frequency = weight + random() % (weight / (2 * capacity) + 1);
          items.emplace_back(i, frequency); count += frequency;
        }
        emit(capacity, items, count);
      }
    }
    // Sparse roots exercise synthetic gap filling, with enough external count
    // that each absorbed frequency satisfies the MG release budget n/k.
    for (auto weight : {std::uint64_t{1}, (std::uint64_t{1} << 52) + 1})
      emit(capacity, {{0, weight}, {1, weight}, {4, weight}}, capacity * weight);
  }
}

void splinesketch::PaperCertificateInspector::selected(std::size_t index, std::uint64_t offset,
                                                     std::uint64_t total, std::size_t capacity) {
  if (initialization_mode)
    std::cout << "{\"kind\":\"initial_offset\",\"index\":" << index
              << ",\"offset\":" << offset << ",\"total\":" << total
              << ",\"capacity\":" << capacity << "}\n";
}

int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--operation-cases") {
    splinesketch::PaperCertificateInspector::operation_cases();
    return 0;
  }
  if (argc == 2 && std::string(argv[1]) == "--initialization-cases") {
    splinesketch::PaperCertificateInspector::initialization_cases();
    return 0;
  }
  unsigned count = 15360, seed = 8;
  std::size_t capacity = 512;
  bool theoretical = false;
  bool consolidate_each = false;
  std::string input, dump_path, pattern = "contracting";
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option == "--quiet") quiet = true;
    else if (option == "--geometry-only") geometry_only = true;
    else if (option == "--theoretical") theoretical = true;
    else if (option == "--consolidate-each") consolidate_each = true;
    else if (option == "--input" && i + 1 < argc) input = argv[++i];
    else if (option == "--dump-input" && i + 1 < argc) dump_path = argv[++i];
    else if (option == "--pattern" && i + 1 < argc) pattern = argv[++i];
    else if (option == "--count" && i + 1 < argc) count = static_cast<unsigned>(std::stoul(argv[++i]));
    else if (option == "--seed" && i + 1 < argc) seed = static_cast<unsigned>(std::stoul(argv[++i]));
    else if (option == "--capacity" && i + 1 < argc) capacity = std::stoull(argv[++i]);
    else throw std::invalid_argument("unknown local-charge test option");
  }
  using Sketch = splinesketch::CertifiedPaperSplineSketch;
  Sketch sketch(capacity, theoretical ? Sketch::BoundPolicy::theoretical : Sketch::BoundPolicy::practical);
  using Raw = splinesketch::PaperSplineSketch;
  Raw raw(capacity, theoretical ? Raw::BoundPolicy::theoretical : Raw::BoundPolicy::practical);
  std::ofstream dump;
  if (!dump_path.empty()) {
    dump.open(dump_path);
    if (!dump) throw std::runtime_error("cannot open dump output");
  }
  std::vector<double> values;
  if (input.empty()) {
    std::mt19937_64 random(seed);
    for (unsigned i = 0; i < count; ++i) {
      if (pattern == "grid") values.push_back(static_cast<double>(random() % (2 * capacity + 1)));
      else if (pattern == "staged-grid") {
        const auto shift = (i / (5 * capacity)) % 8;
        values.push_back(static_cast<double>(2 * capacity - ((random() % (2 * capacity + 1)) >> shift)));
      }
      else if (pattern == "contracting") {
        const double u = static_cast<double>(random() % 1000001) / 1000001;
        values.push_back(1 - std::ldexp(u, -static_cast<int>((i / 2560) % 40)));
      } else throw std::invalid_argument("unknown stream pattern");
    }
  } else {
    std::ifstream file(input);
    if (!file) throw std::runtime_error("cannot open input");
    std::uint64_t value;
    while (file >> value) {
      double x; std::memcpy(&x, &value, sizeof x); values.push_back(x);
    }
  }
  if (!quiet) {
    std::cout << "{\"kind\":\"format\",\"precision\":"
              << std::numeric_limits<long double>::digits << "}\n";
    std::vector<double> support = values;
    std::sort(support.begin(), support.end());
    support.erase(std::unique(support.begin(), support.end()), support.end());
    std::cout << "{\"kind\":\"support\",\"values\":[";
    for (std::size_t i = 0; i < support.size(); ++i) {
      if (i) std::cout << ',';
      std::cout << bits(support[i]);
    }
    std::cout << "]}\n";
  }
  for (double x : values) {
    if (dump.is_open()) dump << bits(x) << '\n';
    sketch.add(x);
    raw.add(x);
    if (consolidate_each) { sketch.consolidate(); raw.consolidate(); }
  }
  sketch.consolidate();
  raw.consolidate();
  if (!quiet && !geometry_only) {
    std::map<double, std::uint64_t> truth;
    for (double x : values) ++truth[x];
    std::vector<double> support;
    std::cout << "{\"kind\":\"query_checks\",\"truth\":[";
    for (const auto& item : truth) {
      if (!support.empty()) std::cout << ',';
      support.push_back(item.first);
      std::cout << '[' << bits(item.first) << ',' << item.second << ']';
    }
    std::cout << "],\"queries\":[";
    std::vector<double> queries{-INFINITY, INFINITY};
    const auto samples = std::min<std::size_t>(128, support.size());
    for (std::size_t i = 0; i < samples; ++i) {
      const double x = support[i * (support.size() - 1) / std::max<std::size_t>(1, samples - 1)];
      queries.push_back(x);
      queries.push_back(std::nextafter(x, -INFINITY));
      queries.push_back(std::nextafter(x, INFINITY));
    }
    std::sort(queries.begin(), queries.end());
    queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
    for (std::size_t i = 0; i < queries.size(); ++i) {
      if (i) std::cout << ',';
      const auto estimate = raw.rank(queries[i]);
      assert(std::isfinite(estimate));
      std::cout << '[' << bits(queries[i]) << ',' << bits(estimate) << ']';
    }
    std::cout << "]}\n";
  }
  std::cout << "{\"kind\":\"summary\",\"violations\":" << violations
            << ",\"passes\":" << passes << ",\"count\":" << sketch.count()
            << ",\"width\":" << sketch.max_rank_uncertainty() << "}\n";
}
