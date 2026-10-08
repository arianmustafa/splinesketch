// Adversarial fixed-capacity search, deletion reduction and history checks.
#define SPLINESKETCH_TESTING
#include <splinesketch/paper_splinesketch.hpp>
#include "../accuracy_workload.hpp"
#include "../../tests/data/paper_width_counterexample.hpp"
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
using Sketch = splinesketch::CertifiedPaperSplineSketch;
namespace splinesketch {
struct PaperCertificateInspector {
  static double factor(const Sketch& sketch) { return static_cast<double>(sketch.factor_); }
  static std::pair<double, double> target(const Sketch& sketch, bool mass) {
    std::pair<double, double> result{-1, 1}; long double best = -1;
    for (std::size_t i = 1; i < sketch.nodes_.size(); ++i) {
      const auto& left = sketch.nodes_[i - 1]; const auto& right = sketch.nodes_[i];
      if (!(std::nextafter(left.x, INFINITY) < right.x)) continue;
      const auto score = mass ? right.mass : static_cast<long double>(right.upper - right.atom - left.lower);
      if (score > best) { best = score; result = {left.x, right.x}; }
    }
    return result;
  }
  static std::vector<double> queries(const Sketch& sketch) {
    std::vector<double> result;
    for (std::size_t i = 0; i < sketch.nodes_.size(); ++i) {
      const double x = sketch.nodes_[i].x;
      result.push_back(x); result.push_back(std::nextafter(x, -INFINITY)); result.push_back(std::nextafter(x, INFINITY));
      if (i) result.push_back(Sketch::midpoint(sketch.nodes_[i - 1].x, x));
    }
    return result;
  }
};
}
using Inspector = splinesketch::PaperCertificateInspector;
struct Metrics { double error = 0; std::uint64_t hash = 0xcbf29ce484222325ULL; std::size_t queries = 0; };
static std::uint64_t bits(double value) { std::uint64_t result; std::memcpy(&result, &value, sizeof result); return result; }
static Metrics measure(const Sketch& sketch, std::vector<double> data) {
  std::sort(data.begin(), data.end());
  auto xs = accuracy_workload::make_queries(data); auto knots = Inspector::queries(sketch);
  xs.insert(xs.end(), knots.begin(), knots.end());
  const auto uniform = sketch.max_rank_error(); Metrics result;
  for (double x : xs) {
    const auto truth = static_cast<std::uint64_t>(std::upper_bound(data.begin(), data.end(), x) - data.begin());
    const auto estimate = sketch.rank_with_error(x); const double error = std::abs(estimate.estimate - static_cast<double>(truth));
    if (estimate.lower_rank > truth || estimate.upper_rank < truth || error > estimate.max_error || error > uniform)
      throw std::runtime_error("invalid certificate in width search");
    result.error = std::max(result.error, error);
    for (auto field : {bits(x), bits(estimate.estimate), estimate.lower_rank, estimate.upper_rank, bits(estimate.max_error)})
      result.hash = (result.hash ^ field) * 0x100000001b3ULL;
    ++result.queries;
  }
  return result;
}
static Sketch replay(const std::vector<double>& data, std::size_t k) {
  Sketch sketch(k); for (double x : data) sketch.add(x); sketch.consolidate(); return sketch;
}
static const char* names[]{"ascending", "descending", "uniform", "moving_cluster", "alternating_clusters",
  "contracting", "expanding", "random_exponents", "distinct_cycles", "nested_left", "nested_right",
  "widest_attack", "mass_attack", "alternating_attack", "binade"};
static std::vector<double> generate(unsigned k, unsigned pattern, unsigned seed, unsigned count,
                                    unsigned& peak_n, std::uint64_t& peak_width, double& peak_factor) {
  Sketch sketch(k); std::mt19937_64 random(seed); std::vector<double> data;
  std::pair<double, double> interval{-1, 1}; double best = -1;
  for (unsigned i = 0; i < count; ++i) {
    if (i % (5 * k) == 0) interval = Inspector::target(sketch, pattern == 12);
    const auto block = i / (5 * k); const double u = static_cast<double>(random() % 1000001) / 1000001;
    double x = 0;
    switch (pattern) {
      case 0: x = i; break;
      case 1: x = -static_cast<double>(i); break;
      case 2: x = u; break;
      case 3: x = block + u; break;
      case 4: x = (block % 2 ? 1 : -1) * (1 + u * 0.001); break;
      case 5: x = std::ldexp(1 + u, -static_cast<int>(block % 1000)); break;
      case 6: x = std::ldexp(1 + u, static_cast<int>(block % 1000)); break;
      case 7: x = std::ldexp(1 + u, static_cast<int>(random() % 2001) - 1000); break;
      case 8: x = i % (k + 1); break;
      case 9: x = std::ldexp(u, -static_cast<int>(block % 1000)); break;
      case 10: x = 1 - std::ldexp(u, -static_cast<int>(block % 50)); break;
      case 11: case 12: case 13: {
        const double t = pattern == 13 && block % 2 ? 0.999 - u * 0.0001 : 0.0001 + u * 0.0001;
        x = interval.first + (interval.second - interval.first) * t;
        if (x == interval.first) x = std::nextafter(x, interval.second);
        break;
      }
      default: x = std::ldexp(1.0, -1021) + (static_cast<double>(random() % 2001) - 1000) * std::numeric_limits<double>::denorm_min(); break;
    }
    sketch.add(x); data.push_back(x);
    if ((i + 1) % (5 * k) == 0 && i + 1 >= 20 * k) {
      const auto width = sketch.max_rank_uncertainty(); const double ratio = static_cast<double>(width) * k / (i + 1);
      if (ratio > best) { best = ratio; peak_n = i + 1; peak_width = width; peak_factor = Inspector::factor(sketch); }
    }
  }
  sketch.consolidate(); const auto width = sketch.max_rank_uncertainty();
  if (static_cast<double>(width) * k / count > best) { peak_n = count; peak_width = width; peak_factor = Inspector::factor(sketch); }
  return data;
}
static std::vector<double> witness() {
  std::vector<double> data;
  for (auto value : paper_width_counterexample) { double x; std::memcpy(&x, &value, sizeof x); data.push_back(x); }
  return data;
}
static bool hypothesis(const std::vector<double>& data) {
  if (data.empty()) return false;
  try { const auto sketch = replay(data, 32); return sketch.max_rank_uncertainty() * 32 >= 12 * data.size(); }
  catch (const std::logic_error&) { return false; }
}
static void reduce() {
  unsigned peak_n = 0; std::uint64_t peak_width = 0; double factor = 0;
  auto data = generate(32, 12, 71, 640, peak_n, peak_width, factor);
  std::size_t pieces = 2, evaluations = 0;
  while (data.size() > 1) {
    bool reduced = false; const auto chunk = (data.size() + pieces - 1) / pieces;
    for (std::size_t first = 0; first < data.size(); first += chunk) {
      auto candidate = data;
      candidate.erase(candidate.begin() + static_cast<std::ptrdiff_t>(first),
                      candidate.begin() + static_cast<std::ptrdiff_t>(std::min(data.size(), first + chunk)));
      ++evaluations;
      if (hypothesis(candidate)) { data = std::move(candidate); pieces = std::max<std::size_t>(2, pieces - 1); reduced = true; break; }
    }
    if (reduced) continue;
    if (pieces >= data.size()) break;
    pieces = std::min(data.size(), pieces * 2);
  }
  if (data != witness()) throw std::runtime_error("reduction differs from saved witness");
  std::cout << "Reduced 640 observations to " << data.size() << " in " << evaluations
            << " evaluations; every single deletion fails k*W/n >= 12\n";
}
static void histories() {
  std::cout << "history,step,k,n,width,uniform,error_percent,queries,query_hash\n";
  const auto base = witness(); auto sketch = replay(base, 32); std::vector<double> data = base;
  const auto emit = [&](const char* history, unsigned step, const Sketch& current, const std::vector<double>& values) {
    const auto metrics = measure(current, values);
    std::cout << history << ',' << step << ',' << current.bucket_capacity() << ',' << current.count() << ','
              << current.max_rank_uncertainty() << ',' << current.max_rank_error() << ','
              << 100 * metrics.error / current.count() << ',' << metrics.queries << ',' << metrics.hash << '\n';
  };
  emit("witness", 0, sketch, data); sketch.resize(1024); emit("grow", 0, sketch, data);
  for (unsigned i = 0; i < 6; ++i) { sketch.resize(6); sketch.resize(1024); emit("shrink_grow", i + 1, sketch, data); }
  for (bool balanced : {true, false}) {
    std::vector<Sketch> parts; for (unsigned i = 0; i < 16; ++i) parts.push_back(replay(base, 32));
    if (balanced) for (unsigned gap = 1; gap < 16; gap *= 2)
      for (unsigned i = 0; i < 16; i += gap * 2) parts[i].merge(parts[i + gap]);
    else for (unsigned i = 1; i < 16; ++i) parts[0].merge(parts[i]);
    data.clear(); for (unsigned i = 0; i < 16; ++i) data.insert(data.end(), base.begin(), base.end());
    emit(balanced ? "balanced_merge16" : "chain_merge16", 4, parts[0], data);
  }
}
int main(int argc, char** argv) {
  std::cout << std::setprecision(17);
  if (argc == 2 && std::string(argv[1]) == "--reduce") { reduce(); return 0; }
  if (argc == 2 && std::string(argv[1]) == "--histories") { histories(); return 0; }
  const bool long_run = argc >= 2 && std::string(argv[1]) == "--long";
  if (argc > 3 || (argc >= 2 && !long_run)) throw std::invalid_argument("usage: width-search [--long [SEED]|--reduce|--histories]");
  const unsigned seed = argc == 3 ? static_cast<unsigned>(std::stoul(argv[2])) : 71;
  const unsigned count = long_run ? 2000000 : 200000;
  std::cout << "k,pattern,seed,peak_n,peak_width,peak_uniform,peak_normalized,peak_factor,peak_error_percent,peak_queries,peak_hash,final_n,final_width,final_uniform,final_error_percent,final_queries,final_hash\n";
  for (unsigned k : {32U, 128U, 512U}) for (unsigned pattern = 0; pattern < 15; ++pattern) {
    unsigned peak_n = 0; std::uint64_t peak_width = 0; double factor = 0;
    const auto data = generate(k, pattern, seed, count, peak_n, peak_width, factor);
    auto peak_data = std::vector<double>(data.begin(), data.begin() + peak_n);
    const auto peak = replay(peak_data, k), final = replay(data, k);
    if (peak.max_rank_uncertainty() != peak_width) throw std::runtime_error("peak replay mismatch");
    const auto a = measure(peak, peak_data), b = measure(final, data);
    std::cout << k << ',' << names[pattern] << ',' << seed << ',' << peak_n << ',' << peak_width << ',' << peak.max_rank_error()
              << ',' << static_cast<double>(peak_width) * k / peak_n << ',' << factor << ',' << 100 * a.error / peak_n << ',' << a.queries << ',' << a.hash
              << ',' << count << ',' << final.max_rank_uncertainty() << ',' << final.max_rank_error() << ',' << 100 * b.error / count << ',' << b.queries << ',' << b.hash << '\n' << std::flush;
  }
}
