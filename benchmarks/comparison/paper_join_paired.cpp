// Replay identical reference-generated streams and query both node sets.
#define SPLINESKETCH_TESTING
#include <splinesketch/splinesketch.hpp>
#define splinesketch previous_splinesketch
#include SPLINESKETCH_PREVIOUS_CORE
#include SPLINESKETCH_PREVIOUS_PAPER
#undef splinesketch
#include <splinesketch/paper_splinesketch.hpp>
#include "../accuracy_workload.hpp"
#include "../../tests/data/paper_width_counterexample.hpp"
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

using Sketch = splinesketch::CertifiedPaperSplineSketch;
static std::uint64_t bits(double value) { std::uint64_t result; std::memcpy(&result, &value, sizeof result); return result; }
static const char* names[]{"ascending", "descending", "uniform", "moving_cluster", "alternating_clusters",
  "contracting", "expanding", "random_exponents", "distinct_cycles", "nested_left", "nested_right",
  "widest_attack", "mass_attack", "alternating_attack", "binade"};
static std::vector<double> witness() {
  std::vector<double> data;
  for (auto value : paper_width_counterexample) { double x; std::memcpy(&x, &value, sizeof x); data.push_back(x); }
  return data;
}
namespace splinesketch {
struct PaperCertificateInspector {
  static std::vector<double> queries(const Sketch& sketch) {
    std::vector<double> xs;
    for (std::size_t i = 0; i < sketch.nodes_.size(); ++i) {
      const auto x = sketch.nodes_[i].x;
      xs.push_back(x); xs.push_back(std::nextafter(x, -INFINITY)); xs.push_back(std::nextafter(x, INFINITY));
      if (i) xs.push_back(Sketch::midpoint(sketch.nodes_[i - 1].x, x));
    }
    return xs;
  }
};
}
using Inspector = splinesketch::PaperCertificateInspector;

namespace previous_splinesketch {
struct PaperCertificateInspector {
  using Sketch = CertifiedPaperSplineSketch;
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
    std::vector<double> xs;
    for (std::size_t i = 0; i < sketch.nodes_.size(); ++i) {
      const auto x = sketch.nodes_[i].x;
      xs.push_back(x); xs.push_back(std::nextafter(x, -INFINITY)); xs.push_back(std::nextafter(x, INFINITY));
      if (i) xs.push_back(Sketch::midpoint(sketch.nodes_[i - 1].x, x));
    }
    return xs;
  }
};
}
using Before = previous_splinesketch::CertifiedPaperSplineSketch;
using BeforeInspector = previous_splinesketch::PaperCertificateInspector;
template<class S> static S replay_pair(const std::vector<double>& data, unsigned k) {
  S sketch(k); for (double x : data) sketch.add(x); sketch.consolidate(); return sketch;
}
template<class S> static double checked_error(const S& sketch, double x, std::uint64_t truth, double bound) {
  const auto result = sketch.rank_with_error(x);
  const double error = std::abs(result.estimate - static_cast<double>(truth));
  if (truth < result.lower_rank || truth > result.upper_rank || error > result.max_error || error > bound)
    throw std::runtime_error("invalid paired certificate");
  return error;
}
static void emit_pair(const std::string& name, unsigned step, const Before& before, const Sketch& after,
                      std::vector<double> data) {
  std::uint64_t input_hash = 0xcbf29ce484222325ULL;
  for (double x : data) input_hash = (input_hash ^ bits(x)) * 0x100000001b3ULL;
  std::sort(data.begin(), data.end());
  auto xs = accuracy_workload::make_queries(data);
  const auto a = BeforeInspector::queries(before), b = Inspector::queries(after);
  xs.insert(xs.end(), a.begin(), a.end()); xs.insert(xs.end(), b.begin(), b.end());
  std::sort(xs.begin(), xs.end()); xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  xs.push_back(-INFINITY); xs.push_back(INFINITY);
  const auto eb = before.max_rank_error(), ea = after.max_rank_error();
  double worst_before = 0, worst_after = 0;
  for (double x : xs) {
    const auto truth = static_cast<std::uint64_t>(std::upper_bound(data.begin(), data.end(), x) - data.begin());
    worst_before = std::max(worst_before, checked_error(before, x, truth, eb));
    worst_after = std::max(worst_after, checked_error(after, x, truth, ea));
  }
  std::cout << name << ',' << step << ',' << before.bucket_capacity() << ',' << data.size() << ',' << input_hash
            << ',' << xs.size() << ',' << before.max_rank_uncertainty() << ',' << after.max_rank_uncertainty()
            << ',' << eb << ',' << ea << ',' << 100 * worst_before / data.size()
            << ',' << 100 * worst_after / data.size() << '\n' << std::flush;
}
static std::vector<double> shared_stream(unsigned k, unsigned pattern, unsigned seed, unsigned count,
                                        unsigned& before_peak, unsigned& after_peak) {
  Before before(k); Sketch after(k); std::mt19937_64 random(seed);
  std::vector<double> data; data.reserve(count); std::pair<double, double> interval{-1, 1};
  double max_before = -1, max_after = -1;
  for (unsigned i = 0; i < count; ++i) {
    if (i % (5 * k) == 0) interval = BeforeInspector::target(before, pattern == 12);
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
      default: x = std::ldexp(1.0, -1021) + (static_cast<double>(random() % 2001) - 1000) * std::numeric_limits<double>::denorm_min();
    }
    data.push_back(x); before.add(x); after.add(x);
    if ((i + 1) % (5 * k) == 0 && i + 1 >= 20 * k) {
      const double rb = static_cast<double>(before.max_rank_uncertainty()) * k / (i + 1);
      const double ra = static_cast<double>(after.max_rank_uncertainty()) * k / (i + 1);
      if (rb > max_before) { max_before = rb; before_peak = i + 1; }
      if (ra > max_after) { max_after = ra; after_peak = i + 1; }
    }
  }
  before.consolidate(); after.consolidate();
  if (static_cast<double>(before.max_rank_uncertainty()) * k / count > max_before) before_peak = count;
  if (static_cast<double>(after.max_rank_uncertainty()) * k / count > max_after) after_peak = count;
  return data;
}
static void compare_stream(unsigned k, unsigned pattern, unsigned seed, unsigned count) {
  unsigned bp = 0, ap = 0;
  const auto data = shared_stream(k, pattern, seed, count, bp, ap);
  std::vector<unsigned> lengths{bp, ap, count};
  std::sort(lengths.begin(), lengths.end()); lengths.erase(std::unique(lengths.begin(), lengths.end()), lengths.end());
  for (auto length : lengths) {
    const std::vector<double> prefix(data.begin(), data.begin() + length);
    emit_pair(std::string(names[pattern]) + "_seed" + std::to_string(seed), length,
              replay_pair<Before>(prefix, k), replay_pair<Sketch>(prefix, k), prefix);
  }
}
int main(int argc, char** argv) {
  const bool long_run = argc == 2 && std::string(argv[1]) == "--long";
  if (argc > 2 || (argc == 2 && !long_run)) throw std::invalid_argument("usage: paired [--long]");
  std::cout << std::setprecision(17)
            << "case,step,k,n,input_hash,queries,before_width,after_width,before_uniform,after_uniform,before_error_percent,after_error_percent\n";
  if (long_run) {
    for (unsigned seed : {97U, 103U}) compare_stream(512, 9, seed, 2000000);
    return 0;
  }
  for (unsigned k : {32U, 128U, 512U})
    for (unsigned pattern = 0; pattern < 15; ++pattern) compare_stream(k, pattern, 103, 60000);
  const auto data = witness();
  auto before = replay_pair<Before>(data, 32); auto after = replay_pair<Sketch>(data, 32);
  emit_pair("witness", 0, before, after, data);
  before.resize(1024); after.resize(1024); emit_pair("grow", 0, before, after, data);
  for (unsigned i = 1; i <= 6; ++i) {
    before.resize(6); after.resize(6); before.resize(1024); after.resize(1024);
    emit_pair("resize", i, before, after, data);
  }
  for (bool balanced : {true, false}) {
    std::vector<Before> bs; std::vector<Sketch> as;
    for (unsigned i = 0; i < 16; ++i) { bs.push_back(replay_pair<Before>(data, 32)); as.push_back(replay_pair<Sketch>(data, 32)); }
    if (balanced) for (unsigned gap = 1; gap < 16; gap *= 2)
      for (unsigned i = 0; i < 16; i += gap * 2) { bs[i].merge(bs[i + gap]); as[i].merge(as[i + gap]); }
    else for (unsigned i = 1; i < 16; ++i) { bs[0].merge(bs[i]); as[0].merge(as[i]); }
    std::vector<double> all; for (unsigned i = 0; i < 16; ++i) all.insert(all.end(), data.begin(), data.end());
    emit_pair(balanced ? "balanced_merge16" : "chain_merge16", 0, bs[0], as[0], all);
  }
}
