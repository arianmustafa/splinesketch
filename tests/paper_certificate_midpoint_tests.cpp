#include <splinesketch/paper_splinesketch.hpp>
#include "../benchmarks/comparison/experiments/paper_certificate_midpoint.hpp"
#include "paper_resize_history_fixture.hpp"
#include <cfenv>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using Sketch = splinesketch::CertifiedPaperSplineSketch;
using namespace splinesketch::experimental;

static void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

static std::uint64_t bits(double value) {
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof result);
  return result;
}

namespace splinesketch {
struct PaperCertificateInspector {
  static std::vector<double> partition(const Sketch& sketch) {
    std::vector<double> queries{-INFINITY, INFINITY};
    const auto append = [&](double key) {
      queries.push_back(key);
      queries.push_back(std::nextafter(key, -INFINITY));
      queries.push_back(std::nextafter(key, INFINITY));
    };
    for (std::size_t i = 0; i < sketch.nodes_.size(); ++i) {
      const auto& node = sketch.nodes_[i];
      require(node.lower <= node.upper && node.upper <= sketch.count_ && node.atom <= node.lower,
              "invalid certificate node");
      if (i) {
        const auto& previous = sketch.nodes_[i - 1];
        require(previous.x < node.x, "grid keys not ordered");
        require(previous.lower <= node.lower && node.atom <= node.lower - previous.lower &&
                previous.upper <= node.upper && node.atom <= node.upper - previous.upper,
                "envelope jump smaller than retained atom");
      }
      append(node.x);
    }
    for (double value : sketch.buffer_) append(value);
    for (const auto& item : sketch.heavy_) append(item.first);
    std::sort(queries.begin(), queries.end());
    queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
    // Both endpoints, and thus the midpoint, are step functions. Keys and
    // neighbors exhaust the representable cells of their common partition;
    // consecutive floating-point keys have no missing open cell to sample.
    return queries;
  }
};
} // namespace splinesketch

static std::size_t check_monotone_partition(const Sketch& sketch) {
  const auto queries = splinesketch::PaperCertificateInspector::partition(sketch);
  const auto uniform = certificate_midpoint_max_error(sketch);
  const auto top = certificate_midpoint(sketch.count(), sketch.count());
  std::uint64_t previous_lower = 0, previous_upper = 0;
  double previous_estimate = 0;
  for (double query : queries) {
    const auto r = certificate_midpoint_rank(sketch, query);
    require(previous_lower <= r.lower_rank && previous_upper <= r.upper_rank,
            "certificate endpoint decreased");
    require(previous_estimate <= r.estimate, "midpoint estimate decreased");
    require(r.upper_rank <= sketch.count() && r.estimate <= top.estimate,
            "midpoint exceeded terminal rank");
    require(r.max_error <= uniform, "uniform bound excludes a certificate cell");
    previous_lower = r.lower_rank; previous_upper = r.upper_rank; previous_estimate = r.estimate;
  }
  require(previous_lower == sketch.count() && previous_upper == sketch.count() &&
          previous_estimate == top.estimate, "terminal rank not exact before conversion");
  const auto bottom = certificate_midpoint_rank(sketch, -INFINITY);
  require(bottom.lower_rank == 0 && bottom.upper_rank == 0 && bottom.estimate == 0,
          "initial rank not zero");
  // IEEE comparisons identify both signs of zero with the same breakpoint.
  const auto negative_zero = certificate_midpoint_rank(sketch, -0.0);
  const auto positive_zero = certificate_midpoint_rank(sketch, 0.0);
  require(negative_zero.lower_rank == positive_zero.lower_rank &&
          negative_zero.upper_rank == positive_zero.upper_rank &&
          bits(negative_zero.estimate) == bits(positive_zero.estimate), "signed zero changed midpoint");
  return queries.size();
}

static void arithmetic_checks() {
  for (std::uint64_t lower = 0; lower <= 64; ++lower)
    for (std::uint64_t upper = lower; upper <= 64; ++upper) {
      const auto r = certificate_midpoint(lower, upper);
      require(r.estimate == (lower + upper) / 2.0 && r.max_error == (upper - lower) / 2.0,
              "small exact midpoint failed");
    }
  const auto odd = certificate_midpoint((std::uint64_t{1} << 53) + 1,
                                       (std::uint64_t{1} << 53) + 2);
  require(odd.estimate == 0x1.0000000000001p53 && odd.max_error == 1,
          "half bit lost by double rounding");
  const auto top = std::numeric_limits<std::uint64_t>::max();
  const auto singleton = certificate_midpoint(top, top);
  require(singleton.estimate == 0x1p64 && singleton.max_error == 1, "2^64 boundary failed");
  const auto full = certificate_midpoint(0, top);
  require(full.estimate == 0x1p63 && full.max_error == 0x1p63, "full integer range failed");
  bool rejected = false;
  try { certificate_midpoint(2, 1); } catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "reversed bounds accepted");
  rejected = false;
  try { certificate_midpoint_uniform_error(2, 1); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "impossible uniform width accepted");
  Sketch empty(8);
  rejected = false;
  try { certificate_midpoint_select(empty, 1); } catch (const std::logic_error&) { rejected = true; }
  require(rejected, "empty selection accepted");
  empty.add(7);
  for (auto target : {std::uint64_t{0}, std::uint64_t{2}}) {
    rejected = false;
    try { certificate_midpoint_select(empty, target); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "invalid selection target accepted");
  }
  for (unsigned bit = 0; bit < 53; ++bit) empty.merge(empty);
  empty.add(7); empty.consolidate();
  require(empty.count() == (std::uint64_t{1} << 53) + 1 &&
          certificate_midpoint_rank(empty, 7).estimate == 0x1p53 &&
          certificate_midpoint_select(empty, empty.count()) == 7,
          "downward rounded terminal target lost");
}

struct Comparison {
  unsigned better = 0, equal = 0, worse = 0;
  std::uint64_t queries = 0, partition_queries = 0, selections = 0;
};

static void check_state(const Sketch& sketch, std::vector<double> data, Comparison& comparison) {
  std::sort(data.begin(), data.end());
  require(sketch.count() == data.size(), "count changed by query adapter");
  comparison.partition_queries += check_monotone_partition(sketch);
  std::vector<double> queries{-INFINITY, INFINITY};
  for (std::size_t i = 0; i < data.size(); ++i) {
    queries.push_back(data[i]);
    queries.push_back(std::nextafter(data[i], -INFINITY));
    queries.push_back(std::nextafter(data[i], INFINITY));
    if (i) queries.push_back(data[i - 1] / 2 + data[i] / 2);
  }
  const auto partition = splinesketch::PaperCertificateInspector::partition(sketch);
  queries.insert(queries.end(), partition.begin(), partition.end());
  std::sort(queries.begin(), queries.end());
  queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
  const auto uniform = certificate_midpoint_max_error(sketch);
  // Independent inverse oracle: walk the complete step-function partition,
  // without bit-order conversion or binary search. All counts here are small.
  std::size_t expected_at = 0;
  double previous_selection = -INFINITY;
  for (std::uint64_t target = 1; target <= sketch.count(); ++target) {
    while (expected_at < partition.size()) {
      const auto bounds = sketch.rank_bounds(partition[expected_at]);
      if (std::isfinite(partition[expected_at]) && (bounds.lower + bounds.upper) / 2.0 >= target) break;
      ++expected_at;
    }
    require(expected_at < partition.size(), "partition oracle found no finite crossing");
    const auto selected = certificate_midpoint_select(sketch, target);
    require(selected == partition[expected_at] && std::isfinite(selected), "selection not first crossing");
    require(selected != 0 || !std::signbit(selected), "selection returned negative zero");
    require(previous_selection <= selected, "selection decreased as target increased");
    previous_selection = selected;
    const auto predecessor = std::nextafter(selected, -INFINITY);
    const auto before = sketch.rank_bounds(predecessor), at = sketch.rank_bounds(selected);
    require((at.lower + at.upper) / 2.0 >= target && (before.lower + before.upper) / 2.0 < target,
            "selection lacks exact midpoint crossing");
    const auto rank_before = static_cast<std::uint64_t>(std::upper_bound(data.begin(), data.end(), predecessor) - data.begin());
    const auto rank_at = static_cast<std::uint64_t>(std::upper_bound(data.begin(), data.end(), selected) - data.begin());
    const auto distance = target < rank_before ? rank_before - target :
                          target > rank_at ? target - rank_at : 0;
    const auto width = sketch.max_rank_uncertainty();
    require(2 * distance <= width && 2 * rank_at + width >= 2 * target &&
            2 * rank_before < 2 * target + width, "quantile rank-bracket guarantee failed");
    ++comparison.selections;
  }
  double projected_worst = 0, midpoint_worst = 0;
  for (double query : queries) {
    const auto truth = static_cast<std::uint64_t>(std::upper_bound(data.begin(), data.end(), query) - data.begin());
    const auto projected = sketch.rank_with_error(query);
    const auto r = certificate_midpoint_rank(sketch, query);
    require(r.lower_rank == projected.lower_rank && r.upper_rank == projected.upper_rank,
            "adapter changed integer certificate");
    require(r.lower_rank <= truth && truth <= r.upper_rank, "certificate excludes truth");
    const double error = std::fabs(r.estimate - static_cast<double>(truth));
    require(error <= r.max_error && error <= uniform, "midpoint allowance too small");
    const auto projected_radius = std::max(std::fabs(projected.estimate - r.lower_rank),
                                         std::fabs(projected.estimate - r.upper_rank));
    require(r.max_error <= projected_radius, "midpoint lost certificate minimax comparison");
    require(r.estimate == (r.lower_rank + r.upper_rank) / 2.0 &&
            r.max_error == (r.upper_rank - r.lower_rank) / 2.0,
            "small reachable certificate midpoint not exact");
    projected_worst = std::max(projected_worst, std::fabs(projected.estimate - truth));
    midpoint_worst = std::max(midpoint_worst, error);
  }
  if (midpoint_worst < projected_worst) ++comparison.better;
  else if (midpoint_worst > projected_worst) ++comparison.worse;
  else ++comparison.equal;
  comparison.queries += queries.size();
}

static void history_checks(bool theoretical, Comparison& comparison) {
  const auto policy = theoretical ? Sketch::BoundPolicy::theoretical : Sketch::BoundPolicy::practical;
  std::mt19937_64 random(223);
  for (unsigned shape = 0; shape < 7; ++shape) {
    Sketch sketch(16, policy);
    std::vector<double> data;
    check_state(sketch, data, comparison);
    for (unsigned i = 0; i < 149; ++i) {
      double value = 0;
      switch (shape) {
        case 0: value = static_cast<double>(random() % 17) - 8; break;
        case 1: value = std::ldexp(static_cast<double>(random() % 31), -20); break;
        case 2: value = static_cast<double>(i * i); break;
        case 3: value = -static_cast<double>(i * i); break;
        case 4: value = (static_cast<double>(random() % 17) - 8) *
                           std::numeric_limits<double>::denorm_min(); break;
        case 5: value = 1 + std::ldexp(static_cast<double>(random() % 17), -52); break;
        case 6: value = (i % 2 ? -1 : 1) * (i % 11 == 0 ? std::numeric_limits<double>::max() :
                           std::ldexp(1.0, static_cast<int>(random() % 2001) - 1000)); break;
      }
      sketch.add(value); data.push_back(value);
      if (i % 29 == 0) check_state(sketch, data, comparison); // Include nonempty exact buffers.
    }
    sketch.consolidate();
    check_state(sketch, data, comparison);
    for (unsigned step = 0; step < 8; ++step) {
      sketch.resize(step % 2 ? 16 : 8);
      check_state(sketch, data, comparison);
    }
    Sketch other(8, policy);
    for (unsigned i = 0; i < 47; ++i) { other.add(i / 4.0); data.push_back(i / 4.0); }
    sketch.merge(other);
    check_state(sketch, data, comparison);
    sketch.finalize();
    check_state(sketch, data, comparison);
  }
  Sketch cluster(16, policy);
  const std::vector<double> cluster_data(paper_test::resize_history_input.begin(),
                                       paper_test::resize_history_input.end());
  for (double x : cluster_data) cluster.add(x);
  cluster.consolidate();
  check_state(cluster, cluster_data, comparison);
  for (unsigned step = 0; step < 24; ++step) {
    cluster.resize(step % 2 ? 16 : 32);
    check_state(cluster, cluster_data, comparison);
  }

  // Both endpoints are reachable with identical active state (the separate
  // resize information test checks that premise). This makes radius 2.5 tight
  // at q=18, while also showing that a midpoint can hurt one actual stream.
  for (bool high : {false, true}) {
    Sketch witness(6, policy);
    std::vector<double> data;
    for (unsigned i = 0; i < 30; ++i) data.push_back(i * i);
    for (unsigned i = 1; i <= 5; ++i) data[i] = (high ? 0 : 18) + 3 * i;
    for (double value : data) witness.add(value);
    witness.consolidate();
    const auto r = certificate_midpoint_rank(witness, 18);
    const auto truth = high ? 6 : 1;
    require(r.lower_rank == 1 && r.upper_rank == 6 && r.estimate == 3.5 && r.max_error == 2.5 &&
            std::fabs(r.estimate - truth) == r.max_error, "tight reachable witness failed");
    if (high) require(std::fabs(witness.rank(18) - truth) < r.max_error,
                      "expected pointwise accuracy tradeoff missing");
    const auto selected = certificate_midpoint_select(witness, 3);
    require(selected == std::numeric_limits<double>::denorm_min() &&
            witness.max_rank_uncertainty() == 5 &&
            std::count_if(data.begin(), data.end(), [selected](double x) { return x <= selected; }) == 1,
            "tight integer rank-bracket witness failed");
    check_state(witness, data, comparison);
  }
  bool rejected = false;
  try { certificate_midpoint_rank(cluster, std::numeric_limits<double>::quiet_NaN()); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "NaN query accepted");
}

static void intervals() {
  std::uint64_t lower, upper;
  while (std::cin >> lower >> upper) {
    std::cout << "{\"lower\":" << lower << ",\"upper\":" << upper << ",\"modes\":[";
    bool first = true;
    for (int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
      require(std::fesetround(mode) == 0, "rounding mode unavailable");
      const auto r = certificate_midpoint(lower, upper);
      if (!first) std::cout << ',';
      first = false;
      std::cout << '[' << bits(r.estimate) << ',' << bits(r.max_error) << ','
                << bits(certificate_midpoint_uniform_error(upper - lower, upper)) << ']';
    }
    const auto midpoint = lower + (upper - lower) / 2;
    std::vector<std::uint64_t> targets{0, lower, upper, midpoint, std::numeric_limits<std::uint64_t>::max()};
    if (midpoint) targets.push_back(midpoint - 1);
    if (midpoint < std::numeric_limits<std::uint64_t>::max()) targets.push_back(midpoint + 1);
    std::sort(targets.begin(), targets.end());
    targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
    std::cout << "],\"targets\":[";
    first = true;
    for (auto target : targets) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << '[' << target << ',' << midpoint_detail::reaches(lower, upper, target) << ']';
    }
    std::cout << "]}\n";
  }
  require(std::cin.eof(), "invalid interval input");
  require(std::fesetround(FE_TONEAREST) == 0, "rounding mode restore failed");
}

static void large_counts(bool emit) {
  for (bool theoretical : {false, true}) {
    Sketch sketch(8, theoretical ? Sketch::BoundPolicy::theoretical : Sketch::BoundPolicy::practical);
    sketch.add(0);
    std::uint64_t zeros = 1, count = 1;
    for (unsigned step = 1; step <= 63; ++step) {
      // Public operations reach count 2^64-1 without materializing the stream.
      sketch.merge(sketch);
      const int value = step % 2;
      sketch.add(value);
      sketch.consolidate();
      zeros = 2 * zeros + (value == 0);
      count = 2 * count + 1;
      require(sketch.count() == count, "large public history count differs");
      check_monotone_partition(sketch);
      for (auto target : {std::uint64_t{1}, zeros, zeros + 1, count})
        require(certificate_midpoint_select(sketch, target) == (target <= zeros ? 0 : 1),
                "large exact weighted selection failed");
      if (step < 51) continue;
      for (int query : {-1, 0, 1}) {
        const auto truth = query < 0 ? 0 : query == 0 ? zeros : count;
        const auto r = certificate_midpoint_rank(sketch, query);
        require(r.lower_rank <= truth && truth <= r.upper_rank, "large certificate excludes truth");
        require(std::isfinite(r.estimate) && std::isfinite(r.max_error), "nonfinite large query");
        if (step == 63 && query == 1)
          require(r.estimate == 0x1p64 && r.max_error == 1, "public UINT64_MAX return allowance failed");
        if (emit)
          std::cout << "{\"step\":" << step << ",\"theoretical\":" << theoretical
                    << ",\"query\":" << query << ",\"count\":" << count << ",\"truth\":" << truth
                    << ",\"lower\":" << r.lower_rank << ",\"upper\":" << r.upper_rank
                    << ",\"estimate_bits\":" << bits(r.estimate) << ",\"radius_bits\":" << bits(r.max_error)
                    << ",\"uniform_bits\":" << bits(certificate_midpoint_max_error(sketch)) << "}\n";
      }
    }
  }
}

using WeightedData = std::vector<std::pair<double, std::uint64_t>>;

static void emit_selection_case(const char* name, const Sketch& sketch, const WeightedData& data,
                                const std::vector<std::uint64_t>& targets) {
  std::cout << "{\"case\":\"" << name << "\",\"count\":" << sketch.count()
            << ",\"width\":" << sketch.max_rank_uncertainty() << ",\"data\":[";
  bool first = true;
  for (auto [value, weight] : data) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '[' << bits(value) << ',' << weight << ']';
  }
  std::cout << "],\"profile\":[";
  first = true;
  for (auto query : splinesketch::PaperCertificateInspector::partition(sketch)) {
    if (!first) std::cout << ',';
    first = false;
    const auto bounds = sketch.rank_bounds(query);
    std::cout << '[' << bits(query) << ',' << bounds.lower << ',' << bounds.upper << ']';
  }
  std::cout << "],\"selections\":[";
  first = true;
  for (auto target : targets) {
    const auto value = certificate_midpoint_select(sketch, target);
    const auto at = sketch.rank_bounds(value), before = sketch.rank_bounds(std::nextafter(value, -INFINITY));
    if (!first) std::cout << ',';
    first = false;
    std::cout << "{\"target\":" << target << ",\"value_bits\":" << bits(value)
              << ",\"lower\":" << at.lower << ",\"upper\":" << at.upper
              << ",\"before_lower\":" << before.lower << ",\"before_upper\":" << before.upper
              << ",\"modes\":[";
    bool first_mode = true;
    for (int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
      require(std::fesetround(mode) == 0, "rounding mode unavailable");
      const auto selected = certificate_midpoint_select(sketch, target);
      require(bits(value) == bits(selected), "selection depends on rounding mode");
      if (!first_mode) std::cout << ',';
      first_mode = false;
      std::cout << bits(selected);
    }
    require(std::fesetround(FE_TONEAREST) == 0, "rounding mode restore failed");
    std::cout << "]}";
  }
  std::cout << "]}\n";
}

static void selection_fixtures() {
  for (bool theoretical : {false, true}) {
    const auto policy = theoretical ? Sketch::BoundPolicy::theoretical : Sketch::BoundPolicy::practical;
    Sketch cluster(16, policy);
    WeightedData cluster_data;
    std::vector<std::uint64_t> targets;
    for (double value : paper_test::resize_history_input) { cluster.add(value); cluster_data.emplace_back(value, 1); }
    cluster.consolidate();
    for (unsigned step = 0; step < 24; ++step) cluster.resize(step % 2 ? 16 : 32);
    for (std::uint64_t target = 1; target <= cluster.count(); ++target) targets.push_back(target);
    emit_selection_case(theoretical ? "cluster-theoretical" : "cluster-practical", cluster, cluster_data, targets);

    Sketch tight(6, policy);
    WeightedData tight_data;
    for (unsigned i = 0; i < 30; ++i) {
      const double value = i >= 1 && i <= 5 ? 3 * i : i * i;
      tight.add(value); tight_data.emplace_back(value, 1);
    }
    tight.consolidate();
    emit_selection_case(theoretical ? "tight-theoretical" : "tight-practical", tight, tight_data, {3, 4});

    Sketch duplicates(8, policy);
    for (unsigned i = 0; i < 99; ++i) duplicates.add(7);
    duplicates.consolidate();
    require(duplicates.max_rank_uncertainty() == 0 && certificate_midpoint_select(duplicates, 1) == 7,
            "duplicate atom selection failed");
    emit_selection_case(theoretical ? "duplicates-theoretical" : "duplicates-practical", duplicates, {{7, 99}}, {1, 50, 99});

    Sketch terminal(8, policy);
    terminal.add(7);
    for (unsigned bit = 0; bit < 53; ++bit) terminal.merge(terminal);
    terminal.add(7); terminal.consolidate();
    require(terminal.count() == (std::uint64_t{1} << 53) + 1 &&
            certificate_midpoint_rank(terminal, 7).estimate == 0x1p53 &&
            certificate_midpoint_select(terminal, terminal.count()) == 7,
            "downward rounded terminal target lost");
    emit_selection_case(theoretical ? "rounded-terminal-theoretical" : "rounded-terminal-practical",
                        terminal, {{7, terminal.count()}}, {1, std::uint64_t{1} << 53, terminal.count()});

    Sketch raw(8, policy);
    WeightedData raw_data;
    const auto tiny = std::numeric_limits<double>::denorm_min();
    for (double value : {-std::numeric_limits<double>::max(), -tiny, -0.0, 0.0, tiny,
                         std::numeric_limits<double>::max()}) { raw.add(value); raw_data.emplace_back(value, 1); }
    emit_selection_case(theoretical ? "raw-theoretical" : "raw-practical", raw, raw_data, {1, 2, 3, 4, 5, 6});
    raw.finalize();
    emit_selection_case(theoretical ? "finalized-theoretical" : "finalized-practical", raw, raw_data, {1, 2, 3, 4, 5, 6});

    Sketch large(8, policy); large.add(0);
    std::uint64_t zeros = 1;
    for (unsigned step = 1; step <= 63; ++step) {
      large.merge(large); large.add(step % 2); large.consolidate();
      zeros = 2 * zeros + (step % 2 == 0);
      if (step != 53 && step != 63) continue;
      targets = {1, zeros - 1, zeros, zeros + 1, large.count() - 1, large.count()};
      emit_selection_case(step == 53 ? (theoretical ? "large-53-theoretical" : "large-53-practical") :
                                     (theoretical ? "large-63-theoretical" : "large-63-practical"),
                          large, {{0, zeros}, {1, large.count() - zeros}}, targets);
    }
  }
}

int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--intervals") { intervals(); return 0; }
  if (argc == 2 && std::string(argv[1]) == "--large-counts") { large_counts(true); return 0; }
  if (argc == 2 && std::string(argv[1]) == "--selection-fixtures") { selection_fixtures(); return 0; }
  require(argc == 1, "usage: midpoint-tests [--intervals|--large-counts|--selection-fixtures]");
  arithmetic_checks();
  large_counts(false);
  Comparison comparison;
  for (bool theoretical : {false, true}) history_checks(theoretical, comparison);
  require(comparison.queries > 10000, "missing history coverage");
  std::cout << "{\"queries\":" << comparison.queries << ",\"partition_queries\":" << comparison.partition_queries
            << ",\"selections\":" << comparison.selections
            << ",\"checkpoints_better\":" << comparison.better
            << ",\"checkpoints_equal\":" << comparison.equal << ",\"checkpoints_worse\":"
            << comparison.worse << "}\n";
}
