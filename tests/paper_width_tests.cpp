#include <splinesketch/paper_splinesketch.hpp>
#include "data/paper_width_counterexample.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <map>
#include <random>

namespace splinesketch {
struct PaperCertificateInspector {
  using Sketch = CertifiedPaperSplineSketch;
  using Node = Sketch::Node;
  using Item = std::pair<double, std::uint64_t>;
  inline static std::uint64_t checks = 0;
  inline static bool export_queries = false;
  static std::uint64_t bits(double value) {
    std::uint64_t result; std::memcpy(&result, &value, sizeof result); return result;
  }

  static std::vector<double> queries(const Sketch& sketch, const std::map<double, std::uint64_t>& truth) {
    std::vector<double> xs{-INFINITY, INFINITY};
    for (auto item : truth) {
      xs.push_back(item.first);
      xs.push_back(std::nextafter(item.first, -INFINITY));
      xs.push_back(std::nextafter(item.first, INFINITY));
    }
    for (std::size_t i = 0; i < sketch.nodes_.size(); ++i) {
      xs.push_back(sketch.nodes_[i].x);
      if (i) xs.push_back(Sketch::midpoint(sketch.nodes_[i - 1].x, sketch.nodes_[i].x));
    }
    return xs;
  }
  static void check(const Sketch& sketch, const std::map<double, std::uint64_t>& truth) {
    const auto bound = sketch.max_rank_error();
    assert(bound <= static_cast<double>(sketch.max_rank_uncertainty()));
    for (std::size_t i = 1; i < sketch.nodes_.size(); ++i) {
      const auto& left = sketch.nodes_[i - 1]; const auto& right = sketch.nodes_[i];
      assert(left.lower + right.atom <= right.lower);
      assert(left.upper + right.atom <= right.upper);
    }
    for (double x : queries(sketch, truth)) {
      std::uint64_t rank = 0;
      for (auto item : truth) if (item.first <= x) rank += item.second;
      const auto answer = sketch.rank_with_error(x);
      assert(answer.lower_rank <= rank && rank <= answer.upper_rank);
      const auto error = std::fabs(static_cast<long double>(answer.estimate) - rank);
      assert(error <= bound && error <= answer.max_error);
      if (export_queries)
        std::cout << bits(answer.estimate) << ',' << rank << ',' << answer.lower_rank << ','
                  << answer.upper_rank << ',' << bits(answer.max_error) << ',' << bits(bound) << '\n';
      ++checks;
    }
  }
  static std::uint64_t heap_width(const Sketch& sketch, const Sketch::HeapRebalance& heap) {
    auto ordered = sketch; ordered.nodes_ = heap.snapshot_nodes();
    return ordered.max_rank_uncertainty();
  }
  static Sketch operators(std::mt19937_64& random, std::map<double, std::uint64_t>& truth,
                          bool heaps, int exponent) {
    Sketch sketch(128, Sketch::BoundPolicy::theoretical);
    for (int i = 0; i < 14; ++i) {
      const double x = std::ldexp(static_cast<double>(i * 8 - 48), exponent);
      const auto weight = random() % 7 + 1;
      sketch.nodes_.emplace_back(x, static_cast<long double>(weight));
      sketch.count_ += weight; truth[x] = weight;
      auto& node = sketch.nodes_.back();
      node.lower = node.upper = sketch.count_; node.atom = weight;
    }
    Sketch::rebuild(sketch.nodes_);
    // Forget some exact cuts first, so subsequent splits inherit nonzero
    // uncertainty. The large theoretical factor makes these joins legal.
    for (unsigned j = 0; j < 5; ++j) {
      const auto at = 1 + random() % (sketch.nodes_.size() - 2);
      const auto previous = sketch.max_rank_uncertainty();
      sketch.join_at(at);
      assert(sketch.max_rank_uncertainty() >= previous);
    }
    const auto old_width = sketch.max_rank_uncertainty();
    std::vector<Item> incoming;
    for (int i : {-50, -31, -7, 13, 29, 59}) {
      const double x = std::ldexp(static_cast<double>(i), exponent);
      const auto weight = random() % 5 + 1;
      incoming.emplace_back(x, weight); sketch.count_ += weight; truth[x] += weight;
    }
    Sketch::Snapshot before{sketch.nodes_, incoming, {}}; before.prepare();
    std::uint64_t added = 0; for (auto item : incoming) added += item.second;
    if (!heaps) {
      sketch.nodes_.insert(sketch.nodes_.begin(), Node{incoming.front().first});
      sketch.nodes_.emplace_back(incoming.back().first);
      sketch.reestimate(before);
      assert(sketch.max_rank_uncertainty() <= old_width + added);
      for (unsigned j = 0; j < 4; ++j) {
        sketch.clear_protection();
        const auto at = 1 + random() % (sketch.nodes_.size() - 2);
        const auto previous = sketch.max_rank_uncertainty();
        sketch.join_at(at);
        assert(sketch.max_rank_uncertainty() >= previous);
        const auto split_at = 1 + random() % (sketch.nodes_.size() - 1);
        const auto joined_width = sketch.max_rank_uncertainty();
        sketch.split_at(split_at, before);
        assert(sketch.max_rank_uncertainty() <= joined_width);
      }
      Sketch::rebuild(sketch.nodes_);
    } else {
      Sketch::HeapRebalance heap(sketch, before, false);
      heap.reestimate(before);
      assert(heap_width(sketch, heap) <= old_width + added);
      for (unsigned j = 0; j < 4; ++j) {
        sketch.clear_protection(); heap.rebuild_heaps();
        auto id = heap.links[heap.head].next;
        const auto previous = heap_width(sketch, heap);
        heap.join(id);
        const auto joined_width = heap_width(sketch, heap);
        assert(joined_width >= previous);
        id = heap.links[heap.head].next;
        heap.split(id);
        assert(heap_width(sketch, heap) <= joined_width);
      }
      heap.finish();
    }
    check(sketch, truth);
    return sketch;
  }
  static void composition() {
    std::mt19937_64 random(0x574944544850524fULL);
    for (unsigned trial = 0; trial < 200; ++trial) {
      std::map<double, std::uint64_t> left_truth, right_truth;
      const int exponent = trial % 3 == 0 ? -500 : trial % 3 == 1 ? 500 : 0;
      auto left = operators(random, left_truth, trial % 2 == 0, exponent);
      auto right = operators(random, right_truth, trial % 2 != 0, exponent);
      auto merged = left;
      merged.merge(right); // union fits, so this checks merge before compression
      assert(merged.max_rank_uncertainty() <= left.max_rank_uncertainty() + right.max_rank_uncertainty());
      auto truth = left_truth;
      for (auto item : right_truth) truth[item.first] += item.second;
      check(merged, truth);
      for (double x : queries(merged, truth)) {
        const auto a = left.rank_bounds(x), b = right.rank_bounds(x), c = merged.rank_bounds(x);
        assert(c.lower == a.lower + b.lower && c.upper == a.upper + b.upper);
      }
    }
  }
  static void noncontracting_split() {
    Sketch sketch(6, Sketch::BoundPolicy::theoretical);
    sketch.count_ = 100; sketch.nodes_ = {{0, 1}, {1, 99}};
    sketch.nodes_[0].lower = sketch.nodes_[0].upper = 1; sketch.nodes_[0].atom = 1;
    sketch.nodes_[1].lower = sketch.nodes_[1].upper = 100; sketch.nodes_[1].atom = 1;
    Sketch::rebuild(sketch.nodes_);
    Sketch::Snapshot before{sketch.nodes_, {}, {}};
    const auto width = sketch.max_rank_uncertainty();
    sketch.split_at(1, before);
    assert(width == 98 && sketch.max_rank_uncertainty() == width);
    assert(sketch.nodes_[1].upper - sketch.nodes_[1].lower == width);
    assert(sketch.nodes_[1].mass < before.nodes[1].mass);
    assert(sketch.nodes_[2].mass < before.nodes[1].mass);
  }
  static void uniform_fallback() {
    Sketch sketch(6); sketch.count_ = 100;
    sketch.nodes_ = {{0, 1}, {1, 14}, {2, 70}, {3, 15}};
    const std::uint64_t lower[]{1, 10, 80, 100}, upper[]{1, 20, 90, 100};
    for (unsigned i = 0; i < 4; ++i) {
      sketch.nodes_[i].lower = lower[i]; sketch.nodes_[i].upper = upper[i];
      sketch.nodes_[i].atom = i == 0 || i == 3 ? 1 : 0;
    }
    Sketch::rebuild(sketch.nodes_);
    const std::map<double, std::uint64_t> truth{{0, 1}, {0.5, 14}, {1.5, 70}, {3, 15}};
    const auto bound = sketch.max_rank_error();
    assert(bound < static_cast<double>(sketch.max_rank_uncertainty()));
    check(sketch, truth);
    sketch.nodes_[2].mass = std::numeric_limits<long double>::quiet_NaN();
    assert(sketch.rank(1.25) == 50); // midpoint of the integer [10,90] interval
    check(sketch, truth);
    sketch.nodes_[2].prefix = -1;
    assert(sketch.max_rank_error() == static_cast<double>(sketch.max_rank_uncertainty()));
    sketch.nodes_[2].prefix = std::numeric_limits<long double>::infinity();
    assert(sketch.max_rank_error() == static_cast<double>(sketch.max_rank_uncertainty()));
    sketch.nodes_[2].prefix = std::numeric_limits<long double>::quiet_NaN();
    assert(sketch.max_rank_error() == static_cast<double>(sketch.max_rank_uncertainty()));
  }
  static void witness() {
    Sketch sketch(32); std::map<double, std::uint64_t> truth;
    for (auto bits : paper_width_counterexample) {
      double x; std::memcpy(&x, &bits, sizeof x);
      sketch.add(x); ++truth[x];
    }
    sketch.consolidate(); check(sketch, truth);
    assert(sketch.count() == 288);
    assert(sketch.max_rank_uncertainty() * 32 >= 12 * sketch.count());
    assert(sketch.factor_ == 3);
    double largest = 0;
    for (double x : queries(sketch, truth)) {
      std::uint64_t rank = 0; for (auto item : truth) if (item.first <= x) rank += item.second;
      largest = std::max(largest, std::abs(sketch.rank(x) - static_cast<double>(rank)));
    }
    assert(largest > 3 * static_cast<double>(sketch.count()) / 32);
    auto repeated = sketch; auto repeated_truth = truth;
    for (unsigned i = 1; i <= 40; ++i) {
      repeated.merge(repeated);
      for (auto& item : repeated_truth) item.second *= 2;
      if (i % 5 == 0) check(repeated, repeated_truth);
    }
    const auto point_width = [&] {
      std::uint64_t width = 0;
      for (const auto& node : sketch.nodes_) width = std::max(width, node.upper - node.lower);
      return width;
    }();
    auto gradual = sketch;
    const auto anchor = [&] {
      const Node* widest = nullptr;
      for (const auto& node : gradual.nodes_)
        if (node.protected_threshold && (!widest || node.upper - node.lower > widest->upper - widest->lower))
          widest = &node;
      assert(widest && widest->upper > widest->lower);
      return *widest;
    }();
    for (unsigned i = 0; i < 12; ++i) {
      gradual.resize(gradual.capacity_ + gradual.capacity_ / 4);
      const auto at = std::find_if(gradual.nodes_.begin(), gradual.nodes_.end(),
                                 [&](const Node& node) { return node.x == anchor.x; });
      assert(at != gradual.nodes_.end() && at->protected_threshold);
      assert(at->lower == anchor.lower && at->upper == anchor.upper);
      assert(gradual.max_rank_uncertainty() >= anchor.upper - anchor.lower);
      check(gradual, truth);
    }
    sketch.resize(1024);
    // This particular large growth retains the widest old point interval.
    // Only protected cuts in non-resetting growth have a general guarantee.
    assert(sketch.max_rank_uncertainty() >= point_width);
    check(sketch, truth);
  }
};
}

int main(int argc, char** argv) {
  using Inspector = splinesketch::PaperCertificateInspector;
  Inspector::export_queries = argc == 2 && std::string(argv[1]) == "--exact-queries";
  if (argc != 1 && !Inspector::export_queries) throw std::invalid_argument("usage: width-tests [--exact-queries]");
  Inspector::noncontracting_split(); Inspector::uniform_fallback();
  Inspector::witness(); Inspector::composition();
  std::cout << "Paper width: " << Inspector::checks << " query checks; split/join/merge inequalities passed\n";
}
