#include <splinesketch/paper_splinesketch.hpp>
#include "data/paper_width_counterexample.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <random>

namespace splinesketch {
struct PaperCertificateInspector {
  using Sketch = CertifiedPaperSplineSketch;
  using Item = std::pair<double, std::uint64_t>;
  inline static std::uint64_t comparisons = 0, transitions = 0;

  static std::uint64_t forgotten(const Sketch& sketch) {
    std::uint64_t known = sketch.buffer_.size();
    for (const auto& item : sketch.heavy_) known += item.second.exact;
    for (const auto& node : sketch.nodes_) known += node.atom;
    assert(known <= sketch.count());
    const auto budget = sketch.count() - known;
    assert(sketch.max_rank_uncertainty() <= budget);
    return budget;
  }

  // Bounds are step functions. Values and adjacent representable queries at
  // every breakpoint exhaust their common partition, including empty cells.
  static std::vector<double> partition(std::initializer_list<const Sketch*> sketches,
                                       const std::vector<Item>& incoming = {}) {
    std::vector<double> points{-INFINITY, INFINITY};
    const auto append = [&](double x) {
      points.push_back(x);
      points.push_back(std::nextafter(x, -INFINITY));
      points.push_back(std::nextafter(x, INFINITY));
    };
    for (const auto* sketch : sketches) {
      for (const auto& node : sketch->nodes_) append(node.x);
      for (double x : sketch->buffer_) append(x);
      for (const auto& item : sketch->heavy_) append(item.first);
    }
    for (const auto& item : incoming) append(item.first);
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());
    return points;
  }

  static void transition(const Sketch& before, const Sketch& after,
                         const std::vector<Item>& incoming = {}) {
    std::uint64_t added = 0;
    for (const auto& item : incoming) added += item.second;
    assert(after.count() == before.count() + added);
    for (double x : partition({&before, &after}, incoming)) {
      std::uint64_t shift = 0;
      for (const auto& item : incoming) if (item.first <= x) shift += item.second;
      const auto old = before.rank_bounds(x), current = after.rank_bounds(x);
      assert(current.lower <= old.lower + shift);
      assert(current.upper >= old.upper + shift);
      assert(current.upper - current.lower >= old.upper - old.lower);
      ++comparisons;
    }
    assert(after.max_rank_uncertainty() >= before.max_rank_uncertainty());
    assert(forgotten(after) >= forgotten(before));
    ++transitions;
  }

  static void merged(const Sketch& left, const Sketch& right, const Sketch& after) {
    assert(after.count() == left.count() + right.count());
    for (double x : partition({&left, &right, &after})) {
      const auto a = left.rank_bounds(x), b = right.rank_bounds(x), c = after.rank_bounds(x);
      assert(c.lower <= a.lower + b.lower);
      assert(c.upper >= a.upper + b.upper);
      assert(c.upper - c.lower >= (a.upper - a.lower) + (b.upper - b.lower));
      ++comparisons;
    }
    assert(after.max_rank_uncertainty() >= std::max(left.max_rank_uncertainty(),
                                                    right.max_rank_uncertainty()));
    assert(forgotten(after) >= forgotten(left) + forgotten(right));
    ++transitions;
  }

  static void projection() {
    // Exhaust the frozen-envelope projection lemma over weighted four-key
    // sources, optional retained cuts, inherited uncertain points, new exact
    // atoms, and every subset of a seven-point target grid. These synthetic
    // nodes test integer projection, without assuming a spline error bound.
    for (unsigned encoded = 0; encoded < 81; ++encoded) {
      unsigned digits = encoded;
      std::uint64_t weights[4];
      for (auto& weight : weights) { weight = digits % 3; digits /= 3; }
      ++weights[0]; ++weights[3];
      for (unsigned retained = 0; retained < 4; ++retained) {
        Sketch source(8, Sketch::BoundPolicy::theoretical);
        std::uint64_t prefix = 0;
        for (unsigned i = 0; i < 4; ++i) {
          prefix += weights[i];
          if (i == 0 || i == 3 || (retained & (1U << (i - 1)))) {
            source.nodes_.emplace_back(static_cast<double>(i));
            auto& node = source.nodes_.back();
            node.lower = node.upper = prefix; node.atom = weights[i];
            node.mass = static_cast<long double>(prefix) -
                (source.nodes_.size() == 1 ? 0 : source.nodes_[source.nodes_.size() - 2].prefix);
            node.prefix = static_cast<long double>(prefix);
          }
        }
        source.count_ = prefix;
        Sketch::rebuild(source.nodes_);
        for (unsigned refined = 0; refined < 2; ++refined) {
          if (refined) {
            auto old = source;
            Sketch::Snapshot frozen{source.nodes_, {}, {}};
            source.split_at(1, frozen);
            transition(old, source);
            // A split inside one frozen source cell exactly preserves its
            // old pointwise width, even at the newly introduced threshold.
            for (double x : partition({&old, &source})) {
              auto a = old.rank_bounds(x), b = source.rank_bounds(x);
              assert(a.lower == b.lower && a.upper == b.upper);
            }
          }
          for (unsigned atoms = 0; atoms < 3; ++atoms) {
            std::vector<Item> incoming{{0.5, atoms}, {2.5, 2 - atoms}};
            Sketch::Snapshot frozen{source.nodes_, incoming, {}};
            frozen.prepare();
            for (unsigned cuts = 0; cuts < 32; ++cuts) {
              Sketch target(8); target.count_ = source.count_ + 2;
              for (unsigned i = 0; i < 7; ++i) {
                if (i != 0 && i != 6 && !(cuts & (1U << (i - 1)))) continue;
                const double x = static_cast<double>(i) / 2;
                target.nodes_.emplace_back(x);
                const auto interval = frozen.bounds(x);
                auto& node = target.nodes_.back();
                node.lower = interval.lower; node.upper = interval.upper;
                node.atom = frozen.atom(x);
              }
              transition(source, target, incoming);
            }
          }
        }
      }
    }
  }

  static void public_histories() {
    std::mt19937_64 random(0x554e434552544149ULL);
    for (auto policy : {Sketch::BoundPolicy::practical, Sketch::BoundPolicy::theoretical}) {
      for (std::size_t capacity : {6U, 32U, 512U}) {
        Sketch sketch(capacity, policy);
        for (std::size_t i = 0; i < 12 * capacity + 17;) {
          const auto before = sketch;
          std::vector<Item> incoming;
          for (unsigned j = 0; j < 17 && i < 12 * capacity + 17; ++j, ++i) {
            const double x = i % 4 == 0 ? static_cast<double>(i % 7) :
                std::ldexp(static_cast<double>(random() % 1024) / 1024,
                           -static_cast<int>(i / (5 * capacity)));
            sketch.add(x); incoming.emplace_back(x, 1);
          }
          transition(before, sketch, incoming);
          if (i % 3 == 0) {
            const auto buffered = sketch;
            sketch.consolidate(); transition(buffered, sketch);
          }
        }
        auto before = sketch;
        sketch.consolidate(); transition(before, sketch);
        before = sketch;
        sketch.resize(2 * capacity + 3); transition(before, sketch);
        before = sketch;
        sketch.resize(capacity); transition(before, sketch);

        Sketch other(capacity + 3, policy);
        for (std::size_t i = 0; i < 8 * capacity + 11; ++i)
          other.add(static_cast<double>(random() % 257) - 128);
        before = sketch;
        sketch.merge(other); merged(before, other, sketch);
        before = sketch;
        sketch.merge(sketch); merged(before, before, sketch);
        assert(sketch.max_rank_uncertainty() >= 2 * before.max_rank_uncertainty());
        before = sketch;
        sketch.finalize(); transition(before, sketch);
        before = sketch;
        sketch.finalize(); transition(before, sketch);
      }
    }
  }

  static void witness() {
    Sketch sketch(32);
    for (auto bits : paper_width_counterexample) {
      double x; std::memcpy(&x, &bits, sizeof x); sketch.add(x);
    }
    sketch.consolidate();
    const auto initial = sketch.max_rank_uncertainty();
    assert(initial > 0);
    auto before = sketch;
    sketch.resize(1024); transition(before, sketch);
    before = sketch;
    sketch.resize(32); transition(before, sketch);
    for (unsigned i = 0; i < 55; ++i) {
      before = sketch;
      sketch.merge(sketch); merged(before, before, sketch);
      assert(sketch.max_rank_uncertainty() >= 2 * before.max_rank_uncertainty());
    }
    assert(sketch.count() > (std::uint64_t{1} << 53));
    // Doubling preserves a lower bound on normalized uncertainty exactly;
    // no floating conversion of a UINT64-sized width or count is involved.
    assert(sketch.max_rank_uncertainty() >= (initial << 55));
  }

  static void numeric_domain() {
    const double largest = std::numeric_limits<double>::max();
    const double tiny = std::numeric_limits<double>::denorm_min();
    const double values[]{-largest, -1e308, -tiny, -0.0, 0.0, tiny,
                          1.0, std::nextafter(1.0, INFINITY), 1e308, largest};
    Sketch sketch(6);
    for (unsigned i = 0; i < 100; ++i) {
      const auto before = sketch;
      const double x = values[i % 10];
      sketch.add(x); transition(before, sketch, {{x, 1}});
      if (i % 11 == 0) {
        const auto buffered = sketch;
        sketch.consolidate(); transition(buffered, sketch);
      }
    }
    const auto before = sketch;
    sketch.consolidate(); transition(before, sketch);
  }
};
}

int main() {
  using Inspector = splinesketch::PaperCertificateInspector;
  Inspector::projection();
  Inspector::public_histories();
  Inspector::witness();
  Inspector::numeric_domain();
  std::cout << "Uncertainty persistence: " << Inspector::transitions << " transitions, "
            << Inspector::comparisons << " partition comparisons passed\n";
}
