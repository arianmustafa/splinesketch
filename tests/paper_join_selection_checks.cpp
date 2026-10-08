// Built by the experiment runner against isolated candidate headers.
#include <splinesketch/paper_splinesketch.hpp>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <random>
#ifndef SPLINESKETCH_JOIN_EXPERIMENT
#error Select width (1) or guarded (2) experiment
#endif
namespace splinesketch {
struct PaperCertificateInspector {
  using Sketch = CertifiedPaperSplineSketch;
  static void run() {
    std::mt19937_64 random(0x4a4f494e50524f4fULL);
    std::size_t comparisons = 0, changed = 0;
    for (unsigned trial = 0; trial < 400; ++trial) {
      Sketch sketch(32, Sketch::BoundPolicy::theoretical);
      for (unsigned i = 0; i < 20; ++i) {
        const auto weight = (random() % 10 + 1) << (trial % 3 == 0 ? 42 : 0);
        sketch.count_ += weight;
        sketch.nodes_.emplace_back(static_cast<double>(i * i), weight);
        auto& node = sketch.nodes_.back();
        node.lower = node.upper = sketch.count_; node.atom = weight;
      }
      Sketch::rebuild(sketch.nodes_);
      for (unsigned i = 0; i < 7; ++i) sketch.join_at(1 + random() % (sketch.nodes_.size() - 2));
      for (auto& node : sketch.nodes_) node.protected_threshold = random() % 4 == 0;
      sketch.nodes_[1].protected_threshold = false;
      Sketch::Snapshot snapshot{sketch.nodes_, {{12.25, 1}, {80.25, 2}, {160.25, 1}}, {}};
      snapshot.prepare(); sketch.count_ += 4;
      auto heap_sketch = sketch;
      Sketch::HeapRebalance heap(heap_sketch, snapshot, false);
      sketch.reestimate(snapshot); heap.reestimate(snapshot);
      heap.verify();
      const auto candidates = sketch.join_candidates();
      // Check a free join and every possible split exclusion. The oracle
      // executes each legal deletion and measures its resulting integer width.
      for (std::size_t split = 0; split < sketch.nodes_.size(); ++split) {
        std::size_t legacy = Sketch::npos;
        long double score = 0;
        std::uint64_t optimum = std::numeric_limits<std::uint64_t>::max();
        for (std::size_t i = 1; i + 1 < sketch.nodes_.size(); ++i) {
          if (i == split || i + 1 == split || sketch.nodes_[i].protected_threshold) continue;
          if (!(sketch.nodes_[i].mass + sketch.nodes_[i + 1].mass <= 0.75L * sketch.bound())) continue;
          auto probe = sketch; probe.join_at(i);
          optimum = std::min(optimum, probe.max_rank_uncertainty());
          const auto cost = sketch.heuristic(i, true);
          if (legacy == Sketch::npos || cost < score) { legacy = i; score = cost; }
        }
        if (legacy == Sketch::npos) continue;
        const auto excluded = split ? split : Sketch::npos;
        auto selected = candidates.nonoverlapping(excluded);
        auto heap_selected = heap.compatible(heap.joins.first_three(),
            split ? static_cast<Sketch::HeapRebalance::Id>(split) : Sketch::HeapRebalance::none);
#if SPLINESKETCH_JOIN_EXPERIMENT == 2
        assert(selected == legacy);
        selected = sketch.certificate_join(selected, excluded);
        heap_selected = heap.certificate_join(heap_selected,
            split ? static_cast<Sketch::HeapRebalance::Id>(split) : Sketch::HeapRebalance::none);
#endif
        assert(selected != Sketch::npos && heap_selected != Sketch::HeapRebalance::none);
        assert(sketch.nodes_[selected].x == heap_sketch.nodes_[heap_selected].x);
        auto probe = sketch; probe.join_at(selected);
        assert(probe.max_rank_uncertainty() == optimum);
#if SPLINESKETCH_JOIN_EXPERIMENT == 2
        auto original = sketch; original.join_at(legacy);
        if (optimum == original.max_rank_uncertainty()) assert(selected == legacy);
#endif
        changed += selected != legacy;
        ++comparisons;
      }
    }
    assert(changed > 0 && comparisons > 4000);
    std::cout << "Join selection: " << comparisons << " deletion/exclusion checks, "
              << changed << " choices differ; scanner and heap agree\n";
  }
};
}
int main() { splinesketch::PaperCertificateInspector::run(); }
