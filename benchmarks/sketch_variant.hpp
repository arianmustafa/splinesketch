#pragma once

#if defined(SPLINESKETCH_PAPER_CERTIFIED)
#include <splinesketch/paper_splinesketch.hpp>
using Sketch = splinesketch::CertifiedPaperSplineSketch;
#elif defined(SPLINESKETCH_PAPER)
#include <splinesketch/paper_splinesketch.hpp>
#ifdef SPLINESKETCH_PAPER_THEORETICAL
struct Sketch : splinesketch::PaperSplineSketch {
  explicit Sketch(std::size_t k = 128) : splinesketch::PaperSplineSketch(k, BoundPolicy::theoretical) {}
};
#else
using Sketch = splinesketch::PaperSplineSketch;
#endif
#elif defined(SPLINESKETCH_CERTIFIED)
#include <splinesketch/certified_splinesketch.hpp>
using Sketch = splinesketch::CertifiedSplineSketch;
#else
#include <splinesketch/splinesketch.hpp>
using Sketch = splinesketch::SplineSketch;
#endif
