// Paired resize experiment: identical inputs, common queries, exact rank truth.
#define SPLINESKETCH_TESTING
#include <splinesketch/splinesketch.hpp>
#define splinesketch previous_splinesketch
#include SPLINESKETCH_PREVIOUS_CORE
#include SPLINESKETCH_PREVIOUS_PAPER
#undef splinesketch
#include <splinesketch/paper_splinesketch.hpp>
#include "../../tests/paper_resize_history_fixture.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

static void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

// Check the local integer loss formula independently at the removed key and
// its two old gaps. The removed key attains the maximum width increase:
// (L_i-L_left)+(U_right-atom_right-U_i). No true-rank error claim follows.
#ifdef SPLINESKETCH_RESIZE_REGRID_ONLY
#define CHECK_REGRID_FLAG(s) require(!(s).rank_resize_, "resize flag escaped regrid")
#else
#define CHECK_REGRID_FLAG(s) ((void)0)
#endif
#define DEFINE_INSPECTOR(NS, CHECK_FLAG) \
namespace NS { \
struct PaperCertificateInspector { \
  static std::vector<double> queries(const CertifiedPaperSplineSketch& s) { \
    CHECK_FLAG; \
    std::vector<double> result; \
    for (std::size_t i = 0; i < s.nodes_.size(); ++i) { \
      const auto& node = s.nodes_[i]; \
      require(node.lower <= node.upper && node.atom <= node.lower, "invalid node bounds"); \
      result.push_back(node.x); \
      result.push_back(std::nextafter(node.x, -INFINITY)); \
      result.push_back(std::nextafter(node.x, INFINITY)); \
      if (i) { \
        const auto& left = s.nodes_[i - 1]; \
        require(left.lower + node.atom <= node.lower, "lower envelope decreased"); \
        require(left.upper + node.atom <= node.upper, "upper envelope decreased"); \
        result.push_back(CertifiedPaperSplineSketch::midpoint(left.x, node.x)); \
      } \
      if (i && i + 1 < s.nodes_.size()) { \
        const auto& left = s.nodes_[i - 1]; const auto& right = s.nodes_[i + 1]; \
        const auto joined = right.upper - right.atom - left.lower; \
        const auto point = node.upper - node.lower; \
        const auto left_gap = node.upper - node.atom - left.lower; \
        const auto right_gap = right.upper - right.atom - node.lower; \
        require(joined >= point && joined >= left_gap && joined >= right_gap, "negative join loss"); \
        const auto loss = joined - point; \
        require(loss == (node.lower - left.lower) + (right.upper - right.atom - node.upper), "loss identity"); \
        require(loss >= joined - left_gap && loss >= joined - right_gap, "point not maximal loss"); \
      } \
    } \
    return result; \
  } \
}; \
}
DEFINE_INSPECTOR(splinesketch, CHECK_REGRID_FLAG(s))
DEFINE_INSPECTOR(previous_splinesketch, (void)0)
#undef DEFINE_INSPECTOR
#undef CHECK_REGRID_FLAG

using Before = previous_splinesketch::CertifiedPaperSplineSketch;
using After = splinesketch::CertifiedPaperSplineSketch;
using BeforePlain = previous_splinesketch::PaperSplineSketch;
using AfterPlain = splinesketch::PaperSplineSketch;

static std::vector<double> input(unsigned shape, unsigned seed, unsigned k) {
  std::mt19937_64 random(seed);
  std::vector<double> values;
  for (unsigned i = 0; i < 8 * k + 13; ++i) {
    double x = 0;
    switch (shape) {
      case 0: x = std::ldexp(static_cast<double>(random() % 100000), -12); break;
      case 1: x = i % 13 ? std::ldexp(static_cast<double>(random() % 31), -20) :
                    (i % 2 ? -1 : 1) * static_cast<double>(random() % 1000 + 1); break;
      case 2: x = static_cast<double>(random() % 17) - 8; break;
      case 3: x = i < k ? static_cast<double>(i) :
                    static_cast<double>(i - k + 1) * std::numeric_limits<double>::denorm_min(); break;
      case 4: x = i; break;
      case 5: x = -static_cast<double>(i); break;
      case 6: x = std::ldexp(1 + static_cast<double>(random() % 100000) / 100000,
                            static_cast<int>(random() % 2001) - 1000); break;
      case 7: x = (i / k % 2 ? 1 : -1) + std::ldexp(static_cast<double>(random() % 31), -20); break;
    }
    values.push_back(x);
  }
  return values;
}

struct Metrics {
  double before = 0, after = 0, before_query = 0, after_query = 0;
  std::size_t queries = 0;
};
static Metrics evaluate(const Before& before, const After& after,
                        const BeforePlain& plain_before, const AfterPlain& plain_after,
                        std::vector<double> data, bool identical = false) {
  std::sort(data.begin(), data.end());
  auto keys = data;
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  std::vector<double> queries{-INFINITY, INFINITY, paper_test::resize_history_query};
  for (std::size_t i = 0; i < keys.size(); ++i) {
    queries.push_back(keys[i]);
    queries.push_back(std::nextafter(keys[i], -INFINITY));
    queries.push_back(std::nextafter(keys[i], INFINITY));
    if (i) for (unsigned j = 1; j < 4; ++j) {
      const long double t = static_cast<long double>(j) / 4;
      queries.push_back(static_cast<double>((1 - t) * keys[i - 1] + t * keys[i]));
    }
  }
  for (const auto& xs : {previous_splinesketch::PaperCertificateInspector::queries(before),
                         splinesketch::PaperCertificateInspector::queries(after)})
    queries.insert(queries.end(), xs.begin(), xs.end());
  std::sort(queries.begin(), queries.end());
  queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
  Metrics result;
  result.queries = queries.size();
  const double before_bound = before.max_rank_error(), after_bound = after.max_rank_error();
  for (double x : queries) {
    const auto truth = static_cast<std::uint64_t>(std::upper_bound(data.begin(), data.end(), x) - data.begin());
    const auto check = [&](const auto& sketch, double bound) {
      const auto r = sketch.rank_with_error(x);
      const auto error = std::fabs(r.estimate - static_cast<double>(truth));
      require(r.lower_rank <= truth && truth <= r.upper_rank, "certificate excludes true rank");
      require(std::isfinite(r.estimate) && r.estimate == sketch.rank(x), "invalid estimate");
      require(error <= r.max_error && error <= bound, "error allowance too small");
      return error;
    };
    const auto before_error = check(before, before_bound), after_error = check(after, after_bound);
    if (before_error > result.before) { result.before = before_error; result.before_query = x; }
    if (after_error > result.after) { result.after = after_error; result.after_query = x; }
    if (identical) {
      const auto a = before.rank_with_error(x);
      const auto b = after.rank_with_error(x);
      require(a.estimate == b.estimate && a.lower_rank == b.lower_rank && a.upper_rank == b.upper_rank &&
              a.max_error == b.max_error, "candidate changed initial query result");
    }
    require(plain_before.rank(x) == plain_after.rank(x), "plain alias changed");
  }
  for (const auto count : {before.count(), after.count(), plain_before.count(), plain_after.count()})
    require(count == data.size(), "count changed");
  return result;
}

template<class S> static S replay(const std::vector<double>& data, unsigned k, bool theory) {
  S s(k, theory ? S::BoundPolicy::theoretical : S::BoundPolicy::practical);
  for (double x : data) s.add(x);
  s.consolidate();
  return s;
}

static std::size_t check_finalize(Before before, After after, BeforePlain plain_before,
                                 AfterPlain plain_after, const std::vector<double>& data, bool identical) {
  before.finalize(); after.finalize(); plain_before.finalize(); plain_after.finalize();
  require(before.finalized() && after.finalized() && plain_before.finalized() && plain_after.finalized(),
          "finalization failed");
  return evaluate(before, after, plain_before, plain_after, data, identical).queries;
}

static void compare(const std::string& group, unsigned shape, unsigned seed, unsigned k,
                    bool theory, bool shrink, std::vector<double> data, bool mixed = false) {
  auto before = replay<Before>(data, k, theory); auto after = replay<After>(data, k, theory);
  auto plain_before = replay<BeforePlain>(data, k, theory); auto plain_after = replay<AfterPlain>(data, k, theory);
  const auto initial = evaluate(before, after, plain_before, plain_after, data, true);
  require(initial.before == initial.after && before.max_rank_uncertainty() == after.max_rank_uncertainty(),
          "candidate changed initial streaming");
  const auto initial_bound = before.max_rank_error();
  std::uint64_t hash = 0xcbf29ce484222325ULL;
  for (double x : data) { std::uint64_t b; std::memcpy(&b, &x, sizeof b); hash = (hash ^ b) * 0x100000001b3ULL; }
  Metrics final = initial;
  double before_ns = 0, after_ns = 0;
  unsigned worse_steps = 0;
  std::size_t query_count = initial.queries + check_finalize(before, after, plain_before, plain_after, data, true);
  for (unsigned step = 0; step < 24; ++step) {
    const auto capacity = step % 2 ? k : shrink ? k / 2 : 2 * k;
    auto start = std::chrono::steady_clock::now(); before.resize(capacity);
    before_ns += std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now(); after.resize(capacity);
    after_ns += std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count();
    plain_before.resize(capacity); plain_after.resize(capacity);
    require(before.bucket_capacity() == capacity && after.bucket_capacity() == capacity &&
            plain_before.bucket_capacity() == capacity && plain_after.bucket_capacity() == capacity,
            "capacity change failed");
    if (mixed && step % 4 == 3) {
      const auto extra = input(shape, seed + step, 4);
      if (step % 8 == 3) {
        for (double x : extra) { before.add(x); after.add(x); plain_before.add(x); plain_after.add(x); }
      } else {
        before.merge(replay<Before>(extra, 16, theory)); after.merge(replay<After>(extra, 16, theory));
        plain_before.merge(replay<BeforePlain>(extra, 16, theory)); plain_after.merge(replay<AfterPlain>(extra, 16, theory));
      }
      data.insert(data.end(), extra.begin(), extra.end());
    }
    final = evaluate(before, after, plain_before, plain_after, data);
    query_count += final.queries;
    worse_steps += final.after > final.before;
    if (group == "exact_discrete")
      require(after.max_rank_uncertainty() == 0 && final.after == 0, "exact discrete resize regressed");
#ifdef SPLINESKETCH_RESIZE_REGRID_ONLY
    // The original sampled prototype exceeds 14 at the first shrink. The
    // scoped candidate keeps this entire practical replay below eight.
    if (group == "released_discrete" && !theory)
      require(final.after <= 8 && after.max_rank_uncertainty() <= 16,
              "released discrete resize regressed");
#endif
    if (group == "cluster_witness") {
      const auto r = after.rank_with_error(paper_test::resize_history_query);
      require(r.lower_rank == 4 && r.upper_rank == 9 && r.estimate == 9, "cluster witness regressed");
    }
  }
  query_count += check_finalize(before, after, plain_before, plain_after, data, false);
  std::cout << group << ',' << shape << ',' << seed << ',' << k << ',' << theory << ',' << shrink << ','
            << data.size() << ',' << hash << ',' << query_count << ',' << initial.before << ',' << initial_bound
            << ',' << final.before << ',' << final.after << ',' << before.max_rank_error() << ','
            << after.max_rank_error() << ',' << before.max_rank_uncertainty() << ',' << after.max_rank_uncertainty()
            << ',' << worse_steps << ',' << before_ns / 24 << ',' << after_ns / 24 << ','
            << final.before_query << ',' << final.after_query << '\n' << std::flush;
}

int main(int argc, char** argv) {
  const bool fresh = argc == 2 && std::string(argv[1]) == "--fresh";
  require(argc == 1 || fresh, "usage: paired [--fresh]");
  std::cout << std::setprecision(17)
            << "group,shape,seed,k,theoretical,shrink,n,input_hash,queries,initial_error,initial_bound,"
               "before_error,after_error,before_bound,after_bound,before_width,after_width,worse_steps,"
               "before_resize_ns,after_resize_ns,before_query,after_query\n";
  // Seeds 113/127/149 informed the narrower scope and are now development
  // data. The fresh seeds were chosen before observing the scoped candidate.
  const auto seeds = fresh ? std::array<unsigned, 3>{179, 191, 211} : std::array<unsigned, 3>{113, 127, 149};
  for (unsigned k : {16, 32, 64}) for (unsigned shape = 0; shape < 8; ++shape)
    for (unsigned seed : seeds) for (bool theory : {false, true}) for (bool shrink : {false, true})
      compare(fresh ? "fresh" : "development", shape, seed, k, theory, shrink, input(shape, seed, k));
  for (bool theory : {false, true}) {
    compare("exact_discrete", 2, 97, 32, theory, true, input(2, 97, 32));
    compare("released_discrete", 2, 113, 16, theory, true, input(2, 113, 16));
    compare("cluster_witness", 1, 0, 16, theory, false,
            {paper_test::resize_history_input.begin(), paper_test::resize_history_input.end()});
    for (unsigned shape = 0; shape < 3; ++shape)
      compare("mixed", shape, 163, 32, theory, true, input(shape, 163, 32), true);
  }
  return 0;
}
