#pragma once

#include "splinesketch.hpp"
#include <array>

namespace splinesketch {

// Section 3 / Appendix B baseline, with Section 4.2's weighted MG batching.
// Kept separate from SplineSketch so comparisons retain the old algorithm.
// Transition checks are enforced in ordinary builds, not just in tests.
template<bool TrackBounds>
class BasicPaperSplineSketch {
#ifdef SPLINESKETCH_TESTING
  friend struct PaperSplineSketchInspector;
  friend struct PaperCertificateInspector;
#endif
  using CoreNode = SplineSketch::Node;
  struct BoundNode : CoreNode {
    using CoreNode::CoreNode;
    std::uint64_t lower = 0, upper = 0, atom = 0;
  };
  using Node = std::conditional_t<TrackBounds, BoundNode, CoreNode>;
  using Heavy = SplineSketch::Heavy;
  using Item = std::pair<double, std::uint64_t>;
  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

 public:
  // Practical follows Section 4.1's adaptive constant. Theoretical fixes 1024,
  // exceeding Appendix B's 45/beta for ideal midpoint arithmetic.
  enum class BoundPolicy { practical, theoretical };
  explicit BasicPaperSplineSketch(std::size_t buckets = 128,
      BoundPolicy policy = BoundPolicy::practical)
      : capacity_(checked_capacity(buckets)), epoch_end_(5 * buckets),
        policy_(policy), factor_(initial_factor()) {
    nodes_.reserve(capacity_ + 2);
    buffer_.reserve(5 * capacity_);
    heavy_.reserve(capacity_);
  }

  BasicPaperSplineSketch(const BasicPaperSplineSketch&) = default;
  BasicPaperSplineSketch(BasicPaperSplineSketch&&) noexcept = default;
  BasicPaperSplineSketch& operator=(BasicPaperSplineSketch&&) noexcept = default;
  BasicPaperSplineSketch& operator=(const BasicPaperSplineSketch& other) {
    if (this != &other) {
      BasicPaperSplineSketch updated(other);
      *this = std::move(updated);
    }
    return *this;
  }

  void add(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("value must be finite");
    if (finalized_) throw std::logic_error("cannot update a finalized sketch");
    if (count_ == std::numeric_limits<std::uint64_t>::max())
      throw std::overflow_error("observation count overflow");
    value = SplineSketch::canonical(value);
    if (buffer_.size() + 1 < 5 * capacity_) {
      buffer_.push_back(value);
      ++count_;
      return;
    }
    BasicPaperSplineSketch updated(*this);
    updated.buffer_.reserve(5 * capacity_);
    updated.buffer_.push_back(value);
    ++updated.count_;
    updated.consolidate_in_place();
    *this = std::move(updated);
  }

  std::uint64_t count() const noexcept { return count_; }
  std::size_t bucket_capacity() const noexcept { return capacity_; }
  std::size_t bucket_count() const noexcept { return nodes_.size(); }
  std::size_t heavy_hitter_count() const noexcept { return heavy_.size(); }
  bool finalized() const noexcept { return finalized_; }

  struct RankBounds { std::uint64_t lower, upper; };
  struct RankEstimate {
    double estimate;
    std::uint64_t lower_rank, upper_rank;
    double max_error;
  };
  template<bool Enabled = TrackBounds, std::enable_if_t<Enabled && TrackBounds, int> = 0>
  RankBounds rank_bounds(double x) const {
    if (std::isnan(x)) throw std::invalid_argument("query must not be NaN");
    auto interval = bucket_bounds(nodes_, x);
    std::uint64_t exact = 0;
    for (double value : buffer_) if (value <= x) ++exact;
    for (const auto& item : heavy_) if (item.first <= x) exact += item.second.exact;
    interval.lower += exact; interval.upper += exact;
    return interval;
  }
  template<bool Enabled = TrackBounds, std::enable_if_t<Enabled && TrackBounds, int> = 0>
  RankEstimate rank_with_error(double x) const {
    const auto result = rank_and_bounds(x);
    const auto interval = result.second;
    const double estimate = certified_estimate(result.first, interval);
    const auto error = std::max(integer_distance(estimate, interval.lower),
                                integer_distance(estimate, interval.upper));
    return {estimate, interval.lower, interval.upper, error};
  }
  // Uniform bound over all returned ranks at the current state. Certificate
  // uncertainty is an exact integer; the return allowance covers binary64
  // conversion even above 2^53. No spline monotonicity premise is needed.
  template<bool Enabled = TrackBounds, std::enable_if_t<Enabled && TrackBounds, int> = 0>
  std::uint64_t max_rank_uncertainty() const noexcept {
    std::uint64_t width = 0;
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
      width = std::max(width, nodes_[i].upper - nodes_[i].lower);
      if (i && std::nextafter(nodes_[i - 1].x, INFINITY) < nodes_[i].x)
        width = std::max(width, nodes_[i].upper - nodes_[i].atom - nodes_[i - 1].lower);
    }
    return width;
  }
  template<bool Enabled = TrackBounds, std::enable_if_t<Enabled && TrackBounds, int> = 0>
  double max_rank_error() const noexcept {
    const auto width = max_rank_uncertainty();
    if (count_ <= (std::uint64_t{1} << 53)) return prefix_error_bound(width);
    const double rounded_count = static_cast<double>(count_);
    const double allowance = std::nextafter(rounded_count, INFINITY) - rounded_count;
    return round_up(outward_integer(width, true) + allowance);
  }

  double rank(double x) const {
    if constexpr (TrackBounds) {
      const auto result = rank_and_bounds(x);
      return certified_estimate(result.first, result.second);
    } else return raw_rank(x);
  }

 private:
  double raw_rank(double x) const {
    if (std::isnan(x)) throw std::invalid_argument("query must not be NaN");
    long double value = spline_rank(nodes_, x);
    for (double v : buffer_) if (v <= x) value += 1;
    for (const auto& item : heavy_) if (item.first <= x) value += item.second.exact;
    return static_cast<double>(std::clamp(value, 0.0L, static_cast<long double>(count_)));
  }
  // The certified estimate and integer contribution share one exact-data scan;
  // floating additions retain precisely the original rank() order.
  std::pair<double, RankBounds> rank_and_bounds(double x) const {
    if (std::isnan(x)) throw std::invalid_argument("query must not be NaN");
    long double value = spline_rank(nodes_, x);
    std::uint64_t exact = 0;
    for (double v : buffer_) if (v <= x) { value += 1; ++exact; }
    for (const auto& item : heavy_) if (item.first <= x) {
      value += item.second.exact; exact += item.second.exact;
    }
    auto interval = bucket_bounds(nodes_, x);
    interval.lower += exact; interval.upper += exact;
    return {static_cast<double>(std::clamp(value, 0.0L, static_cast<long double>(count_))), interval};
  }

 public:
  double quantile(double q) const {
    if (!(q >= 0 && q <= 1)) throw std::invalid_argument("q must be in [0, 1]");
    if (!count_) throw std::logic_error("quantile of empty sketch");
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    if (!nodes_.empty()) { lo = nodes_.front().x; hi = nodes_.back().x; }
    for (double v : buffer_) { lo = std::min(lo, v); hi = std::max(hi, v); }
    for (const auto& v : heavy_) { lo = std::min(lo, v.first); hi = std::max(hi, v.first); }
    if (q == 0) return lo;
    if (q == 1) return hi;
    const auto target = detail::quantile_rank_target(q, count_);
    const auto certified_reaches = [&](double estimate, RankBounds interval) {
      // Integer certificates settle crossings which the rounded display can
      // move across a target. Inside the interval retain the spline policy.
      if (interval.lower >= target) return true;
      if (interval.upper < target) return false;
      return detail::rank_reaches_target(certified_estimate(estimate, interval), target);
    };
    const auto query_reaches = [&](double x) {
      if constexpr (TrackBounds) {
        const auto result = rank_and_bounds(x);
        return certified_reaches(result.first, result.second);
      } else return detail::rank_reaches_target(rank(x), target);
    };
    const std::size_t exact_size = buffer_.size() + heavy_.size();
    auto low = SplineSketch::ordered_bits(lo), high = SplineSketch::ordered_bits(hi);
    if (exact_size <= 16) {
      while (low < high) {
        const auto mid = low + (high - low) / 2;
        if (query_reaches(SplineSketch::from_ordered_bits(mid))) high = mid;
        else low = mid + 1;
      }
      return SplineSketch::from_ordered_bits(low);
    }
    Snapshot exact{{}, {}, {}};
    exact.items.reserve(exact_size);
    for (double x : buffer_) exact.items.emplace_back(x, 1);
    for (const auto& item : heavy_) exact.items.emplace_back(item.first, item.second.exact);
    exact.prepare();
    const long double slack = static_cast<long double>(count_) *
        (4 * (static_cast<long double>(exact_size) + 1) * std::numeric_limits<long double>::epsilon() +
         std::numeric_limits<double>::epsilon());
    const auto reaches = [&](double x) {
      const auto value = spline_rank(nodes_, x) + exact.rank(x);
      if (std::fabs(value - target) <= slack) return query_reaches(x);
      const auto estimate = static_cast<double>(std::clamp(value, 0.0L, static_cast<long double>(count_)));
      if constexpr (TrackBounds) {
        auto interval = bucket_bounds(nodes_, x);
        const auto added = exact.integer_rank(x);
        interval.lower += added; interval.upper += added;
        const auto certified = certified_estimate(estimate, interval);
        if (std::fabs(certified - target) <= slack) return query_reaches(x);
        return certified_reaches(estimate, interval);
      } else return detail::rank_reaches_target(estimate, target);
    };
    while (low < high) {
      const auto mid = low + (high - low) / 2;
      if (reaches(SplineSketch::from_ordered_bits(mid))) high = mid;
      else low = mid + 1;
    }
    return SplineSketch::from_ordered_bits(low);
  }

  void consolidate() {
    if (buffer_.empty()) return;
    BasicPaperSplineSketch updated(*this);
    updated.buffer_.reserve(5 * capacity_);
    updated.consolidate_in_place();
    *this = std::move(updated);
  }

  void resize(std::size_t buckets) {
    checked_capacity(buckets);
    if (finalized_) throw std::logic_error("cannot resize a finalized sketch");
    if (buckets == capacity_) return;
    BasicPaperSplineSketch updated(*this);
    updated.consolidate_in_place();
    const auto old = updated.capacity_;
    // A small shrink can also invalidate the practical extrema reserve.
    // Treat that resize as a protection boundary before incorporating data.
    const bool reserve_exhausted = policy_ == BoundPolicy::practical && buckets < old &&
        static_cast<std::size_t>(std::count_if(updated.nodes_.begin(), updated.nodes_.end(),
            [](const Node& node) { return node.protected_threshold; })) > buckets - 2;
    const bool reset_protection = reserve_exhausted ||
        (buckets > old ? buckets - old : old - buckets) > old / 4;
    if (reset_protection)
      updated.clear_protection();
    updated.capacity_ = buckets;
    updated.nodes_.reserve(buckets + 2);
    updated.buffer_.reserve(5 * buckets);
    updated.heavy_.reserve(buckets);
    std::vector<Item> released;
    updated.shrink_heavy(released);
    updated.incorporate(released);
    updated.rebalance({}, false, buckets > old);
    // Finish the protection boundary after refinement at the new capacity.
    if (reset_protection) updated.clear_protection();
    updated.trim_storage();
    *this = std::move(updated);
  }

  void merge(const BasicPaperSplineSketch& other) {
    if (finalized_ || other.finalized_)
      throw std::logic_error("cannot merge a finalized sketch");
    if (count_ > std::numeric_limits<std::uint64_t>::max() - other.count_)
      throw std::overflow_error("observation count overflow");
    const auto& larger = count_ >= other.count_ ? *this : other;
    BasicPaperSplineSketch result(larger.capacity_, larger.policy_);
    result.count_ = count_ + other.count_;
    result.epoch_end_ = larger.epoch_end_;
    result.factor_ = std::max(factor_, other.factor_);
    // A sorted threshold union; do not flush either input's raw buffer first.
    std::size_t a = 0, b = 0;
    while (a < nodes_.size() || b < other.nodes_.size()) {
      const double x = b == other.nodes_.size() ||
          (a < nodes_.size() && nodes_[a].x < other.nodes_[b].x)
          ? nodes_[a].x : other.nodes_[b].x;
      result.nodes_.emplace_back(x);
      if (a < nodes_.size() && nodes_[a].x == x) ++a;
      if (b < other.nodes_.size() && other.nodes_[b].x == x) ++b;
    }
    std::size_t at = 0;
    long double previous = 0;
    for (auto& node : result.nodes_) {
      // Rounded cubic samples can decrease at neighboring cuts. Carry the
      // previous target forward so a later recovery does not add that dip
      // again as mass. The last target still equals the two endpoint totals.
      const long double value = std::max(previous,
          spline_rank(nodes_, node.x) + spline_rank(other.nodes_, node.x));
      node.mass = value - previous;
      if constexpr (TrackBounds) {
        const auto left = bucket_bounds(nodes_, node.x), right = bucket_bounds(other.nodes_, node.x);
        node.lower = left.lower + right.lower; node.upper = left.upper + right.upper;
        node.atom = bucket_atom(nodes_, node.x) + bucket_atom(other.nodes_, node.x);
      }
      previous = value;
      while (at < larger.nodes_.size() && larger.nodes_[at].x < node.x) ++at;
      if (at < larger.nodes_.size() && larger.nodes_[at].x == node.x)
        node.protected_threshold = larger.nodes_[at].protected_threshold;
    }
    rebuild(result.nodes_);
    result.advance_epoch();
    std::vector<Item> released;
    // Combine equal MG keys before reducing the summary.
    for (const auto& item : heavy_) result.heavy_.emplace(item);
    for (const auto& item : other.heavy_) {
      auto& h = result.heavy_[item.first];
      h.residual += item.second.residual;
      h.exact += item.second.exact;
    }
    result.shrink_heavy(released);
    result.buffer_.insert(result.buffer_.end(), buffer_.begin(), buffer_.end());
    result.buffer_.insert(result.buffer_.end(), other.buffer_.begin(), other.buffer_.end());
    result.filter_buffer(released);
    result.incorporate(released);
    result.rebalance({}, false);
    result.trim_storage();
    *this = std::move(result);
  }

  // Section 4.2's final storage reduction. A finalized sketch is query-only.
  void finalize() {
    if (finalized_) return;
    BasicPaperSplineSketch updated(*this);
    updated.consolidate_in_place();
    std::vector<Item> released;
    for (auto it = updated.heavy_.begin(); it != updated.heavy_.end();) {
      if (static_cast<long double>(it->second.exact) <
          static_cast<long double>(count_) / (2 * capacity_)) {
        released.emplace_back(it->first, it->second.exact);
        it = updated.heavy_.erase(it);
      } else ++it;
    }
    updated.incorporate(released);
    updated.capacity_ = std::max({capacity_ - updated.heavy_.size(), capacity_ / 2, std::size_t{6}});
    updated.clear_protection(); // final resizing may be substantial
    updated.rebalance({}, false);
    updated.finalized_ = true;
    std::vector<double>().swap(updated.buffer_);
    updated.nodes_.shrink_to_fit();
    updated.heavy_.rehash(updated.heavy_.size());
    *this = std::move(updated);
  }

 private:
  struct Snapshot {
    std::vector<Node> nodes;
    std::vector<Item> items;
    std::vector<std::uint64_t> sums;
    void prepare() {
      std::sort(items.begin(), items.end());
      std::size_t out = 0;
      for (auto item : items) {
        if (out && items[out - 1].first == item.first) items[out - 1].second += item.second;
        else items[out++] = item;
      }
      items.resize(out);
      std::uint64_t sum = 0;
      sums.reserve(items.size());
      for (auto item : items) { sum += item.second; sums.push_back(sum); }
    }
    long double rank(double x) const {
      const auto it = std::upper_bound(items.begin(), items.end(), x,
          [](double value, const Item& item) { return value < item.first; });
      const auto at = static_cast<std::size_t>(it - items.begin());
      return spline_rank(nodes, x) + (at ? sums[at - 1] : 0);
    }
    std::uint64_t integer_rank(double x) const {
      const auto it = std::upper_bound(items.begin(), items.end(), x,
          [](double value, const Item& item) { return value < item.first; });
      const auto at = static_cast<std::size_t>(it - items.begin());
      return at ? sums[at - 1] : 0;
    }
    RankBounds bounds(double x) const {
      auto interval = bucket_bounds(nodes, x);
      const auto exact = integer_rank(x);
      interval.lower += exact; interval.upper += exact;
      return interval;
    }
    std::uint64_t atom(double x) const {
      const auto it = std::lower_bound(items.begin(), items.end(), x,
          [](const Item& item, double value) { return item.first < value; });
      return bucket_atom(nodes, x) + (it != items.end() && it->first == x ? it->second : 0);
    }
  };
  std::size_t capacity_;
  std::uint64_t count_ = 0, epoch_end_;
  BoundPolicy policy_;
  long double factor_;
  bool finalized_ = false;
  std::vector<Node> nodes_;
  std::vector<double> buffer_;
  std::unordered_map<double, Heavy> heavy_;

  static RankBounds bucket_bounds(const std::vector<Node>& nodes, double x) {
    if constexpr (TrackBounds) {
      if (nodes.empty() || x < nodes.front().x) return {0, 0};
      if (x >= nodes.back().x) return {nodes.back().upper, nodes.back().upper};
      const auto it = std::lower_bound(nodes.begin(), nodes.end(), x,
          [](const Node& node, double value) { return node.x < value; });
      if (it->x == x) return {it->lower, it->upper};
      return {(it - 1)->lower, it->upper - it->atom};
    } else return {0, 0};
  }
  static std::uint64_t bucket_atom(const std::vector<Node>& nodes, double x) {
    if constexpr (TrackBounds) {
      const auto it = std::lower_bound(nodes.begin(), nodes.end(), x,
          [](const Node& node, double value) { return node.x < value; });
      return it != nodes.end() && it->x == x ? it->atom : 0;
    } else return 0;
  }
  static double outward_integer(std::uint64_t value, bool up) noexcept {
    const double rounded = static_cast<double>(value);
    if (value <= (std::uint64_t{1} << 53)) return rounded;
    return std::nextafter(rounded, up ? INFINITY : -INFINITY);
  }
  static double round_up(double value) noexcept {
    return value == 0 ? 0 : std::nextafter(value, INFINITY);
  }
  static double integer_distance(double value, std::uint64_t endpoint) noexcept {
    // At and above 2^52 every binary64 value is an integer. Compute the
    // distance in uint64 instead of losing low count bits in a conversion.
    // The certified estimate is always in [0,2^64]; that final endpoint
    // needs a separate branch because conversion to uint64 would overflow.
    if (value >= 0x1p52) {
      if (value == 0x1p64) {
        if (!endpoint) return value;
        return outward_integer(std::numeric_limits<std::uint64_t>::max() - endpoint + 1, true);
      }
      const auto integer = static_cast<std::uint64_t>(value);
      return outward_integer(integer >= endpoint ? integer - endpoint : endpoint - integer, true);
    }
    return round_up(std::max(std::fabs(value - outward_integer(endpoint, false)),
                             std::fabs(outward_integer(endpoint, true) - value)));
  }
  static double certified_estimate(double value, RankBounds interval) noexcept {
    if (!std::isfinite(value))
      value = static_cast<double>(interval.lower + (interval.upper - interval.lower) / 2);
    return std::clamp(value, static_cast<double>(interval.lower), static_cast<double>(interval.upper));
  }
  // Intersect the integer-width guarantee with the actual spline-prefix
  // bracket. This does not assume monotone cubic evaluation: finite values
  // are clamped to their endpoint prefixes, and a nonfinite value uses the
  // integer midpoint. Exact-data additions and return conversion have a
  // conservative ULP allowance. Unsupported prefixes retain the width bound.
  double prefix_error_bound(std::uint64_t width) const noexcept {
    if constexpr (!TrackBounds) return static_cast<double>(width);
    else {
      const double fallback = static_cast<double>(width);
      if (!width || nodes_.empty()) return fallback;
      const auto additions = buffer_.size() + heavy_.size();
      if (additions > (std::uint64_t{1} << 48) - 3 ||
          !std::isfinite(nodes_.back().prefix) || nodes_.back().prefix < 0) return fallback;
      const double top = std::nextafter(static_cast<double>(nodes_.back().prefix), INFINITY);
      const double sum = round_up(top + static_cast<double>(count_));
      const double ceiling = round_up(2 * sum);
      if (!std::isfinite(ceiling)) return fallback;
      const double ulp = std::nextafter(ceiling, INFINITY) - ceiling;
      // h additions, one return conversion, one prefix conversion and one
      // endpoint-distance subtraction, each bounded by this common ULP.
      const double rounding = round_up(static_cast<double>(additions + 3) * ulp);
      if (!std::isfinite(rounding)) return fallback;
      double bound = 0, previous_prefix = 0;
      for (std::size_t i = 0; i < nodes_.size(); ++i) {
        const auto& node = nodes_[i];
        if (!std::isfinite(node.prefix) || node.prefix < 0 ||
            (i && node.prefix < nodes_[i - 1].prefix)) return fallback;
        const double prefix = static_cast<double>(node.prefix);
        const auto point_width = node.upper - node.lower;
        if (point_width) {
          const double point = std::max(std::fabs(prefix - static_cast<double>(node.upper)),
                                        std::fabs(prefix - static_cast<double>(node.lower)));
          bound = std::max(bound, std::min(static_cast<double>(point_width), round_up(point + rounding)));
        }
        if (i && std::nextafter(nodes_[i - 1].x, INFINITY) < node.x) {
          const auto lower = nodes_[i - 1].lower, upper = node.upper - node.atom;
          const auto interval_width = upper - lower;
          const double extent = std::max(std::fabs(previous_prefix - static_cast<double>(upper)),
                                         std::fabs(prefix - static_cast<double>(lower)));
          const double finite = std::min(static_cast<double>(interval_width), round_up(extent + rounding));
          const double nonfinite = static_cast<double>(interval_width / 2 + interval_width % 2);
          bound = std::max(bound, std::max(finite, nonfinite));
        }
        if (bound >= fallback) return fallback;
        previous_prefix = prefix;
      }
      return std::min(fallback, bound);
    }
  }
  // Reuse the ordered reestimation cursors instead of searching each key.
  static void set_bounds(Node& node, const Snapshot& before, std::size_t old, std::size_t pending) {
    if constexpr (TrackBounds) {
      RankBounds interval;
      std::uint64_t atom = 0;
      if (old < before.nodes.size() && before.nodes[old].x == node.x) {
        interval = {before.nodes[old].lower, before.nodes[old].upper};
        atom = before.nodes[old].atom;
      } else if (old == 0) interval = {0, 0};
      else if (old == before.nodes.size())
        interval = {before.nodes.back().upper, before.nodes.back().upper};
      else interval = {before.nodes[old - 1].lower, before.nodes[old].upper - before.nodes[old].atom};
      const auto added = pending ? before.sums[pending - 1] : 0;
      node.lower = interval.lower + added; node.upper = interval.upper + added;
      node.atom = atom + (pending && before.items[pending - 1].first == node.x
                         ? before.items[pending - 1].second : 0);
    }
  }

  static std::size_t checked_capacity(std::size_t k) {
    if (k < 6 || k > std::numeric_limits<std::size_t>::max() / 10)
      throw std::invalid_argument("bucket count must be at least 6 and fit the buffer");
    return k;
  }
  long double initial_factor() const { return policy_ == BoundPolicy::theoretical ? 1024 : 3; }
  long double bound() const { return factor_ * static_cast<long double>(count_) / capacity_; }
  void trim_storage() {
    // Merges and shrinks must not retain arrays sized for an old larger sketch.
    std::vector<Node> compacted;
    compacted.reserve(capacity_ + 2);
    compacted.insert(compacted.end(), nodes_.begin(), nodes_.end());
    nodes_.swap(compacted);
    std::vector<double> raw;
    raw.reserve(5 * capacity_);
    buffer_.swap(raw); // these callers have already flushed both raw buffers
    heavy_.rehash(capacity_);
  }
  void clear_protection() {
    for (auto& node : nodes_) node.protected_threshold = false;
  }
  void advance_epoch() {
    if (count_ <= epoch_end_) return;
    clear_protection();
    factor_ = initial_factor();
    while (epoch_end_ < count_) {
      const auto step = std::max<std::uint64_t>(1, epoch_end_ / 4 + (epoch_end_ % 4 != 0));
      if (epoch_end_ > std::numeric_limits<std::uint64_t>::max() - step) {
        epoch_end_ = std::numeric_limits<std::uint64_t>::max();
        break;
      }
      epoch_end_ += step;
    }
  }
  static double midpoint(double a, double b) {
    return a < 0 && b > 0 ? a / 2 + b / 2 : a + (b - a) / 2;
  }
  static long double fraction(double a, double b, double x) {
    long double h = static_cast<long double>(b) - a;
    if (std::isfinite(h)) return (static_cast<long double>(x) - a) / h;
    return (static_cast<long double>(x) / 2 - static_cast<long double>(a) / 2) /
        (static_cast<long double>(b) / 2 - static_cast<long double>(a) / 2);
  }
  static void rebuild(std::vector<Node>& nodes) {
    SplineSketch::rebuild(nodes);
  }
  static long double spline_rank(const std::vector<Node>& nodes, double x) {
    if (nodes.empty() || x < nodes.front().x) return 0;
    if (x >= nodes.back().x) return nodes.back().prefix;
    auto it = std::lower_bound(nodes.begin(), nodes.end(), x,
        [](const Node& n, double v) { return n.x < v; });
    if (it->x == x) return it->prefix;
    const auto& l = *(it - 1);
    const auto& r = *it;
    const long double t = fraction(l.x, r.x, x);
    const long double h = static_cast<long double>(r.x) - l.x;
    const long double dl = h * l.slope, dr = h * r.slope;
    // Fallback is decided per whole interval, preserving monotonicity.
    if (!std::isfinite(dl) || !std::isfinite(dr))
      return l.prefix + t * r.mass;
    const long double value = l.prefix + (3 * t * t - 2 * t * t * t) * r.mass +
        t * (1 - t) * ((1 - t) * dl - t * dr);
    return std::clamp(value, l.prefix, r.prefix);
  }
  void filter_buffer(std::vector<Item>& released) {
    std::sort(buffer_.begin(), buffer_.end());
    std::vector<Item> incoming;
    for (double x : buffer_) {
      if (!incoming.empty() && incoming.back().first == x) ++incoming.back().second;
      else incoming.emplace_back(x, 1);
    }
    // Stage 1: all keys already held by MG receive their whole batch frequency.
    for (auto& item : incoming) {
      const auto it = heavy_.find(item.first);
      if (it != heavy_.end()) {
        it->second.residual += item.second;
        it->second.exact += item.second;
        item.second = 0;
      }
    }
    // Stage 2: weighted Misra-Gries on new distinct keys, in value order.
    for (auto item : incoming) {
      if (!item.second) continue;
      Heavy h{item.second, item.second};
      while (h.residual && heavy_.size() >= capacity_ - 1) {
        std::uint64_t delta = h.residual;
        for (const auto& held : heavy_) delta = std::min(delta, held.second.residual);
        for (auto it = heavy_.begin(); it != heavy_.end();) {
          it->second.residual -= delta;
          if (!it->second.residual) {
            released.emplace_back(it->first, it->second.exact);
            it = heavy_.erase(it);
          } else ++it;
        }
        h.residual -= delta;
      }
      if (h.residual) heavy_.emplace(item.first, h);
      else released.emplace_back(item.first, h.exact);
    }
    buffer_.clear();
  }
  void shrink_heavy(std::vector<Item>& released) {
    while (heavy_.size() >= capacity_) {
      std::uint64_t delta = std::numeric_limits<std::uint64_t>::max();
      for (const auto& h : heavy_) delta = std::min(delta, h.second.residual);
      for (auto it = heavy_.begin(); it != heavy_.end();) {
        it->second.residual -= delta;
        if (!it->second.residual) {
          released.emplace_back(it->first, it->second.exact);
          it = heavy_.erase(it);
        } else ++it;
      }
    }
  }
  void consolidate_in_place() {
    advance_epoch();
    std::vector<Item> released;
    filter_buffer(released);
    incorporate(released);
  }
  void initialize(const Snapshot& before) {
    if (before.items.size() <= capacity_) {
      for (auto item : before.items) nodes_.emplace_back(item.first, item.second);
    } else {
      const auto total = before.sums.back();
      for (std::size_t i = 0; i < capacity_; ++i) {
        // Rounded weighted positions can miss the maximum even in extended
        // precision. Select both endpoints with exact integer offsets.
        std::uint64_t offset = 0;
        if (i == capacity_ - 1) offset = total - 1;
        else if (i != 0) {
          const long double position = static_cast<long double>(i) * (total - 1) / (capacity_ - 1);
          const auto rounded = std::floor(position + 0.5L);
          // Clamp before conversion: binary64 positions can round to 2^64.
          offset = rounded >= static_cast<long double>(total - 1)
              ? total - 1 : static_cast<std::uint64_t>(rounded);
        }
        const auto at = std::upper_bound(before.sums.begin(), before.sums.end(), offset);
        const auto index = static_cast<std::size_t>(at - before.sums.begin());
        if (nodes_.empty() || nodes_.back().x != before.items[index].first)
          nodes_.emplace_back(before.items[index].first);
      }
      // MG may leave fewer distinct quantile selections than the capacity.
      for (auto item : before.items) {
        if (nodes_.size() == capacity_) break;
        const auto it = std::lower_bound(nodes_.begin(), nodes_.end(), item.first,
            [](const Node& node, double x) { return node.x < x; });
        if (it == nodes_.end() || it->x != item.first) nodes_.insert(it, Node{item.first});
      }
      reestimate(before);
    }
    rebuild_prefixes();
    if constexpr (TrackBounds) {
      std::size_t pending = 0;
      for (auto& node : nodes_) {
        while (pending < before.items.size() && before.items[pending].first <= node.x) ++pending;
        set_bounds(node, before, 0, pending);
      }
    }
    // Complete initialization with exact buffer ranks. These thresholds are
    // initial selections, so they carry no split-history protection.
    while (nodes_.size() < capacity_) {
      std::size_t split = npos;
      long double width = -1;
      for (std::size_t i = 1; i < nodes_.size(); ++i) {
        const long double h = static_cast<long double>(nodes_[i].x) - nodes_[i - 1].x;
        if (can_split(i) && h > width) { split = i; width = h; }
      }
      if (split == npos) break;
      split_at(split, before);
    }
    rebuild(nodes_);
    clear_protection();
  }
  void reestimate(const Snapshot& before) {
    std::size_t old = 0, pending = 0;
    long double previous = 0;
    for (auto& node : nodes_) {
      while (old < before.nodes.size() && before.nodes[old].x < node.x) ++old;
      while (pending < before.items.size() && before.items[pending].first <= node.x) ++pending;
      const long double old_rank = old < before.nodes.size() && before.nodes[old].x == node.x
          ? before.nodes[old].prefix : spline_rank(before.nodes, node.x);
      const long double value = old_rank + (pending ? before.sums[pending - 1] : 0);
      node.mass = std::max(0.0L, value - previous);
      set_bounds(node, before, old, pending);
      previous = value;
    }
    rebuild_prefixes();
  }
  void incorporate(std::vector<Item>& released) {
    if (released.empty()) return;
    Snapshot all{{}, std::move(released), {}};
    all.prepare();
    if (nodes_.empty()) {
      initialize(all);
      return;
    }
    // Keep stable slots and heaps across all proof-required batches. Rebuilding
    // them for every small batch costs more than the scans they replace.
    if (use_heap_backend()) {
      incorporate_heaps(all);
      return;
    }
    // Appendix B bounds batch *occurrences*, not the number of distinct keys.
    const auto batch_limit = occurrence_batch_limit();
    std::size_t at = 0;
    std::uint64_t used = 0;
    while (at < all.items.size()) {
      Snapshot batch{nodes_, {}, {}};
      std::uint64_t size = 0;
      while (at < all.items.size() && size < batch_limit) {
        const auto take = std::min(all.items[at].second - used, batch_limit - size);
        batch.items.emplace_back(all.items[at].first, take);
        size += take;
        used += take;
        if (used == all.items[at].second) { ++at; used = 0; }
      }
      batch.prepare();
      if (batch.items.front().first < nodes_.front().x)
        nodes_.insert(nodes_.begin(), Node{batch.items.front().first});
      if (batch.items.back().first > nodes_.back().x)
        nodes_.emplace_back(batch.items.back().first);
      reestimate(batch);
      rebalance(batch, true);
    }
  }
  std::uint64_t occurrence_batch_limit() const {
    // floor(C_b*n/(2k)), saturated at n, without a widened integer type or
    // floating-point count conversion. Initial 3/1024 factors remain exact
    // integers after doubling or taking the maximum of merge inputs.
    const auto maximum = std::numeric_limits<std::uint64_t>::max();
    if (factor_ >= static_cast<long double>(maximum)) return std::max<std::uint64_t>(1, count_);
    const auto multiplier = static_cast<std::uint64_t>(factor_);
    const auto denominator = 2 * static_cast<std::uint64_t>(capacity_);
    if (multiplier >= denominator) return std::max<std::uint64_t>(1, count_);
    const auto whole = (count_ / denominator) * multiplier;
    const auto r = count_ % denominator;
    std::uint64_t quotient = 0, remainder = 0;
    unsigned bits = 0;
    for (auto m = multiplier; m; m >>= 1) ++bits;
    while (bits) {
      --bits;
      quotient *= 2;
      if (remainder >= denominator - remainder) {
        remainder -= denominator - remainder;
        ++quotient;
      } else remainder *= 2;
      if ((multiplier >> bits) & 1U) {
        if (remainder >= denominator - r) {
          remainder -= denominator - r;
          ++quotient;
        } else remainder += r;
      }
    }
    return std::max<std::uint64_t>(1, whole + quotient);
  }
  bool can_split(std::size_t i) const {
    if (!i || i >= nodes_.size()) return false;
    const auto mid = midpoint(nodes_[i - 1].x, nodes_[i].x);
    return mid > nodes_[i - 1].x && mid < nodes_[i].x;
  }
  static std::size_t split_protections(const Node& left, const Node& right) {
    return 1 + static_cast<std::size_t>(!left.protected_threshold) +
        static_cast<std::size_t>(!right.protected_threshold);
  }
  bool protection_budget_allows(std::size_t protected_count,
                               const Node& left, const Node& right) const {
    // P <= k-2 reserves enough slack plus removable cuts for two new extrema,
    // including when the grid has fewer than k nodes. Retain every old barrier.
    return policy_ == BoundPolicy::theoretical ||
        protected_count + split_protections(left, right) <= capacity_ - 2;
  }
  long double heuristic(std::size_t i, bool joined = false) const {
    const std::size_t right = i + static_cast<std::size_t>(joined);
    const long double h = static_cast<long double>(nodes_[right].x) - nodes_[i - 1].x;
    const long double hp = i > 1 ? static_cast<long double>(nodes_[i - 1].x) - nodes_[i - 2].x : h;
    const long double hn = right + 1 < nodes_.size() ?
        static_cast<long double>(nodes_[right + 1].x) - nodes_[right].x : h;
    const long double mass = nodes_[i].mass + (joined ? nodes_[right].mass : 0);
    const long double next = right + 1 < nodes_.size() ? nodes_[right + 1].mass : 0;
    const long double value = std::max(std::fabs(mass / h - nodes_[i - 1].mass / hp) / (h + hp),
        std::fabs(next / hn - mass / h) / (h + hn)) * h * h;
    // Extreme binary64 arithmetic may overflow. Ordering remains deterministic;
    // candidate legality never depends on this optional heuristic.
    return std::isfinite(value) ? value : std::numeric_limits<long double>::max();
  }
  struct JoinCandidates {
    using Candidate = std::pair<long double, std::size_t>;
    std::array<Candidate, 3> best{{{0, npos}, {0, npos}, {0, npos}}};
    std::size_t count = 0;
    std::size_t nonoverlapping(std::size_t split) const {
      for (const auto& candidate : best)
        if (candidate.second != npos && candidate.second != split && candidate.second + 1 != split)
          return candidate.second;
      return npos;
    }
  };
  JoinCandidates join_candidates() const {
    JoinCandidates result;
    const auto limit = 0.75L * bound();
    for (std::size_t i = 1; i + 1 < nodes_.size(); ++i) {
      if (nodes_[i].protected_threshold || !(nodes_[i].mass + nodes_[i + 1].mass <= limit)) continue;
      ++result.count;
      const typename JoinCandidates::Candidate candidate{heuristic(i, true), i};
      for (std::size_t j = 0; j < result.best.size(); ++j)
        if (result.best[j].second == npos || candidate < result.best[j]) {
          for (std::size_t at = result.best.size() - 1; at > j; --at) result.best[at] = result.best[at - 1];
          result.best[j] = candidate;
          break;
        }
    }
    // A split excludes only the two adjacent join pairs. The cheapest three
    // therefore contain exactly the first legal non-overlapping sorted join.
    return result;
  }
  void raise_bound(bool requires_removal) {
    if (policy_ == BoundPolicy::theoretical)
      throw std::logic_error("paper baseline cannot find a legal join");
    bool unprotected = false;
    for (std::size_t i = 1; i + 1 < nodes_.size(); ++i)
      unprotected |= !nodes_[i].protected_threshold;
    // A full grid can stop needing a mandatory split after the bound rises.
    // An oversized grid must still remove a cut, so protection can block it.
    if (requires_removal && !unprotected)
      throw std::logic_error("paper baseline has no unprotected join threshold");
    factor_ *= 2;
  }
  void join_at(std::size_t i) {
    if (nodes_[i].protected_threshold || nodes_[i].mass + nodes_[i + 1].mass > 0.75L * bound())
      throw std::logic_error("illegal paper join");
    nodes_[i + 1].mass += nodes_[i].mass;
    nodes_.erase(nodes_.begin() + static_cast<std::ptrdiff_t>(i));
    rebuild_prefixes();
  }
  void split_at(std::size_t i, const Snapshot& before) {
    const auto mid = midpoint(nodes_[i - 1].x, nodes_[i].x);
    const auto left = nodes_[i - 1].prefix, right = nodes_[i].prefix;
    const long double value = std::clamp(before.rank(mid), left, right);
    nodes_.insert(nodes_.begin() + static_cast<std::ptrdiff_t>(i), Node{mid, value - left, 0, 0, true});
    if constexpr (TrackBounds) {
      const auto interval = before.bounds(mid);
      nodes_[i].lower = interval.lower; nodes_[i].upper = interval.upper;
      nodes_[i].atom = before.atom(mid);
    }
    nodes_[i + 1].mass = right - value;
    nodes_[i - 1].protected_threshold = nodes_[i + 1].protected_threshold = true;
    rebuild_prefixes();
  }
  void rebuild_prefixes() {
    long double sum = 0;
    for (auto& node : nodes_) { sum += node.mass; node.prefix = sum; }
  }
  // Transient stable slots: nodes_ is repurposed while this private checked
  // copy is being consolidated. Its prefix/slope fields cache heap scores;
  // the frozen Snapshot remains the only interpolation source until finish().
  struct HeapRebalance {
    using Id = std::uint32_t;
    static constexpr Id none = std::numeric_limits<Id>::max();
    enum Kind : unsigned { mass_heap, split_heap, join_heap };
    struct Links {
      Id prev = none, next = none;
      std::array<Id, 3> position{{none, none, none}};
      bool dirty = false;
    };
    struct Heap {
      HeapRebalance& state;
      Kind kind;
      std::vector<Id> entries;
      long double key(Id id) const {
        const auto& node = state.sketch.nodes_[id];
        return kind == mass_heap ? node.mass : kind == split_heap ? node.prefix : node.slope;
      }
      bool better(Id a, Id b) const {
        const auto x = key(a), y = key(b);
        if (x != y) return kind == join_heap ? x < y : x > y;
        return state.sketch.nodes_[a].x < state.sketch.nodes_[b].x;
      }
      void exchange(std::size_t a, std::size_t b) {
        std::swap(entries[a], entries[b]);
        state.links[entries[a]].position[kind] = static_cast<Id>(a);
        state.links[entries[b]].position[kind] = static_cast<Id>(b);
      }
      std::size_t up(std::size_t at) {
        while (at && better(entries[at], entries[(at - 1) / 2])) {
          const auto parent = (at - 1) / 2;
          exchange(at, parent); at = parent;
        }
        return at;
      }
      void down(std::size_t at) {
        for (;;) {
          auto child = 2 * at + 1;
          if (child >= entries.size()) return;
          if (child + 1 < entries.size() && better(entries[child + 1], entries[child])) ++child;
          if (!better(entries[child], entries[at])) return;
          exchange(at, child); at = child;
        }
      }
      void erase(Id id) {
        const auto at = state.links[id].position[kind];
        if (at == none) return;
        exchange(at, entries.size() - 1);
        entries.pop_back();
        state.links[id].position[kind] = none;
        if (at < entries.size()) down(up(at));
      }
      void insert(Id id) {
        state.links[id].position[kind] = static_cast<Id>(entries.size());
        entries.push_back(id);
        up(entries.size() - 1);
      }
      void build() {
        for (std::size_t i = 0; i < entries.size(); ++i)
          state.links[entries[i]].position[kind] = static_cast<Id>(i);
        for (std::size_t i = entries.size() / 2; i; --i) down(i - 1);
      }
      Id top() const { return entries.empty() ? none : entries.front(); }
      std::array<Id, 3> first_three() const {
        std::array<Id, 3> result{{none, none, none}};
        if (entries.empty()) return result;
        std::array<std::size_t, 4> frontier{{0, 0, 0, 0}};
        std::size_t count = 1;
        for (std::size_t i = 0; i < result.size() && count; ++i) {
          std::size_t best = 0;
          for (std::size_t j = 1; j < count; ++j)
            if (better(entries[frontier[j]], entries[frontier[best]])) best = j;
          const auto at = frontier[best];
          result[i] = entries[at];
          frontier[best] = frontier[--count];
          if (i + 1 == result.size()) break;
          if (2 * at + 1 < entries.size()) frontier[count++] = 2 * at + 1;
          if (2 * at + 2 < entries.size()) frontier[count++] = 2 * at + 2;
        }
        return result;
      }
    };
    BasicPaperSplineSketch& sketch;
    const Snapshot* before;
    std::vector<Links> links;
    std::vector<Id> dirty;
    Heap masses{*this, mass_heap, {}}, splits{*this, split_heap, {}}, joins{*this, join_heap, {}};
    Id head = 0, tail, free = none;
    std::size_t size;
    std::size_t protected_count = 0;
    bool growing;

    HeapRebalance(BasicPaperSplineSketch& s, const Snapshot& snapshot, bool grow)
        : sketch(s), before(&snapshot), tail(static_cast<Id>(s.nodes_.size() - 1)),
          size(s.nodes_.size()), growing(grow && size < s.capacity_) {
      const auto slots = std::max(size, s.capacity_) + 2;
      links.reserve(slots);
      links.resize(size);
      dirty.reserve(slots);
      masses.entries.reserve(slots); splits.entries.reserve(slots); joins.entries.reserve(slots);
      sketch.nodes_.reserve(slots);
      for (std::size_t i = 0; i < size; ++i) {
        protected_count += static_cast<std::size_t>(sketch.nodes_[i].protected_threshold);
        links[i].prev = i ? static_cast<Id>(i - 1) : none;
        links[i].next = i + 1 < size ? static_cast<Id>(i + 1) : none;
      }
      for (Id id = head; id != none; id = links[id].next) cache(id);
      rebuild_heaps();
      verify();
    }
    HeapRebalance(const HeapRebalance&) = delete;
    HeapRebalance& operator=(const HeapRebalance&) = delete;
    bool active(Id id) const { return id != none && links[id].prev != id; }
    bool splittable(Id id) const {
      const auto left = links[id].prev;
      if (left == none) return false;
      const auto mid = midpoint(sketch.nodes_[left].x, sketch.nodes_[id].x);
      return mid > sketch.nodes_[left].x && mid < sketch.nodes_[id].x;
    }
    long double score(Id id, bool joined) const {
      const auto left = links[id].prev;
      const auto right = joined ? links[id].next : id;
      const auto previous = links[left].prev, next = links[right].next;
      const auto& nodes = sketch.nodes_;
      const long double h = static_cast<long double>(nodes[right].x) - nodes[left].x;
      const long double hp = previous != none ? static_cast<long double>(nodes[left].x) - nodes[previous].x : h;
      const long double hn = next != none ? static_cast<long double>(nodes[next].x) - nodes[right].x : h;
      const long double mass = nodes[id].mass + (joined ? nodes[right].mass : 0);
      const long double following = next != none ? nodes[next].mass : 0;
      const long double value = std::max(std::fabs(mass / h - nodes[left].mass / hp) / (h + hp),
          std::fabs(following / hn - mass / h) / (h + hn)) * h * h;
      return std::isfinite(value) ? value : std::numeric_limits<long double>::max();
    }
    bool heuristic_eligible(Id id) const {
      return splittable(id) && sketch.nodes_[id].mass >
          (growing ? 0 : 2 * static_cast<long double>(sketch.count_) / sketch.capacity_);
    }
    bool join_eligible(Id id) const {
      return links[id].prev != none && links[id].next != none &&
          !sketch.nodes_[id].protected_threshold && sketch.nodes_[id].mass +
          sketch.nodes_[links[id].next].mass <= 0.75L * sketch.bound();
    }
    void cache(Id id) {
      auto& node = sketch.nodes_[id];
      node.prefix = links[id].prev == none ? 0 : score(id, false);
      node.slope = links[id].prev == none || links[id].next == none ? 0 : score(id, true);
    }
    void rebuild_heaps() {
      masses.entries.clear(); splits.entries.clear(); joins.entries.clear();
      for (Id id = head; id != none; id = links[id].next) {
        links[id].position = {{none, none, none}};
        if (splittable(id) && sketch.nodes_[id].mass > 0) masses.entries.push_back(id);
        if (heuristic_eligible(id)) splits.entries.push_back(id);
        if (join_eligible(id)) joins.entries.push_back(id);
      }
      masses.build(); splits.build(); joins.build();
    }
    void refresh(Id id) {
      if (!active(id)) return;
      cache(id);
      if (splittable(id) && sketch.nodes_[id].mass > 0) masses.insert(id);
      if (heuristic_eligible(id)) splits.insert(id);
      if (join_eligible(id)) joins.insert(id);
    }
    std::array<Id, 7> detach_neighbours(Id center) {
      std::array<Id, 7> affected{{none, none, none, none, none, none, none}};
      auto at = center;
      for (unsigned i = 0; i < 2 && links[at].prev != none; ++i) at = links[at].prev;
      for (auto& id : affected) {
        if (at == none) break;
        id = at; at = links[at].next;
        masses.erase(id); splits.erase(id); joins.erase(id);
      }
      return affected;
    }
    void join(Id id) {
      if (!join_eligible(id)) throw std::logic_error("illegal paper heap join");
      const auto affected = detach_neighbours(id);
      const auto left = links[id].prev, right = links[id].next;
      sketch.nodes_[right].mass += sketch.nodes_[id].mass;
      links[left].next = right; links[right].prev = left;
      links[id].prev = id; links[id].next = free; free = id;
      --size;
      for (Id changed : affected) refresh(changed);
    }
    void split(Id right) {
      if (!splittable(right)) throw std::logic_error("unsplittable paper heap bucket");
      const auto left = links[right].prev;
      const auto mid = midpoint(sketch.nodes_[left].x, sketch.nodes_[right].x);
      // Ordered summation deliberately matches the scanning implementation.
      // An associative tree would change floating-point inverse boundaries.
      long double low = 0, high = 0;
      for (Id at = head; at != none; at = links[at].next) {
        high += sketch.nodes_[at].mass;
        if (at == left) low = high;
        if (at == right) break;
      }
      const auto value = std::clamp(before->rank(mid), low, high);
      const auto affected = detach_neighbours(right);
      Id id;
      if (free != none) { id = free; free = links[id].next; }
      else {
        id = static_cast<Id>(links.size());
        links.emplace_back(); sketch.nodes_.emplace_back(mid);
      }
      sketch.nodes_[id] = Node{mid, value - low, 0, 0, true};
      if constexpr (TrackBounds) {
        const auto interval = before->bounds(mid);
        sketch.nodes_[id].lower = interval.lower; sketch.nodes_[id].upper = interval.upper;
        sketch.nodes_[id].atom = before->atom(mid);
      }
      sketch.nodes_[right].mass = high - value;
      protected_count += split_protections(sketch.nodes_[left], sketch.nodes_[right]);
      sketch.nodes_[left].protected_threshold = sketch.nodes_[right].protected_threshold = true;
      links[id] = Links{}; links[id].prev = left; links[id].next = right;
      links[left].next = id; links[right].prev = id;
      ++size;
      for (Id changed : affected) refresh(changed);
      refresh(id);
      if (growing && size == sketch.capacity_) { growing = false; rebuild_heaps(); }
    }
    bool overlaps(Id join_id, Id split_id) const {
      return join_id == split_id || links[join_id].next == split_id;
    }
    Id compatible(const std::array<Id, 3>& candidates, Id split_id) const {
      for (Id id : candidates) if (id != none && !overlaps(id, split_id)) return id;
      return none;
    }
    void raise(bool requires_removal) {
      if (sketch.policy_ == BoundPolicy::theoretical)
        throw std::logic_error("paper baseline cannot find a legal join");
      bool unprotected = false;
      for (Id id = head; id != none; id = links[id].next)
        if (links[id].prev != none && links[id].next != none)
          unprotected |= !sketch.nodes_[id].protected_threshold;
      if (requires_removal && !unprotected)
        throw std::logic_error("paper baseline has no unprotected join threshold");
      sketch.factor_ *= 2;
      rebuild_heaps();
    }
    void consider(Id split_id, Id join_id, long double& gain, Id& split, Id& joined) const {
      if (join_id == none) return;
      const auto value = sketch.nodes_[split_id].prefix - 1.5L * sketch.nodes_[join_id].slope;
      if (value > gain || (value == gain && value > 0 &&
          (split == none || sketch.nodes_[split_id].x < sketch.nodes_[split].x))) {
        gain = value; split = split_id; joined = join_id;
      }
    }
    void tied_gains(std::size_t at, Id cheapest, long double cost, long double& gain,
                    Id& split, Id& joined) const {
      if (at >= splits.entries.size()) return;
      const auto id = splits.entries[at];
      const auto value = sketch.nodes_[id].prefix - cost;
      if (value < gain) return; // all descendants have no larger heuristic
      if (!overlaps(cheapest, id)) consider(id, cheapest, gain, split, joined);
      tied_gains(2 * at + 1, cheapest, cost, gain, split, joined);
      tied_gains(2 * at + 2, cheapest, cost, gain, split, joined);
    }
    void best_gain(const std::array<Id, 3>& candidates, Id& split, Id& joined) const {
      const auto cheapest = candidates[0];
      long double gain = 0;
      for (Id id : splits.first_three())
        if (id != none && !overlaps(cheapest, id)) {
          consider(id, cheapest, gain, split, joined);
          break;
        }
      if (gain > 0) tied_gains(0, cheapest, 1.5L * sketch.nodes_[cheapest].slope, gain, split, joined);
      // Only the two splits touching the cheapest join need a different join.
      for (Id id : {cheapest, links[cheapest].next})
        if (links[id].position[split_heap] != none) consider(id, compatible(candidates, id), gain, split, joined);
    }
    void mark(Id id) {
      if (!active(id) || links[id].dirty) return;
      masses.erase(id); splits.erase(id); joins.erase(id);
      links[id].dirty = true; dirty.push_back(id);
    }
    void mark_neighbours(Id center) {
      auto at = center;
      for (unsigned i = 0; i < 2 && links[at].prev != none; ++i) at = links[at].prev;
      for (unsigned i = 0; i < 7 && at != none; ++i) {
        mark(at); at = links[at].next;
      }
    }
    Id new_slot(double x) {
      Id id;
      if (free != none) { id = free; free = links[id].next; }
      else { id = static_cast<Id>(links.size()); links.emplace_back(); sketch.nodes_.emplace_back(x); }
      sketch.nodes_[id] = Node{x}; links[id] = Links{};
      return id;
    }
    void reestimate(const Snapshot& snapshot) {
      before = &snapshot;
      if (snapshot.items.front().first < sketch.nodes_[head].x) {
        mark_neighbours(head);
        const auto id = new_slot(snapshot.items.front().first);
        links[id].next = head; links[head].prev = id; head = id; ++size;
        mark(id);
      }
      if (snapshot.items.back().first > sketch.nodes_[tail].x) {
        mark_neighbours(tail);
        const auto id = new_slot(snapshot.items.back().first);
        links[id].prev = tail; links[tail].next = id; tail = id; ++size;
        mark(id);
      }
      std::size_t old = 0, pending = 0;
      long double previous = 0;
      for (Id id = head; id != none; id = links[id].next) {
        const auto x = sketch.nodes_[id].x;
        while (old < snapshot.nodes.size() && snapshot.nodes[old].x < x) ++old;
        while (pending < snapshot.items.size() && snapshot.items[pending].first <= x) ++pending;
        const long double old_rank = old < snapshot.nodes.size() && snapshot.nodes[old].x == x
            ? snapshot.nodes[old].prefix : spline_rank(snapshot.nodes, x);
        const long double value = old_rank + (pending ? snapshot.sums[pending - 1] : 0);
        const long double mass = std::max(0.0L, value - previous);
        set_bounds(sketch.nodes_[id], snapshot, old, pending);
        if (mass != sketch.nodes_[id].mass) {
          mark_neighbours(id);
          sketch.nodes_[id].mass = mass;
        }
        previous = value;
      }
      for (Id id : dirty) { links[id].dirty = false; refresh(id); }
      dirty.clear();
      verify();
    }
    std::vector<Node> snapshot_nodes() const {
      std::vector<Node> result;
      result.reserve(size);
      for (Id id = head; id != none; id = links[id].next) {
        const auto& node = sketch.nodes_[id];
        result.emplace_back(node.x, node.mass, 0, 0, node.protected_threshold);
        if constexpr (TrackBounds) {
          result.back().lower = node.lower; result.back().upper = node.upper;
          result.back().atom = node.atom;
        }
      }
      rebuild(result);
      return result;
    }
    void finish() {
      auto& order = masses.entries;
      order.clear();
      for (Id id = head; id != none; id = links[id].next) order.push_back(id);
      // Reuse links as a permutation map; no second Node array is allocated.
      for (std::size_t i = 0; i < links.size(); ++i) {
        links[i].position[0] = links[i].prev = static_cast<Id>(i);
      }
      for (std::size_t i = 0; i < order.size(); ++i) {
        const auto wanted = order[i], at = links[wanted].position[0], occupant = links[i].prev;
        std::swap(sketch.nodes_[i], sketch.nodes_[at]);
        std::swap(links[i].prev, links[at].prev);
        links[wanted].position[0] = static_cast<Id>(i);
        links[occupant].position[0] = at;
      }
      sketch.nodes_.resize(size);
      rebuild(sketch.nodes_);
    }
    void verify() const {
#ifdef SPLINESKETCH_VERIFY_PAPER_HEAPS
      const auto require = [](bool valid) {
        if (!valid) throw std::logic_error("paper heap invariant failed");
      };
      std::size_t live = 0;
      std::size_t protected_live = 0;
      Id previous = none;
      for (Id id = head; id != none; id = links[id].next) {
        require(id < links.size() && active(id) && links[id].prev == previous && !links[id].dirty);
        require(++live <= links.size());
        const auto& node = sketch.nodes_[id];
        protected_live += static_cast<std::size_t>(node.protected_threshold);
        require(std::isfinite(node.mass) && node.mass >= 0);
        if (previous != none) require(sketch.nodes_[previous].x < node.x);
        require(node.prefix == (previous == none ? 0 : score(id, false)));
        require(node.slope == (previous == none || links[id].next == none ? 0 : score(id, true)));
        require((links[id].position[mass_heap] != none) == (splittable(id) && node.mass > 0));
        require((links[id].position[split_heap] != none) == heuristic_eligible(id));
        require((links[id].position[join_heap] != none) == join_eligible(id));
        previous = id;
      }
      require(live == size && previous == tail);
      require(protected_live == protected_count);
      for (const auto* heap : {&masses, &splits, &joins}) {
        for (std::size_t i = 0; i < heap->entries.size(); ++i) {
          const auto id = heap->entries[i];
          require(id < links.size() && active(id) && links[id].position[heap->kind] == i);
          if (i) require(!heap->better(id, heap->entries[(i - 1) / 2]));
        }
      }
#endif
    }
    void run(bool materialize = true) {
      while (size > sketch.capacity_) {
        if (joins.top() == none) { raise(true); continue; }
        join(joins.top());
        verify();
      }
      for (;;) {
        auto split_id = masses.top();
        const bool mandatory = split_id != none && sketch.nodes_[split_id].mass > sketch.bound();
        Id join_id = none;
        const auto candidates = joins.first_three();
        if (!mandatory) {
          split_id = none;
          if (growing) split_id = splits.top();
          else if (size == sketch.capacity_ && joins.entries.size() >= (sketch.capacity_ + 2) / 3 + 2)
            best_gain(candidates, split_id, join_id);
        }
        if (split_id == none) break;
        if (!sketch.protection_budget_allows(protected_count,
                sketch.nodes_[links[split_id].prev], sketch.nodes_[split_id])) {
          // Mandatory splits can be avoided by raising the practical bound.
          // Optional refinement stops before consuming the extrema reserve.
          if (!mandatory) break;
          raise(false); continue;
        }
        if (size < sketch.capacity_) { split(split_id); verify(); continue; }
        if (join_id == none) join_id = compatible(candidates, split_id);
        if (join_id == none) { raise(false); continue; }
        join(join_id); split(split_id);
        verify();
      }
      verify();
      if (materialize) finish();
    }
  };
  void incorporate_heaps(const Snapshot& all) {
    const auto batch_limit = occurrence_batch_limit();
    Snapshot batch{nodes_, {}, {}};
    HeapRebalance state(*this, batch, false);
    std::size_t at = 0;
    std::uint64_t used = 0;
    bool first = true;
    while (at < all.items.size()) {
      if (!first) batch.nodes = state.snapshot_nodes();
      first = false;
      batch.items.clear(); batch.sums.clear();
      std::uint64_t size = 0;
      while (at < all.items.size() && size < batch_limit) {
        const auto take = std::min(all.items[at].second - used, batch_limit - size);
        batch.items.emplace_back(all.items[at].first, take);
        size += take; used += take;
        if (used == all.items[at].second) { ++at; used = 0; }
      }
      batch.prepare();
      state.reestimate(batch);
      state.run(false);
    }
    state.finish();
  }
  bool use_heap_backend() const {
    if (std::max(capacity_, nodes_.size()) >= HeapRebalance::none - 2) return false;
#ifdef SPLINESKETCH_PAPER_FORCE_HEAPS
    return true;
#else
    // At smaller capacities the simpler scans beat heap bookkeeping locally.
    return capacity_ >= 512 || nodes_.size() >= 512;
#endif
  }
  void rebalance(const Snapshot& snapshot, bool supplied, bool growing = false) {
    if (nodes_.empty()) return;
    const Snapshot local{supplied ? std::vector<Node>{} : nodes_, {}, {}};
    const Snapshot& before = supplied ? snapshot : local;
    if (use_heap_backend()) {
      HeapRebalance(*this, before, growing).run();
      return;
    }
    // Merge/resize/new extrema reduction: every join still obeys Definition 1.
    while (nodes_.size() > capacity_) {
      const auto joins = join_candidates();
      if (!joins.count) { raise_bound(true); continue; }
      join_at(joins.best.front().second);
    }
    std::size_t protected_count = 0;
    for (const auto& node : nodes_)
      protected_count += static_cast<std::size_t>(node.protected_threshold);
    for (;;) {
      std::size_t split = npos;
      long double score = -1;
      bool mandatory = false;
      for (std::size_t i = 1; i < nodes_.size(); ++i) {
        if (nodes_[i].mass <= bound() || !can_split(i)) continue;
        if (nodes_[i].mass > score) { split = i; score = nodes_[i].mass; }
      }
      mandatory = split != npos;
      const auto joins = nodes_.size() >= capacity_ ? join_candidates() :
          JoinCandidates{};
      std::size_t join = npos;
      if (!mandatory) {
        if (growing && nodes_.size() < capacity_) {
          for (std::size_t i = 1; i < nodes_.size(); ++i)
            if (can_split(i) && nodes_[i].mass > 0 && heuristic(i) > score) {
              split = i; score = heuristic(i);
            }
        } else if (nodes_.size() == capacity_ && joins.count >= (capacity_ + 2) / 3 + 2) {
          long double gain = 0;
          for (std::size_t i = 1; i < nodes_.size(); ++i) {
            if (nodes_[i].mass <= 2 * static_cast<long double>(count_) / capacity_ || !can_split(i)) continue;
            for (const auto& candidate : joins.best) {
              if (candidate.second == npos) break;
              if (candidate.second == i || candidate.second + 1 == i) continue;
              const long double value = heuristic(i) - 1.5L * candidate.first;
              if (value > gain) { gain = value; split = i; join = candidate.second; }
              break;
            }
          }
        }
      }
      if (split == npos) break;
      if (!protection_budget_allows(protected_count, nodes_[split - 1], nodes_[split])) {
        if (!mandatory) break;
        raise_bound(false); continue;
      }
      const auto protections = split_protections(nodes_[split - 1], nodes_[split]);
      if (nodes_.size() < capacity_) {
        split_at(split, before); protected_count += protections; continue;
      }
      if (join == npos) join = joins.nonoverlapping(split);
      if (join == npos) { raise_bound(false); continue; }
      // Join first so vector capacity never needs more than k slots here.
      join_at(join);
      if (join < split) --split;
      split_at(split, before);
      protected_count += protections;
    }
    // Edits read masses/prefixes and the frozen snapshot, never live slopes.
    // Restore slopes before another batch snapshot or public query uses them.
    rebuild(nodes_);
  }
};

using PaperSplineSketch = BasicPaperSplineSketch<false>;
using CertifiedPaperSplineSketch = BasicPaperSplineSketch<true>;

static_assert(std::is_nothrow_move_assignable<CertifiedPaperSplineSketch>::value,
              "certified paper sketch commits must not throw");
static_assert(std::is_nothrow_move_assignable<PaperSplineSketch>::value,
              "PaperSplineSketch commits must not throw");

} // namespace splinesketch
