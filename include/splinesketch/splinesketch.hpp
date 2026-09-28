#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace splinesketch {

// A buffered SplineSketch for finite doubles. rank(x) estimates the number of
// inserted values <= x. Heavy hitters are held exactly in a Misra-Gries table.
// The implementation follows Sections 3 and 4 of arXiv:2504.01206v3.
class SplineSketch {
  static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                std::numeric_limits<double>::is_iec559,
                "SplineSketch requires IEEE-754 binary64 doubles");
 public:
  explicit SplineSketch(std::size_t buckets = 128)
      : capacity_(buckets), epoch_end_(4 * buckets) {
    if (buckets < 6 || buckets > std::numeric_limits<std::size_t>::max() / 4)
      throw std::invalid_argument("bucket count must be at least 6");
    nodes_.reserve(buckets + 2);
    pending_.reserve(buckets + 2);
    heavy_.reserve(buckets);
  }

  void add(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("value must be finite");
    if (count_ == std::numeric_limits<std::uint64_t>::max())
      throw std::overflow_error("observation count overflow");
    value = canonical(value);
    ++count_;
    auto found = heavy_.find(value);
    if (found != heavy_.end()) {
      ++found->second.residual;
      ++found->second.exact;
      return;
    }
    if (heavy_.size() < capacity_ - 1) {
      heavy_.emplace(value, Heavy{1, 1});
      return;
    }
    // A full MG table consumes one unit from each counter and the new item.
    for (auto it = heavy_.begin(); it != heavy_.end();) {
      if (--it->second.residual == 0) {
        forward(it->first, it->second.exact);
        it = heavy_.erase(it);
      } else {
        ++it;
      }
    }
    forward(value, 1);
    if (pending_.size() >= capacity_) consolidate();
  }

  std::uint64_t count() const noexcept { return count_; }
  std::size_t bucket_capacity() const noexcept { return capacity_; }
  std::size_t bucket_count() const noexcept { return nodes_.size(); }
  std::size_t heavy_hitter_count() const noexcept { return heavy_.size(); }

  double rank(double value) const {
    if (std::isnan(value)) throw std::invalid_argument("query must not be NaN");
    if (value == -std::numeric_limits<double>::infinity()) return 0;
    if (value == std::numeric_limits<double>::infinity())
      return static_cast<double>(count_);
    long double result = spline_rank(nodes_, value);
    for (const auto& item : pending_)
      if (item.first <= value) result += item.second;
    for (const auto& item : heavy_)
      if (item.first <= value) result += item.second.exact;
    return static_cast<double>(std::clamp(result, 0.0L,
                                         static_cast<long double>(count_)));
  }

  // Inverse of the estimated CDF. The answer can lie between input values.
  double quantile(double q) const {
    if (!(q >= 0 && q <= 1) || std::isnan(q))
      throw std::invalid_argument("q must be in [0, 1]");
    if (count_ == 0) throw std::logic_error("quantile of empty sketch");
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    if (!nodes_.empty()) {
      lo = nodes_.front().x;
      hi = nodes_.back().x;
    }
    for (const auto& item : pending_) {
      lo = std::min(lo, item.first);
      hi = std::max(hi, item.first);
    }
    for (const auto& item : heavy_) {
      lo = std::min(lo, item.first);
      hi = std::max(hi, item.first);
    }
    if (q == 0) return lo;
    if (q == 1) return hi;
    if (lo == hi) return lo;
    const long double target = std::ceil(q * static_cast<long double>(count_));
    // Preparing exact-value prefix counts pays off across the inverse search,
    // but keep small tables allocation-free. All scratch storage is local so
    // const queries do not mutate the sketch or introduce a shared cache.
    const std::size_t exact_size = pending_.size() + heavy_.size();
    const bool prepare_exact = exact_size > 16;
    Snapshot exact{{}, {}, {}};
    if (prepare_exact) {
      exact.pending.reserve(exact_size);
      exact.pending.insert(exact.pending.end(), pending_.begin(), pending_.end());
      for (const auto& item : heavy_)
        exact.pending.emplace_back(item.first, item.second.exact);
      exact.prepare();
    }
    // Grouping integer counts changes floating-point addition order. Near the
    // target, use rank() itself to preserve its rounding and inverse semantics.
    const long double rounding_slack = static_cast<long double>(count_) *
        (4 * (static_cast<long double>(exact_size) + 1) *
             std::numeric_limits<long double>::epsilon() +
         std::numeric_limits<double>::epsilon());
    const auto reaches_target = [&](double value) {
      if (!prepare_exact) return rank(value) >= target;
      const long double estimate = spline_rank(nodes_, value) + exact.rank(value);
      if (std::fabs(estimate - target) <= rounding_slack)
        return rank(value) >= target;
      return static_cast<double>(std::clamp(estimate, 0.0L,
          static_cast<long double>(count_))) >= target;
    };
    // Binary search the ordered IEEE-754 bit pattern. This also works when
    // the endpoints are 1e308 apart or only one ULP apart.
    std::uint64_t low = ordered_bits(lo), high = ordered_bits(hi);
    while (low < high) {
      const std::uint64_t middle = low + (high - low) / 2;
      if (reaches_target(from_ordered_bits(middle))) high = middle;
      else low = middle + 1;
    }
    return from_ordered_bits(low);
  }

  void consolidate() {
    if (pending_.empty()) return;
    advance_epoch();
    Snapshot before{nodes_, pending_, {}};
    before.prepare();
    if (nodes_.empty()) {
      initialize(before);
    } else {
      const double min_value = before.pending.front().first;
      const double max_value = before.pending.back().first;
      if (min_value < nodes_.front().x)
        nodes_.insert(nodes_.begin(), Node{min_value});
      if (max_value > nodes_.back().x)
        nodes_.push_back(Node{max_value});
      reestimate(before);
    }
    pending_.clear();
    reduce_to_capacity();
    rebalance(before);
  }

  void resize(std::size_t new_capacity) {
    if (new_capacity < 6 || new_capacity > std::numeric_limits<std::size_t>::max() / 4)
      throw std::invalid_argument("bucket count must be at least 6");
    consolidate();
    if (new_capacity == capacity_) return;
    if (new_capacity > capacity_ + capacity_ / 4 ||
        capacity_ > new_capacity + new_capacity / 4)
      for (auto& node : nodes_) node.protected_threshold = false;
    capacity_ = new_capacity;
    if (nodes_.capacity() < capacity_ + 2) nodes_.reserve(capacity_ + 2);
    pending_.reserve(capacity_ + 2);
    heavy_.reserve(capacity_);
    shrink_heavy();
    consolidate();
    reduce_to_capacity();
    if (nodes_.size() < capacity_) {
      Snapshot before{nodes_, {}, {}};
      while (nodes_.size() < capacity_) {
        const auto split = best_split(false);
        if (split == npos) break;
        split_at(split, before);
      }
    }
    const auto room = std::numeric_limits<std::uint64_t>::max() - count_;
    const auto next_epoch = count_ + std::min<std::uint64_t>(room, count_ / 4 + 1);
    epoch_end_ = std::max<std::uint64_t>(next_epoch,
                                          static_cast<std::uint64_t>(4 * capacity_));
    rebuild(nodes_);
  }

  void merge(const SplineSketch& other) {
    if (this == &other) {
      SplineSketch copy(other);
      merge(copy);
      return;
    }
    if (std::numeric_limits<std::uint64_t>::max() - count_ < other.count_)
      throw std::overflow_error("observation count overflow");
    SplineSketch a(*this), b(other);
    a.consolidate();
    b.consolidate();
    SplineSketch result(capacity_);
    result.count_ = a.count_ + b.count_;
    const SplineSketch& larger = a.count_ >= b.count_ ? a : b;
    result.epoch_end_ = larger.epoch_end_;
    result.bucket_bound_factor_ = std::max(a.bucket_bound_factor_, b.bucket_bound_factor_);
    for (const auto& node : a.nodes_) result.nodes_.push_back(Node{node.x});
    for (const auto& node : b.nodes_) result.nodes_.push_back(Node{node.x});
    std::sort(result.nodes_.begin(), result.nodes_.end(),
              [](const Node& x, const Node& y) { return x.x < y.x; });
    result.nodes_.erase(std::unique(result.nodes_.begin(), result.nodes_.end(),
        [](const Node& x, const Node& y) { return x.x == y.x; }), result.nodes_.end());
    for (auto& node : result.nodes_) {
      const auto it = std::lower_bound(larger.nodes_.begin(), larger.nodes_.end(), node.x,
          [](const Node& candidate, double x) { return candidate.x < x; });
      if (it != larger.nodes_.end() && it->x == node.x)
        node.protected_threshold = it->protected_threshold;
    }
    result.advance_epoch();
    long double previous = 0;
    for (auto& node : result.nodes_) {
      const long double current = spline_rank(a.nodes_, node.x) +
                                  spline_rank(b.nodes_, node.x);
      node.mass = std::max(0.0L, current - previous);
      previous = current;
    }
    rebuild(result.nodes_);
    for (const auto& item : a.heavy_)
      result.insert_summary(item.first, item.second);
    for (const auto& item : b.heavy_)
      result.insert_summary(item.first, item.second);
    result.reduce_to_capacity();
    result.consolidate();
    *this = std::move(result);
  }

 private:
  struct Node {
    double x = 0;
    long double mass = 0;
    long double prefix = 0;
    long double slope = 0;
    bool protected_threshold = false;
  };
  struct Heavy {
    std::uint64_t residual = 0;
    std::uint64_t exact = 0;
  };
  struct Snapshot {
    std::vector<Node> nodes;
    std::vector<std::pair<double, std::uint64_t>> pending;
    std::vector<long double> sums;

    void prepare() {
      std::sort(pending.begin(), pending.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
      std::size_t out = 0;
      for (const auto& item : pending) {
        if (out && pending[out - 1].first == item.first)
          pending[out - 1].second += item.second;
        else pending[out++] = item;
      }
      pending.resize(out);
      sums.reserve(out);
      long double sum = 0;
      for (const auto& item : pending) {
        sum += item.second;
        sums.push_back(sum);
      }
    }
    long double rank(double x) const {
      const auto it = std::upper_bound(pending.begin(), pending.end(), x,
          [](double x, const auto& item) { return x < item.first; });
      const std::size_t at = static_cast<std::size_t>(it - pending.begin());
      return SplineSketch::spline_rank(nodes, x) + (at ? sums[at - 1] : 0);
    }
  };

  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
  std::size_t capacity_;
  std::uint64_t count_ = 0;
  std::uint64_t epoch_end_;
  long double bucket_bound_factor_ = 3;
  std::vector<Node> nodes_;
  std::vector<std::pair<double, std::uint64_t>> pending_;
  std::unordered_map<double, Heavy> heavy_;

  static double canonical(double x) { return x == 0 ? 0 : x; }
  static std::uint64_t ordered_bits(double x) {
    std::uint64_t bits;
    std::memcpy(&bits, &x, sizeof bits);
    return bits & (std::uint64_t{1} << 63) ? ~bits : bits ^ (std::uint64_t{1} << 63);
  }
  static double from_ordered_bits(std::uint64_t bits) {
    bits = bits & (std::uint64_t{1} << 63) ? bits ^ (std::uint64_t{1} << 63) : ~bits;
    double x;
    std::memcpy(&x, &bits, sizeof x);
    return x;
  }
  static double midpoint(double a, double b) {
    return static_cast<double>((static_cast<long double>(a) + b) / 2);
  }
  static long double span(double a, double b) {
    return static_cast<long double>(b) - a;
  }
  static long double endpoint_slope(long double h0, long double h1,
                                    long double d0, long double d1) {
    long double slope = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
    if (slope <= 0) return 0;
    if (d1 == 0 || slope > 3 * d0) slope = 3 * d0;
    return slope;
  }
  static void rebuild(std::vector<Node>& nodes) {
    long double sum = 0;
    for (auto& node : nodes) {
      sum += node.mass;
      node.prefix = sum;
      node.slope = 0;
    }
    const auto size = nodes.size();
    if (size < 2) return;
    long double previous_h = span(nodes[0].x, nodes[1].x);
    long double previous_delta = nodes[1].mass / previous_h;
    if (size == 2) {
      nodes[0].slope = nodes[1].slope = previous_delta;
      return;
    }
    const long double second_h = span(nodes[1].x, nodes[2].x);
    const long double second_delta = nodes[2].mass / second_h;
    nodes[0].slope = endpoint_slope(previous_h, second_h,
                                    previous_delta, second_delta);
    for (std::size_t i = 1; i + 1 < size; ++i) {
      const long double next_h = span(nodes[i].x, nodes[i + 1].x);
      const long double next_delta = nodes[i + 1].mass / next_h;
      if (previous_delta > 0 && next_delta > 0) {
        const long double w1 = 2 * next_h + previous_h;
        const long double w2 = next_h + 2 * previous_h;
        nodes[i].slope = (w1 + w2) /
            (w1 / previous_delta + w2 / next_delta);
      }
      previous_h = next_h;
      previous_delta = next_delta;
    }
    const long double penultimate_h = span(nodes[size - 3].x, nodes[size - 2].x);
    const long double penultimate_delta = nodes[size - 2].mass / penultimate_h;
    nodes.back().slope = endpoint_slope(previous_h, penultimate_h,
                                        previous_delta, penultimate_delta);
  }
  // Structural edits need current prefixes for subsequent splits. Their
  // callers rebuild slopes once the batch is complete, before any snapshot
  // or public query can use the edited spline.
  void rebuild_prefixes() {
    long double sum = 0;
    for (auto& node : nodes_) {
      sum += node.mass;
      node.prefix = sum;
    }
  }
  static long double spline_rank(const std::vector<Node>& nodes, double x) {
    if (nodes.empty() || x < nodes.front().x) return 0;
    if (x >= nodes.back().x) return nodes.back().prefix;
    const auto it = std::lower_bound(nodes.begin(), nodes.end(), x,
        [](const Node& node, double x) { return node.x < x; });
    if (it->x == x) return it->prefix;
    const Node& right = *it;
    const Node& left = *(it - 1);
    const long double h = span(left.x, right.x);
    const long double t = span(left.x, x) / h;
    const long double t2 = t * t, t3 = t2 * t;
    const long double value = (2 * t3 - 3 * t2 + 1) * left.prefix +
        (t3 - 2 * t2 + t) * h * left.slope +
        (-2 * t3 + 3 * t2) * right.prefix +
        (t3 - t2) * h * right.slope;
    return std::clamp(value, left.prefix, right.prefix);
  }
  void forward(double x, std::uint64_t weight) {
    if (weight) pending_.emplace_back(x, weight);
  }
  void advance_epoch() {
    if (count_ <= epoch_end_) return;
    for (auto& node : nodes_) node.protected_threshold = false;
    bucket_bound_factor_ = 3;
    while (epoch_end_ < count_) {
      const auto increment = std::max<std::uint64_t>(1, epoch_end_ / 4);
      if (epoch_end_ > std::numeric_limits<std::uint64_t>::max() - increment) {
        epoch_end_ = std::numeric_limits<std::uint64_t>::max();
        break;
      }
      epoch_end_ += increment;
    }
  }
  void initialize(const Snapshot& before) {
    const auto& items = before.pending;
    const std::size_t desired = std::min(capacity_, items.size());
    if (desired == 0) return;
    if (items.size() <= capacity_) {
      for (const auto& item : items) nodes_.push_back(Node{item.first, static_cast<long double>(item.second)});
    } else {
      const long double total = before.sums.back();
      for (std::size_t i = 0; i < desired; ++i) {
        const long double target = std::ceil(static_cast<long double>(i) * (total - 1) /
                                             (desired - 1));
        const auto it = std::upper_bound(before.sums.begin(), before.sums.end(), target);
        const auto index = std::min<std::size_t>(items.size() - 1, it - before.sums.begin());
        nodes_.push_back(Node{items[index].first});
      }
      std::sort(nodes_.begin(), nodes_.end(),
                [](const Node& a, const Node& b) { return a.x < b.x; });
      nodes_.erase(std::unique(nodes_.begin(), nodes_.end(),
          [](const Node& a, const Node& b) { return a.x == b.x; }), nodes_.end());
      // Fill duplicate quantile selections with distinct input values.
      for (const auto& item : items) {
        if (nodes_.size() == desired) break;
        const auto it = std::lower_bound(nodes_.begin(), nodes_.end(), item.first,
            [](const Node& node, double x) { return node.x < x; });
        if (it == nodes_.end() || it->x != item.first) nodes_.insert(it, Node{item.first});
      }
      reestimate(before);
    }
    rebuild(nodes_);
  }
  void reestimate(const Snapshot& before) {
    long double previous = 0;
    std::size_t old_at = 0, pending_at = 0;
    for (auto& node : nodes_) {
      // All three sequences are sorted. Existing thresholds can reuse their
      // snapshot prefix directly; new endpoints use the ordinary spline path.
      while (old_at < before.nodes.size() && before.nodes[old_at].x < node.x)
        ++old_at;
      while (pending_at < before.pending.size() &&
             before.pending[pending_at].first <= node.x)
        ++pending_at;
      const long double old_rank =
          old_at < before.nodes.size() && before.nodes[old_at].x == node.x
              ? before.nodes[old_at].prefix : spline_rank(before.nodes, node.x);
      const long double current = old_rank + (pending_at ? before.sums[pending_at - 1] : 0);
      node.mass = std::max(0.0L, current - previous);
      previous = current;
    }
    rebuild(nodes_);
  }
  long double bound() const {
    return bucket_bound_factor_ * static_cast<long double>(count_) / capacity_;
  }
  long double heuristic(std::size_t i) const {
    if (i == 0 || i >= nodes_.size()) return 0;
    const long double h = span(nodes_[i - 1].x, nodes_[i].x);
    const long double density = nodes_[i].mass / h;
    const long double hp = i > 1 ? span(nodes_[i - 2].x, nodes_[i - 1].x) : h;
    const long double dp = nodes_[i - 1].mass / hp;
    const long double hn = i + 1 < nodes_.size() ?
        span(nodes_[i].x, nodes_[i + 1].x) : h;
    const long double dn = i + 1 < nodes_.size() ? nodes_[i + 1].mass / hn : 0;
    return std::max(std::fabs(density - dp) / (h + hp),
                    std::fabs(dn - density) / (h + hn)) * h * h;
  }
  long double joined_heuristic(std::size_t i) const {
    const long double h = span(nodes_[i - 1].x, nodes_[i + 1].x);
    const long double density = (nodes_[i].mass + nodes_[i + 1].mass) / h;
    const long double hp = i > 1 ? span(nodes_[i - 2].x, nodes_[i - 1].x) : h;
    const long double dp = nodes_[i - 1].mass / hp;
    const long double hn = i + 2 < nodes_.size() ?
        span(nodes_[i + 1].x, nodes_[i + 2].x) : h;
    const long double dn = i + 2 < nodes_.size() ? nodes_[i + 2].mass / hn : 0;
    return std::max(std::fabs(density - dp) / (h + hp),
                    std::fabs(dn - density) / (h + hn)) * h * h;
  }
  std::size_t best_join(std::size_t avoid = npos, bool enforce_bound = true,
                        bool enforce_protection = true) const {
    std::size_t best = npos;
    long double score = std::numeric_limits<long double>::infinity();
    for (std::size_t i = 1; i + 1 < nodes_.size(); ++i) {
      if (i == avoid || (avoid != npos && i + 1 == avoid)) continue;
      if (enforce_protection && nodes_[i].protected_threshold) continue;
      if (enforce_bound && nodes_[i].mass + nodes_[i + 1].mass > 0.75L * bound()) continue;
      const long double candidate = joined_heuristic(i);
      if (candidate < score) { score = candidate; best = i; }
    }
    return best;
  }
  void join_at(std::size_t i) {
    nodes_[i + 1].mass += nodes_[i].mass;
    nodes_.erase(nodes_.begin() + static_cast<std::ptrdiff_t>(i));
    rebuild_prefixes();
  }
  bool can_split(std::size_t i) const {
    if (i == 0 || i >= nodes_.size()) return false;
    const double mid = midpoint(nodes_[i - 1].x, nodes_[i].x);
    return mid > nodes_[i - 1].x && mid < nodes_[i].x;
  }
  void split_at(std::size_t i, const Snapshot& before) {
    const double mid = midpoint(nodes_[i - 1].x, nodes_[i].x);
    const long double left_rank = nodes_[i - 1].prefix;
    const long double right_rank = nodes_[i].prefix;
    const long double middle_rank = std::clamp(before.rank(mid), left_rank, right_rank);
    nodes_.insert(nodes_.begin() + static_cast<std::ptrdiff_t>(i),
                  Node{mid, middle_rank - left_rank});
    nodes_[i + 1].mass = right_rank - middle_rank;
    nodes_[i - 1].protected_threshold = true;
    nodes_[i].protected_threshold = true;
    nodes_[i + 1].protected_threshold = true;
    rebuild_prefixes();
  }
  std::size_t best_split(bool oversized_only) const {
    std::size_t best = npos;
    long double best_score = -1;
    for (std::size_t i = 1; i < nodes_.size(); ++i) {
      if (!can_split(i)) continue;
      if (oversized_only && nodes_[i].mass <= bound()) continue;
      const long double score = oversized_only ? nodes_[i].mass / bound() : heuristic(i);
      if (score > best_score) { best = i; best_score = score; }
    }
    return best;
  }
  void reduce_to_capacity() {
    const auto old_size = nodes_.size();
    while (nodes_.size() > capacity_) {
      auto join = best_join(npos, false, true);
      if (join == npos) join = best_join(npos, false, false);
      if (join == npos) break;
      join_at(join);
    }
    // resize() and merge() may take a snapshot immediately after reduction.
    if (nodes_.size() != old_size) rebuild(nodes_);
  }
  void rebalance(const Snapshot& before) {
    bool changed = false;
    // Each split protects three thresholds, so at most O(k) passes occur.
    for (std::size_t pass = 0; pass < capacity_ * 2; ++pass) {
      const auto split = best_split(true);
      if (split == npos) break;
      if (nodes_.size() < capacity_) {
        split_at(split, before);
        changed = true;
        continue;
      }
      auto join = best_join(split);
      if (join == npos) {
        bucket_bound_factor_ *= 2;
        continue;
      }
      const double left = nodes_[split - 1].x;
      const double right = nodes_[split].x;
      join_at(join);
      changed = true;
      const auto it = std::lower_bound(nodes_.begin(), nodes_.end(), right,
          [](const Node& node, double x) { return node.x < x; });
      const std::size_t new_split = static_cast<std::size_t>(it - nodes_.begin());
      if (new_split > 0 && nodes_[new_split - 1].x == left) split_at(new_split, before);
    }
    for (std::size_t pass = 0; pass < capacity_; ++pass) {
      if (nodes_.size() < 6) break;
      // A split bucket overlaps at most two join pairs. The cheapest three
      // joinable pairs therefore contain the best non-overlapping choice.
      std::array<std::pair<long double, std::size_t>, 3> cheapest{{
          {std::numeric_limits<long double>::infinity(), npos},
          {std::numeric_limits<long double>::infinity(), npos},
          {std::numeric_limits<long double>::infinity(), npos}}};
      std::size_t joinable = 0;
      for (std::size_t j = 1; j + 1 < nodes_.size(); ++j) {
        if (nodes_[j].protected_threshold ||
            nodes_[j].mass + nodes_[j + 1].mass > 0.75L * bound()) continue;
        ++joinable;
        const auto candidate = std::make_pair(joined_heuristic(j), j);
        if (candidate < cheapest[2]) {
          cheapest[2] = candidate;
          std::sort(cheapest.begin(), cheapest.end());
        }
      }
      if (joinable < capacity_ / 3 + 2) break;
      std::size_t best_split_index = npos, best_join_index = npos;
      long double largest_gain = 0;
      for (std::size_t i = 1; i < nodes_.size(); ++i) {
        if (!can_split(i) || nodes_[i].mass <= 2 * static_cast<long double>(count_) / capacity_) continue;
        auto candidate = cheapest.begin();
        while (candidate != cheapest.end() &&
               (candidate->second == i || candidate->second + 1 == i)) ++candidate;
        if (candidate == cheapest.end() || candidate->second == npos) continue;
        const long double gain = heuristic(i) - 1.5L * candidate->first;
        if (gain > largest_gain) {
          largest_gain = gain;
          best_split_index = i;
          best_join_index = candidate->second;
        }
      }
      if (best_split_index == npos) break;
      const double left = nodes_[best_split_index - 1].x;
      const double right = nodes_[best_split_index].x;
      join_at(best_join_index);
      changed = true;
      const auto it = std::lower_bound(nodes_.begin(), nodes_.end(), right,
          [](const Node& node, double x) { return node.x < x; });
      const std::size_t new_split = static_cast<std::size_t>(it - nodes_.begin());
      if (new_split > 0 && nodes_[new_split - 1].x == left) split_at(new_split, before);
    }
    if (changed) rebuild(nodes_);
  }
  void insert_summary(double x, Heavy entry) {
    auto found = heavy_.find(x);
    if (found != heavy_.end()) {
      found->second.residual += entry.residual;
      found->second.exact += entry.exact;
      return;
    }
    while (entry.residual && heavy_.size() >= capacity_ - 1) {
      std::uint64_t minimum = entry.residual;
      for (const auto& item : heavy_)
        minimum = std::min(minimum, item.second.residual);
      for (auto it = heavy_.begin(); it != heavy_.end();) {
        it->second.residual -= minimum;
        if (it->second.residual == 0) {
          forward(it->first, it->second.exact);
          it = heavy_.erase(it);
        } else ++it;
      }
      entry.residual -= minimum;
    }
    if (entry.residual) heavy_.emplace(x, entry);
    else forward(x, entry.exact);
  }
  void shrink_heavy() {
    while (heavy_.size() >= capacity_) {
      std::uint64_t minimum = std::numeric_limits<std::uint64_t>::max();
      for (const auto& item : heavy_)
        minimum = std::min(minimum, item.second.residual);
      for (auto it = heavy_.begin(); it != heavy_.end();) {
        it->second.residual -= minimum;
        if (it->second.residual == 0) {
          forward(it->first, it->second.exact);
          it = heavy_.erase(it);
        } else ++it;
      }
    }
  }
};

}  // namespace splinesketch
