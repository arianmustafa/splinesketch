#pragma once

#include <splinesketch/splinesketch.hpp>
#include <variant>

namespace splinesketch {
namespace detail {

// An integer rank envelope. Entries stand for retained, ordered occurrences;
// equal keys remain separate. See docs/certified-rank-proofs.md.
class RankEnvelope {
 public:
  struct Entry {
    double value;
    std::uint64_t delta, prefix;  // gap = prefix - previous prefix
  };
  struct Bounds { std::uint64_t lower, upper; };
  struct Insertion { std::size_t index; Entry entry; };

  explicit RankEnvelope(std::size_t capacity) : capacity_(capacity) {}
  std::uint64_t count() const noexcept { return count_; }
  std::size_t capacity() const noexcept { return capacity_; }
  std::size_t size() const noexcept {
    return std::visit([](const auto& entries) { return entries.size(); }, storage_);
  }
  std::uint64_t threshold() const noexcept {
    return std::max<std::uint64_t>(1, count_ / capacity_);
  }
  std::uint64_t uncertainty() const noexcept { return uncertainty_; }
  Entry entry(std::size_t index) const noexcept {
    return std::visit([&](const auto& entries) { return widen(entries[index]); }, storage_);
  }

  // Allocate before changing the spline. Promotion changes only representation;
  // a subsequent spline failure leaves observations and certificates unchanged.
  Insertion prepare_add(double value) {
    if (count_ >= std::numeric_limits<std::uint32_t>::max() && storage_.index() == 0)
      promote();
    return std::visit([&](auto& entries) -> Insertion {
      const std::size_t index = upper(entries, value);
      const std::uint64_t previous = index ? entries[index - 1].prefix : 0;
      const std::uint64_t prefix = previous + 1;
      const std::uint64_t delta = index && index < entries.size()
          ? entries[index].prefix - previous + entries[index].delta - 1 : 0;
      if (entries.size() == entries.capacity()) {
        const auto limit = entries.max_size();
        if (entries.size() == limit) throw std::length_error("rank envelope full");
        const auto spare = std::max<std::size_t>(8, entries.capacity() / 4);
        entries.reserve(entries.capacity() + std::min(limit - entries.capacity(), spare));
      }
      return {index, {value, delta, prefix}};
    }, storage_);
  }

  void commit_add(Insertion insertion) noexcept {
    std::visit([&](auto& entries) {
      using Stored = typename std::decay_t<decltype(entries)>::value_type;
      const Stored inserted = narrow<Stored>(insertion.entry);
      entries.push_back(inserted);
      // Move the suffix and shift its integer prefixes in one pass.
      for (std::size_t i = entries.size() - 1; i > insertion.index; --i) {
        auto moved = entries[i - 1];
        ++moved.prefix;
        entries[i] = moved;
      }
      entries[insertion.index] = inserted;
    }, storage_);
    ++count_;
    // Insertion changes only the pairs adjacent to the inserted tag (and
    // the old first tag when inserting a new minimum). A threshold increase
    // can affect all pairs. Skip a full scan only when no deletion is possible.
    const bool threshold_changed = count_ > capacity_ && count_ % capacity_ == 0;
    if (threshold_changed) compress();
    else compress_near(insertion.index);
  }

  Bounds bounds(double value) const noexcept {
    return std::visit([&](const auto& entries) {
      return before(entries, upper(entries, value));
    }, storage_);
  }
  void resize(std::size_t capacity) noexcept {
    capacity_ = std::min(capacity_, capacity);
    compress();
  }

  RankEnvelope merged(const RankEnvelope& other) const {
    if (std::numeric_limits<std::uint64_t>::max() - count_ < other.count_)
      throw std::overflow_error("observation count overflow");
    RankEnvelope result(std::min(capacity_, other.capacity_));
    result.count_ = count_ + other.count_;
    if (result.count_ > std::numeric_limits<std::uint32_t>::max())
      result.storage_.emplace<1>();
    std::visit([&](const auto& left, const auto& right, auto& output) {
      using Stored = typename std::decay_t<decltype(output)>::value_type;
      const auto limit = output.max_size();
      if (right.size() > limit || left.size() > limit - right.size())
        throw std::length_error("rank envelope full");
      output.reserve(left.size() + right.size());
      std::size_t a = 0, b = 0;
      while (a < left.size() || b < right.size()) {
        // Every A occurrence precedes B occurrences at equal values. The
        // union cursors locate B's strict rank and A's inclusive rank directly.
        const bool take_a = b == right.size() ||
            (a < left.size() && left[a].value <= right[b].value);
        const Entry tag = take_a ? widen(left[a++]) : widen(right[b++]);
        const Bounds extra = take_a ? other.before(right, b) : before(left, a);
        const std::uint64_t minimum = tag.prefix + extra.lower;
        const std::uint64_t maximum = tag.prefix + tag.delta + extra.upper;
        output.push_back(narrow<Stored>({tag.value, maximum - minimum, minimum}));
      }
    }, storage_, other.storage_, result.storage_);
    result.compress();
    return result;
  }

 private:
  struct SmallEntry { double value; std::uint32_t delta, prefix; };
  using Storage = std::variant<std::vector<SmallEntry>, std::vector<Entry>>;
  template<class Stored>
  static Entry widen(const Stored& entry) noexcept {
    return {entry.value, entry.delta, entry.prefix};
  }
  template<class Stored>
  static Stored narrow(Entry entry) noexcept {
    // The selected storage can represent all positions through count_.
    using Integer = decltype(Stored::prefix);
    return {entry.value, static_cast<Integer>(entry.delta),
                         static_cast<Integer>(entry.prefix)};
  }
  void promote() {
    const auto& small = std::get<0>(storage_);
    std::vector<Entry> wide;
    wide.reserve(small.capacity());
    for (const auto& entry : small) wide.push_back(widen(entry));
    storage_.emplace<1>(std::move(wide));
  }
  template<class Entries>
  static std::size_t upper(const Entries& entries, double value) noexcept {
    return static_cast<std::size_t>(std::upper_bound(entries.begin(), entries.end(), value,
        [](double key, const auto& entry) { return key < entry.value; }) - entries.begin());
  }
  template<class Entries>
  Bounds before(const Entries& entries, std::size_t next) const noexcept {
    return {next == 0 ? 0 : static_cast<std::uint64_t>(entries[next - 1].prefix),
            next == entries.size() ? count_
              : static_cast<std::uint64_t>(entries[next].prefix) + entries[next].delta - 1};
  }
  void compress_near(std::size_t inserted) noexcept {
    const auto limit = threshold();
    std::visit([&](auto& entries) {
      using Stored = typename std::decay_t<decltype(entries)>::value_type;
      if (entries.size() < 3) return;
      const auto end = std::min(inserted + 1, entries.size() - 2);
      std::size_t write = end + 1, read = end, begin;
      for (;;) {
        const auto tag = entries[read];
        const bool same_previous = read > 0 && entries[read - 1].value == tag.value;
        const bool same_successor = entries[write].value == tag.value;
        const std::uint64_t previous = read ? entries[read - 1].prefix : 0;
        const auto gap = static_cast<std::uint64_t>(entries[write].prefix) - previous;
        const bool remove = read && same_previous == same_successor &&
            entries[write].delta <= limit && gap <= limit - entries[write].delta;
        if (remove) {
          if (entries[read - 1].value < entries[write].value)
            uncertainty_ = std::max(uncertainty_, gap + entries[write].delta - 1);
        } else {
          entries[--write] = tag;
        }
        // Keeping an old tag left of the insertion restores the old successor
        // of every remaining pair. They were already ineligible for deletion.
        if ((!remove && read < inserted) || read == 0) {
          begin = read;
          break;
        }
        --read;
      }
      const auto removed = write - begin;
      if (removed) {
        std::memmove(entries.data() + begin, entries.data() + write,
                     (entries.size() - write) * sizeof(Stored));
        entries.resize(entries.size() - removed);
      }
    }, storage_);
  }
  void compress() noexcept {
    const std::uint64_t limit = threshold();
    std::visit([&](auto& entries) {
      // Keep extrema and duplicate-run endpoints. The surviving prefixes
      // stay fixed; removed gaps are implicit in their new differences.
      std::size_t write = entries.size();
      for (std::size_t read = entries.size(); read > 0;) {
        const auto tag = entries[--read];
        const bool same_previous = read > 0 && entries[read - 1].value == tag.value;
        const bool same_successor = write < entries.size() && entries[write].value == tag.value;
        const bool tie_boundary = same_previous != same_successor;
        const std::uint64_t previous = read ? entries[read - 1].prefix : 0;
        if (read != 0 && !tie_boundary && write < entries.size() &&
            entries[write].delta <= limit &&
            entries[write].prefix - previous <= limit - entries[write].delta) {
          // Absorbed into the successor without changing its prefix or delta.
        } else {
          entries[--write] = tag;
        }
      }
      const std::size_t size = entries.size() - write;
      uncertainty_ = 0;
      for (std::size_t i = 0; i < size; ++i) {
        entries[i] = entries[write + i];
        // No query falls between equal-valued tags. Calculate maximum width
        // during compaction rather than scanning a third time.
        if (i && entries[i - 1].value < entries[i].value)
          uncertainty_ = std::max(uncertainty_,
              static_cast<std::uint64_t>(entries[i].prefix) - entries[i - 1].prefix
                + entries[i].delta - 1);
      }
      entries.resize(size);
    }, storage_);
  }
  std::size_t capacity_;
  std::uint64_t count_ = 0, uncertainty_ = 0;
  Storage storage_;
};
}  // namespace detail

// Opt-in deterministic rank certificates; additional storage grows with the
// envelope. The compact SplineSketch remains available without this cost.
class CertifiedSplineSketch {
 public:
  struct RankEstimate {
    double estimate;
    std::uint64_t lower_rank, upper_rank;
    double max_error;  // Outward-rounded absolute error, in observations.
  };
  explicit CertifiedSplineSketch(std::size_t buckets = 128)
      : spline_(buckets), envelope_(buckets) {}

  CertifiedSplineSketch(const CertifiedSplineSketch&) = default;
  CertifiedSplineSketch(CertifiedSplineSketch&&) noexcept = default;
  CertifiedSplineSketch& operator=(CertifiedSplineSketch&&) noexcept = default;
  CertifiedSplineSketch& operator=(const CertifiedSplineSketch& other) {
    if (this != &other) {
      CertifiedSplineSketch updated(other);
      *this = std::move(updated);
    }
    return *this;
  }

  void add(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("value must be finite");
    if (count() == std::numeric_limits<std::uint64_t>::max())
      throw std::overflow_error("observation count overflow");
    if (value == 0) value = 0;  // Canonical positive zero.
    const auto insertion = envelope_.prepare_add(value);
    spline_.add(value);
    envelope_.commit_add(insertion);
  }
  std::uint64_t count() const noexcept { return spline_.count(); }
  std::size_t bucket_capacity() const noexcept { return spline_.bucket_capacity(); }
  std::size_t bucket_count() const noexcept { return spline_.bucket_count(); }
  std::size_t heavy_hitter_count() const noexcept { return spline_.heavy_hitter_count(); }
  std::size_t certificate_capacity() const noexcept { return envelope_.capacity(); }
  std::size_t certificate_size() const noexcept { return envelope_.size(); }
  std::uint64_t max_rank_uncertainty() const noexcept { return envelope_.uncertainty(); }

  double rank(double value) const {
    if (std::isnan(value)) throw std::invalid_argument("query must not be NaN");
    return estimate(value, envelope_.bounds(value));
  }
  RankEstimate rank_with_error(double value) const {
    if (std::isnan(value)) throw std::invalid_argument("query must not be NaN");
    const auto interval = envelope_.bounds(value);
    const double result = estimate(value, interval);
    if (interval.lower == interval.upper && interval.upper <= exact_integer_limit)
      return {result, interval.lower, interval.upper, 0};
    const double lower = outward_integer(interval.lower, false);
    const double upper = outward_integer(interval.upper, true);
    const double error = std::max(std::fabs(result - lower), std::fabs(result - upper));
    return {result, interval.lower, interval.upper,
            std::nextafter(error, std::numeric_limits<double>::infinity())};
  }
  // A uniform bound, including conversion rounding above 2^53 observations.
  double max_rank_error() const noexcept {
    double error = outward_integer(max_rank_uncertainty(), true);
    if (count() > exact_integer_limit) {
      const double n = static_cast<double>(count());
      error += std::nextafter(n, std::numeric_limits<double>::infinity()) - n;
      error = std::nextafter(error, std::numeric_limits<double>::infinity());
    }
    return error;
  }
  // Inverse of the certified estimated CDF; no value-error guarantee is made.
  double quantile(double q) const {
    if (!(q >= 0 && q <= 1)) throw std::invalid_argument("q must be in [0, 1]");
    if (!count()) throw std::logic_error("quantile of empty sketch");
    const double minimum = envelope_.entry(0).value;
    const double maximum = envelope_.entry(envelope_.size() - 1).value;
    if (q == 0) return minimum;
    if (q == 1) return maximum;
    const long double target = std::ceil(q * static_cast<long double>(count()));
    auto low = ordered_bits(minimum), high = ordered_bits(maximum);
    while (low < high) {
      const auto middle = low + (high - low) / 2;
      if (rank(from_ordered_bits(middle)) >= target) high = middle;
      else low = middle + 1;
    }
    return from_ordered_bits(low);
  }
  void consolidate() { spline_.consolidate(); }
  void resize(std::size_t capacity) {
    if (capacity < 6 || capacity > std::numeric_limits<std::size_t>::max() / 4)
      throw std::invalid_argument("bucket count must be at least 6");
    CertifiedSplineSketch updated(*this);
    updated.spline_.resize(capacity);
    updated.envelope_.resize(capacity);
    *this = std::move(updated);
  }
  void merge(const CertifiedSplineSketch& other) {
    if (std::numeric_limits<std::uint64_t>::max() - count() < other.count())
      throw std::overflow_error("observation count overflow");
    CertifiedSplineSketch updated(*this);
    updated.spline_.merge(other.spline_);
    updated.envelope_ = envelope_.merged(other.envelope_);
    *this = std::move(updated);
  }

 private:
  static constexpr std::uint64_t exact_integer_limit = std::uint64_t{1} << 53;
  static double outward_integer(std::uint64_t value, bool up) noexcept {
    const double rounded = static_cast<double>(value);
    if (value <= exact_integer_limit) return rounded;
    return std::nextafter(rounded, up ? std::numeric_limits<double>::infinity()
                                     : -std::numeric_limits<double>::infinity());
  }
  double estimate(double value, detail::RankEnvelope::Bounds interval) const {
    double result = spline_.rank(value);
    if (!std::isfinite(result))
      result = static_cast<double>(interval.lower + (interval.upper - interval.lower) / 2);
    return std::clamp(result, static_cast<double>(interval.lower),
                              static_cast<double>(interval.upper));
  }
  static std::uint64_t ordered_bits(double value) noexcept {
    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits & (std::uint64_t{1} << 63) ? ~bits : bits ^ (std::uint64_t{1} << 63);
  }
  static double from_ordered_bits(std::uint64_t ordered) noexcept {
    const auto bits = ordered & (std::uint64_t{1} << 63)
        ? ordered ^ (std::uint64_t{1} << 63) : ~ordered;
    double value;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
  }
  SplineSketch spline_;
  detail::RankEnvelope envelope_;
};
static_assert(std::is_nothrow_move_assignable<CertifiedSplineSketch>::value,
              "certified sketch commit must not throw");
}  // namespace splinesketch
