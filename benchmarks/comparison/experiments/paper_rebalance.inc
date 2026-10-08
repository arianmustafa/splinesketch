  void rebalance(const Snapshot& before) {
    // Section 4.1's iteration-based Consolidate: choose disjoint edits,
    // rebuild once per iteration, and use the original snapshot for counters.
    using Candidate = std::pair<long double, std::size_t>;
    std::vector<Candidate> joins, splits;
    std::vector<unsigned char> edit;
    std::vector<Node> next;
    joins.reserve(capacity_);
    splits.reserve(capacity_);
    edit.reserve(capacity_);
    next.reserve(capacity_ + 2);
    for (;;) {
      joins.clear();
      splits.clear();
      edit.assign(nodes_.size(), 0); // 1: split, 2: left join, 3: right join
      std::size_t split_count = 0, join_count = 0;
      for (std::size_t i = 1; i < nodes_.size(); ++i) {
        if (can_split(i)) {
          if (nodes_[i].mass > bound()) {
            edit[i] = 1;
            ++split_count;
          }
          if (nodes_[i].mass > 0.01L * bound())
            splits.emplace_back(heuristic(i), i);
        }
        if (i + 1 < nodes_.size() && !nodes_[i].protected_threshold &&
            nodes_[i].mass + nodes_[i + 1].mass <= 0.75L * bound())
          joins.emplace_back(joined_heuristic(i), i);
      }
      std::sort(joins.begin(), joins.end());
      std::sort(splits.begin(), splits.end(), [](const Candidate& a, const Candidate& b) {
        return a.first > b.first || (a.first == b.first && a.second < b.second);
      });
      const std::size_t room = capacity_ - nodes_.size();
      const std::size_t required = split_count > room ? split_count - room : 0;
      std::size_t j = 0;
      for (; j < joins.size() && join_count < required; ++j) {
        const auto i = joins[j].second;
        if (edit[i] || edit[i + 1]) continue;
        edit[i] = 2;
        edit[i + 1] = 3;
        ++join_count;
      }
      if (join_count < required) {
        bucket_bound_factor_ *= 2;
        continue;
      }
      // The prototype keeps k/4+2 join candidates in reserve, and accepts
      // heuristic splits above 0.01 of the current bucket bound.
      std::size_t s = 0;
      while (s < splits.size() && split_count - join_count < room) {
        const auto i = splits[s++].second;
        if (edit[i]) continue;
        edit[i] = 1;
        ++split_count;
      }
      const std::size_t reserve = capacity_ / 4 + 2;
      const std::size_t join_limit = joins.size() > reserve ? joins.size() - reserve : 0;
      while (s < splits.size() && j < join_limit) {
        const auto i = splits[s].second;
        const auto q = joins[j].second;
        if (splits[s].first <= 1.5L * joins[j].first) break;
        if (!edit[i] && !edit[q] && !edit[q + 1] && i != q && i != q + 1) {
          edit[i] = 1;
          edit[q] = 2;
          edit[q + 1] = 3;
          ++split_count;
          ++join_count;
        }
        ++s;
        ++j;
      }
      if (split_count == 0 && join_count == 0) break;
      next.clear();
      for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (edit[i] == 2) continue;
        if (edit[i] == 1) {
          next.back().protected_threshold = true;
          next.emplace_back(midpoint(nodes_[i - 1].x, nodes_[i].x), 0, 0, 0, true);
        }
        next.push_back(nodes_[i]);
        if (edit[i] == 1) next.back().protected_threshold = true;
      }
      nodes_.swap(next);
      reestimate(before);
    }
  }
