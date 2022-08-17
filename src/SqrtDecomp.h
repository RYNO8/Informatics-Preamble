#pragma once

namespace DS {
    // O(Q + N log N), where Q is the size of `queries_` and N is the total size of the interval to be considered
    // Mo's algorithm
    template<typename T> std::vector<T> chunkQueries(std::vector<std::pair<int, int>> &queries_, std::function<void(int)> add, std::function<void(int)> rem, std::function<T(int, int)> answer) {
        std::vector<std::pair<std::pair<int, int>, int>> queries;
        for (int i = 0; i < (int)queries_.size(); ++i) queries.push_back({ queries_[i], i });
        sort(queries.begin(), queries.end(), [](std::pair<std::pair<int, int>, int> a, std::pair<std::pair<int, int>, int> b) {
            return std::make_pair(a.first.first / SQRT_MAXN, a.first.second) < std::make_pair(b.first.first / SQRT_MAXN, b.first.second);
        });

        std::vector<T> ans(queries_.size());
        int l = 1, r = 1;
        add(1);
        for (auto &query : queries) {
            while (l > query.first.first) {
                l--;
                add(l);
            }
            while (l < query.first.first) {
                rem(l);
                l++;
            }

            while (r > query.first.second) {
                rem(r);
                r--;
            }
            while (r < query.first.second) {
                ++r;
                add(r);
            }
            ans[query.second] = answer(l, r);
        }

        return ans;
    }
};
