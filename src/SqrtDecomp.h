#pragma once

namespace DS {
    // O(Q + N log N), where Q is the size of `queries_` and N is the total size of the interval to be considered
    // Mo's algorithm
    bool sqrtCmp(std::pair<std::pair<int, int>, int> a, std::pair<std::pair<int, int>, int> b) {
        return std::make_pair(a.first.first / SQRT_MAXN, a.first.second) < std::make_pair(b.first.first / SQRT_MAXN, b.first.second);
    }

    template<typename IndexType, typename AnsType, typename It>
    std::vector<AnsType> chunkQueries(
        It begin,
        It end,
        std::function<void(IndexType)> add,
        std::function<void(IndexType)> rem,
        std::function<AnsType(IndexType, IndexType)> answer
    ) {
        std::vector<std::pair<std::pair<IndexType, IndexType>, int>> queries;
        int Q = 0;
        for (auto val = begin; val != end; ++val, ++Q) queries.push_back({ *val, Q });
        sort(queries.begin(), queries.end(), sqrtCmp);

        std::vector<AnsType> ans(Q);
        IndexType l = 1, r = 1;
        add(1);
        for (auto &query : queries) {
            while (r > query.first.second) rem(r--);
            while (r < query.first.second) add(++r);
            while (l > query.first.first) add(--l);
            while (l < query.first.first) rem(l++);

            ans[query.second] = answer(l, r);
        }

        return ans;
    }
};
