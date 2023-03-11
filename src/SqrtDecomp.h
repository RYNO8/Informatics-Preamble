#ifndef SQRTDECOMP_H
#define SQRTDECOMP_H
#include "Constants.h"

namespace DS {
    // O(Q sqrt N + N sqrt N), where Q is the size of `queries_` and N is the total size of the interval to be considered
    // Mo's algorithm
    // note: answer should accept an inclusive inclusive range
    template<typename IndexType, typename AnsType, typename It>
    std::vector<AnsType> chunkQueries(
        It begin,
        It end,
        const std::function<void(IndexType)> &add,
        const std::function<void(IndexType)> &rem,
        const std::function<AnsType(IndexType, IndexType)> &answer
    ) {
        struct Query {
            IndexType l, r;
            size_t out_i;
            bool operator<(const Query &b) const {
                return std::make_pair(l / SQRT_MAXN, r) < std::make_pair(b.l / SQRT_MAXN, b.r);
            }
        };

        std::vector<Query> queries;
        size_t Q = 0;
        for (auto val = begin; val != end; ++val, ++Q) queries.push_back({ val->first, val->second, Q });

        std::vector<AnsType> ans(Q);
        if (queries.empty()) return ans;
        
        sort(queries.begin(), queries.end());

        IndexType l = queries[0].l, r = queries[0].l;
        add(queries[0].l);
        for (auto &query : queries) {
            while (r > query.r) rem(r--);
            while (r < query.r) add(++r);
            while (l > query.l) add(--l);
            while (l < query.l) rem(l++);

            ans[query.out_i] = answer(l, r);
        }

        return ans;
    }
};

#endif