#pragma once

// O(Q + N log N), where Q is the size of `queries_` and N is the total size of the interval to be considered
// Mo's algorithm
template<typename T> vector<T> chunkQueries(vector<pair<int, int>> &queries_, function<void(int)> add, function<void(int)> rem, function<T(int, int)> answer) {
    vector<pair<pair<int, int>, int>> queries;
    for (int i = 0; i < (int)queries_.size(); ++i) queries.push_back({ queries_[i], i });
    sort(queries.begin(), queries.end(), [](pair<pair<int, int>, int> a, pair<pair<int, int>, int> b) {
        return make_pair(a.first.first / SQRT_MAXN, a.first.second) < make_pair(b.first.first / SQRT_MAXN, b.first.second);
    });

    vector<T> ans(queries_.size());
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
