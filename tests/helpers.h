#include <algorithm>
#include <vector>

// requires obj1::operator[], obj2::operator[], obj1.size(), obj2.size(), obj1::item == obj2::item
template<typename T1, typename T2>
bool unordered_eq(const T1 &obj1, T2 obj2) {
    std::vector<bool> seen(obj2.size());
    for (size_t i1 = 0; i1 < obj1.size(); ++i1) {
        for (size_t i2 = 0; i2 < obj2.size(); ++i2) {
            if (!seen[i2] && obj1[i1] == obj2[i2]) {
                seen[i2] = true;
                break; // break from 1 loop
            }
        }
    }
    return *std::min_element(seen.begin(), seen.end()) == true;
}