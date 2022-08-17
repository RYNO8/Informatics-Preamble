#pragma once

// @returns A mapping between original coordinates and compressed coordinates
template<typename Iterator>
map<typename iterator_traits<Iterator>::value_type, typename iterator_traits<Iterator>::value_type> coordCompressMap(Iterator begin, Iterator end, int startI = 0) {
    using value_type = typename iterator_traits<Iterator>::value_type;
    map<int, int> compressed;

    int i = startI;
    for (auto val = begin; val != end; ++val, ++i) {
        compressed[*val] = i;
    }
    return compressed;
}

// @returns A mapping between original coordinates and compressed coordinates
template<typename T> map<T, T> coordCompressMap(T* begin, T* end, int startI = 0) {
    map<int, int> compressed;

    int i = startI;
    for (auto val = begin; val != end; ++val, ++i) {
        compressed[*val] = i;
    }
    return compressed;
}

// @returns A mapping between original coordinates and compressed coordinates
template<typename T> map<T, T> coordCompressMap(set<T> toCompress, int startI = 0) {
    return coordCompressMap(toCompress.begin(), toCompress.end(), startI);
}

// @returns A mapping between original coordinates and compressed coordinates
template<typename T> map<T, ll> coordCompressMap(vector<T> points, int startI = 0) {
    return coordCompressMap(points.begin(), points.end(), startI);
}


// Compresses coordinates in place
template<typename Iterator> void coordCompressTransform(Iterator begin, Iterator end, int startI = 0){
    using value_type = typename iterator_traits<Iterator>::value_type;
    map<value_type, value_type> m = coordCompressMap(begin, end, startI);
    
    for (Iterator i = begin; i != end; i++) *i = m[*i];
}

// Compresses coordinates in place
template<typename T> void coordCompressTransform(T* begin, T* end, int startI = 0) {
    map<T, T> m = coordCompressMap(begin, end, startI);

    for (T *i = begin; i != end; i++) *i = m[*i];
}

// Compresses coordinates in place
template<typename T> void coordCompressTransform(vector<T> points, int startI = 0) {
    return coordCompressTransform(points.begin(), points.end(), startI);
}