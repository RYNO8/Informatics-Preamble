#pragma once

namespace DS {
    // @returns A mapping between original coordinates and compressed coordinates
    template<typename Iterator>
    std::map<typename std::iterator_traits<Iterator>::value_type, typename std::iterator_traits<Iterator>::value_type> coordCompressMap(Iterator begin, Iterator end, int startI = 0) {
        using value_type = typename std::iterator_traits<Iterator>::value_type;
        std::map<int, int> compressed;

        int i = startI;
        for (auto val = begin; val != end; ++val, ++i) {
            compressed[*val] = i;
        }
        return compressed;
    }

    // @returns A mapping between original coordinates and compressed coordinates
    template<typename T> std::map<T, T> coordCompressMap(T* begin, T* end, int startI = 0) {
        std::map<int, int> compressed;

        int i = startI;
        for (auto val = begin; val != end; ++val, ++i) {
            compressed[*val] = i;
        }
        return compressed;
    }

    // @returns A mapping between original coordinates and compressed coordinates
    template<typename T> std::map<T, T> coordCompressMap(std::set<T> toCompress, int startI = 0) {
        return coordCompressMap(toCompress.begin(), toCompress.end(), startI);
    }

    // @returns A mapping between original coordinates and compressed coordinates
    template<typename T> std::map<T, ll> coordCompressMap(std::vector<T> points, int startI = 0) {
        return coordCompressMap(points.begin(), points.end(), startI);
    }


    // Compresses coordinates in place
    template<typename Iterator> void coordCompressTransform(Iterator begin, Iterator end, int startI = 0){
        using value_type = typename std::iterator_traits<Iterator>::value_type;
        std::map<value_type, value_type> m = coordCompressMap(begin, end, startI);
        
        for (Iterator i = begin; i != end; i++) *i = m[*i];
    }

    // Compresses coordinates in place
    template<typename T> void coordCompressTransform(T* begin, T* end, int startI = 0) {
        std::map<T, T> m = coordCompressMap(begin, end, startI);

        for (T *i = begin; i != end; i++) *i = m[*i];
    }

    // Compresses coordinates in place
    template<typename T> void coordCompressTransform(std::vector<T> points, int startI = 0) {
        return coordCompressTransform(points.begin(), points.end(), startI);
    }
};