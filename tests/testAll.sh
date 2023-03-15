#!/bin/bash

# for FILE in *.test.cpp &&
# do g++ $FILE -o temp && ./temp &&
# done

echo "Misc" &&
g++ misc.test.cpp -O3 -std=c++17 -Wall -o temp && ./temp && 

echo "BIT" &&
g++ bit.test.cpp -O3 -std=c++17 -Wall -o temp && ./temp &&

echo "CHT" && g++ cht.test.cpp -O3 -std=c++17 -Wall -o temp && ./temp &&

# echo "Geometry" &&

echo "Graph" && g++ graph.test.cpp -O3 -std=c++17 -Wall -o temp && ./temp &&

# echo "Grid" &&

# echo "Matrix" &&

echo "ModInt" &&
g++ modint.test.cpp -O3 -std=c++17 -Wall -o temp && ./temp &&

# echo "Polynomial" &&

echo "Ranges" &&
g++ ranges.test.cpp -O3 -std=c++17 -Wall -o temp && ./temp &&

echo "Segtree" &&
g++ segtree.test.cpp -O3 -std=c++17 -Wall -o temp && ./temp &&

echo "Sqrt Decomp" &&
g++ sqrtdecomp.test.cpp -O3 -std=c++17 -Wall -o temp && ./temp &&

# echo "Tree" &&

echo "Util" &&
g++ util.test.cpp -O3 -std=c++17 -Wall -o temp && ./temp && 

echo "Finished! Yay!"

rm temp