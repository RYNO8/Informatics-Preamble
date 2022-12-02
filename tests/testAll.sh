#!/bin/bash

# for FILE in *.test.cpp;
# do g++ $FILE -o temp && ./temp;
# done

echo "Util";
g++ util.test.cpp -O3 -Wall -o temp && ./temp; 
echo "Misc";
g++ misc.test.cpp -O3 -Wall -o temp && ./temp; 
echo "CHT";
g++ cht.test.cpp -O3 -Wall -o temp && diff <(echo "4
-1 10 -20
2 2 3 4" | ./temp) <(echo 9);
echo "Sqrt Decomp";
g++ sqrtdecomp.test.cpp -O3 -Wall -o temp && diff <(echo "22 6
87 79 87 95 83 85 67 72 95 68 79 71 69 95 71 79 68 95 72 67 85 83
4 22
9 9
10 17
10 11
9 17
10 18" | ./temp) <(echo "YES
YES
NO
NO
YES
YES");
echo "BIT";
g++ bit.test.cpp -O3 -Wall -o temp && ./temp;
echo "Ranges";
g++ ranges.test.cpp -O3 -Wall -o temp && ./temp;
echo "Segtree";
g++ segtree.test.cpp -O3 -Wall -o temp && ./temp;
echo "ModInt";
g++ modint.test.cpp -O3 -Wall -o temp && ./temp;