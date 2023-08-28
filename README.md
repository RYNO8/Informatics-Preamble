# Informatics-Preamble
Helpful datastructures and templates when doing informatics problems in C++. See `Preamble.cpp`  for example usage. There are also some helpful templates to help with standard formats of questions, such as Google Kickstart/Codejam.

My goal for this is to abstract the nitty gritty of implementations, so that informatics can be focused towards solving interesting problems instead. Hopefully this also makes informatics more accessible to newbies.

I guess my secondary goal is now to learn C++ by trying to use it properly.

Please let me know if there are any logic errors or footguns here.

Functionalities include
 - lazy persistent segtrees
 - directed / undirected graphs: graph traversal, MST, SSC, Toposort, bridges detection, ...
 - trees: traversals + LCA with jump pointers
 - mod arithmetic integers
 - matricies
 - 2D grids of values
 - convex hull trick


# todo list (no particular order)
- [ ] combine `Grid.h` and `Matrix.h`
- [ ] implement coordinate geometry (eww)
- [ ] dynamic tree for `Tree.h` (look at references)
- [ ] fast factorisation, fast isPrime, how?? Pollard rho was dodgy
- [ ] use RAII and smart pointers for segtree
- [ ] docstrings for ranges
- [ ] Graph with self edges and multiple edges? (This is so rare though, I feel like support isn't worth it)
- [ ] only works on C++17
- [ ] cartesian product
- [ ] remove const value returns & reference returns 
- [ ] use the `noexpect` qualifier
- [ ] safe integer https://www.forrestthewoods.com/blog/perfect_prevention_of_int_overflows/

# Notes
 - It is probably not the fastest or most memory efficient, because I have prioritised verbosity over constant factor optimisations.
 - I hope that once its completed, its size would still be less than the maximum file size allowed
 - Only tested on gcc C++17

# C++pack.py
For such a large preamble, its helpful to `#include` it rather than copy pasting it at the
top of your code. `c++pack.py` helps to bundle all code into a single file, which is
helpful when you want to submit (or reduce the submission file size). Some slightly hacky
lexical analysis stuff.
