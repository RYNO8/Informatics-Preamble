# Informatics-Preamble
Helpful datastructures and templates when doing informatics problems in C++. See `Preamble.cpp`  for example usage. There are also some helpful templates to help with standard formats of questions, such as Google Kickstart/Codejam.

My goal for this is to abstract the nitty gritty of implementations, so that informatics can be focused towards solving interesting problems instead. Hopefully this also makes informatics more accessible to newbies.

Please let me know if there are any errors.

Functionalities include
 - lazy persistent segtrees
 - directed / undirected graphs: graph traversal, MST, SSC, Toposort, bridges detection, ...
 - trees: traversals + LCA with jump pointers
 - mod arithmetic integers
 - matricies
 - 2D grids of values
 - convex hull trick


# todo list (no particular order)
[] implement coordinate geometry (eww)
[] vector .find .erase
[] wrapper for sqrt decomp
[] dynamic tree for Tree.h (look at references)
[] fast factorisation, fast isPrime, how?? Pollard rho was dodgy
[] change arguments to references if possible
[] write + fix tests
[] docstrings for ranges
[] segtree iterator
[] use `requires std::integral<T> || std::floating_point<T>` where possible

# Notes
 - It is probably not the fastest or most memory efficient, because I have prioritised verbosity over constant factor optimisations.
 - I hope that once its completed, its size would still be less than the maximum file size allowed
 - Only tested on gcc C++17

# C++pack.py
For such a large preamble, its helpful to `#include` it rather than copy pasting it at the
top of your code. `C++pack.py` helps to bubdle all code into a single file, which is
helpful when you want to submit (or reduce the submission file size). Some slightly hacky
lexical analysis stuff.
