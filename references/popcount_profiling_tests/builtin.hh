static unsigned bitcount (TestType n)
{
    if(sizeof(n)==sizeof(int)) return __builtin_popcount(n);
    if(sizeof(n)==sizeof(long)) return __builtin_popcountl(n);
    return __builtin_popcountll(n);
}
