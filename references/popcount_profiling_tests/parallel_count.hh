static unsigned bitcount (TestType n)
{
    for(unsigned m = 0; (1u << m) < TEST_BITS; ++m)
    {
        TestType divisor = (TestType(1) << (1u << m)) + 1; // 3,5,17,257,...
        TestType mask = (~(TestType)0) / divisor;

        n = (n & mask) + ((n >> (1u << m)) & mask);
    }
    return n;
}
