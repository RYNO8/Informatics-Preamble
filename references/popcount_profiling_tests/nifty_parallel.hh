static unsigned bitcount(TestType n)
{
    for(unsigned m = 0; m < 3; ++m)
    {
        TestType divisor = (1 << (1 << m)) + 1; // 3,5,17,257,...
        TestType mask = (~(TestType)0) / divisor;

        //printf("divisor=%d, n = (n & 0x%08X) + ((n >> %d) & 0x%08X)\n",
        //    divisor, mask,1<<m, mask);

        n = (n & mask) + ((n >> (1 << m)) & mask);
    }
    //printf("n=%lld; n%%15=%d, n%%31=%d, n%%63=%d, n%%255=%d, n%%65535=%d\n",
    //    (long long)n, n%15, n%31, n%63, n%255, n%65535);

    TestType h01 = (~(TestType)0) / 255u;

    return (n * h01) >> (TEST_BITS-8); // n%255
    //return n % 255;
}
