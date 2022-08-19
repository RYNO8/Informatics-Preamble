static unsigned bitcount (TestType n)
{
    TestType m1 = (~(TestType)0) / 3u;
    TestType m2 = (~(TestType)0) / 5u;
    TestType m4 = (~(TestType)0) / 17u;
    TestType h01 = (~(TestType)0) / 255u;

    n -= (n >> 1) & m1;             //put count of each 2 bits into those 2 bits
    n = (n & m2) + ((n >> 2) & m2); //put count of each 4 bits into those 4 bits 
    n = (n + (n >> 4)) & m4;        //put count of each 8 bits into those 8 bits 

    return (n * h01) >> (TEST_BITS-8); // n%255
}
