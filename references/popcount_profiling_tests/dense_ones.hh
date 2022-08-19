static unsigned bitcount (TestType n)
{
	unsigned count = TEST_BITS;
	n = ~n;
	while (n)
	{
		--count;
		n &= (n - 1); // Zero the lowest-order one-bit
	}
	return count;
}