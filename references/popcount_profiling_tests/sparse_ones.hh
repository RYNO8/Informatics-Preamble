static unsigned bitcount (TestType n)
{
	unsigned count = 0;
	while (n)
	{
		++count;
		n &= (n - 1); // Zero the lowest-order one-bit
	}
	return count;
}