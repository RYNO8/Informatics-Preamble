static unsigned bitcount (TestType n)
{
	unsigned count = 0;
	while (n) {
		count += (n & 1);
		n >>= 1;
	}
	return count;
}
