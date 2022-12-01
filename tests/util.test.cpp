void testPopcount(ll trials = 1e9) {
	for (ll i = 0; i < trials; ++i) {
		uint32_t x = uint_dis(rng);
		assert(popcount(x) == __builtin_popcount(x));
		uint64_t y = ull_dis(rng);
		assert(popcount(y) == __builtin_popcountll(y));
	}
}


int main() {
    testPopcount(1e4);
}