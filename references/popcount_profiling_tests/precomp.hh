/* Compile-time bit-counter */

template<typename T>
unsigned constexpr BitsIn(T n)
{
    return ParallelCount_CompileTime(n);
}

#define a(n) n+0, n+1, n+1, n+2          // 2 bits
#define b(n) a(n+0),a(n+1),a(n+1),a(n+2) // 4 bits
#define c(n) b(n+0),b(n+1),b(n+1),b(n+2) // 6 bits
#define d(n) c(n+0),c(n+1),c(n+1),c(n+2) // 8 bits
#define e(n) d(n+0),d(n+1),d(n+1),d(n+2) //10 bits
#define f(n) e(n+0),e(n+1),e(n+1),e(n+2) //12 bits
#define g(n) f(n+0),f(n+1),f(n+1),f(n+2) //14 bits
#define h(n) g(n+0),g(n+1),g(n+1),g(n+2) //16 bits
#define i(n) h(n+0),h(n+1),h(n+1),h(n+2) //18 bits
#define j(n) i(n+0),i(n+1),i(n+1),i(n+2) //20 bits
#define k(n) j(n+0),j(n+1),j(n+1),j(n+2) //22 bits

static const unsigned char bits_in[0x400000] alignas(64) = {
    // k(0)
    #include "precomp.inc"
};
static char bits_in_dummy[ sizeof(bits_in) ];
#undef k
#undef j
#undef i
#undef h
#undef g
#undef f
#undef e
#undef d
#undef c
#undef b
#undef a

template<unsigned PreCompBits, unsigned a=0, bool over = a >= TEST_BITS>
struct BitCountHelper
{
    // Using template recursion instead of a for-loop
    // will ensure that the compiler will completely
    // and totally unroll the loop.
    static unsigned bitcount(TestType n)
    {
        unsigned index = (n >> a) & ((1 << PreCompBits) - 1);
        return bits_in[ index ]
             + BitCountHelper<PreCompBits, a + PreCompBits>::bitcount(n);
    }

    static unsigned bitcount_clflush(TestType n)
    {
        unsigned index = (n >> a) & ((1 << PreCompBits) - 1);
        _mm_clflush( &bits_in[ index ] );
        return BitCountHelper<PreCompBits, a + PreCompBits>::bitcount_clflush(n);
    }

    static unsigned bitcount_clflush_dummy(TestType n)
    {
        unsigned index = (n >> a) & ((1 << PreCompBits) - 1);
        _mm_clflush( &bits_in_dummy[ index ] );
        return BitCountHelper<PreCompBits, a + PreCompBits>::bitcount_clflush_dummy(n);
    }
};
template<unsigned PreCompBits, unsigned a>
struct BitCountHelper<PreCompBits,a,true>
{
    static unsigned bitcount(TestType) { return 0; }
    static unsigned bitcount_clflush(TestType n)
        { return BitCountHelper<PreCompBits>::bitcount(n); }
    static unsigned bitcount_clflush_dummy(TestType n)
        { return BitCountHelper<PreCompBits>::bitcount(n); }
};
