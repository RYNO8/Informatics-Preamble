/* COPYRIGHT (C) 2013 Joel Yliluoma - http://iki.fi/bisqwit/ */
/* LICENSE: MIT */
/*
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <stdio.h>
#include <random>
#include <initializer_list>
#include <type_traits>
#include <algorithm>
#include <string>
#include <chrono>

#include <sched.h>
#include <unistd.h>

#ifndef TestType
typedef uint32_t TestType;
#endif

#ifndef CHAR_BITS
constexpr unsigned CHAR_BITS = 8;
#endif
constexpr unsigned TEST_BITS = CHAR_BITS * sizeof(TestType);

template<typename T>
struct Test_noclflush: public T
{
    //static unsigned bitcount_clflush(TestType n) { return T::bitcount(n); }
    //static unsigned bitcount_clflush_dummy(TestType n) { return T::bitcount(n); }
    static const bool No_clflush = true;
};

struct Test_Naive_base {
#include "naive.hh"
};
using Test_Naive = Test_noclflush<Test_Naive_base>;

struct Test_DenseOnes_base {
#include "dense_ones.hh"
};
using Test_DenseOnes = Test_noclflush<Test_DenseOnes_base>;

struct Test_SparseOnes_base {
#include "sparse_ones.hh"
};
using Test_SparseOnes = Test_noclflush<Test_SparseOnes_base>;

struct Test_NiftyParallel_base {
#include "nifty_parallel.hh"
};
using Test_NiftyParallel = Test_noclflush<Test_NiftyParallel_base>;

struct Test_Parallel_base {
#include "parallel_count.hh"
};
using Test_Parallel = Test_noclflush<Test_Parallel_base>;

struct Test_Builtin_base {
#include "builtin.hh"
};
using Test_Builtin = Test_noclflush<Test_Builtin_base>;

struct Test_WP2_base {
#include "wp2.hh"
};
using Test_WP2 = Test_noclflush<Test_WP2_base>;

struct Test_WP3_base {
#include "wp3.hh"
};
using Test_WP3 = Test_noclflush<Test_WP3_base>;

#include "precomp.hh"

using Test_Precomp2  = BitCountHelper<2>;
using Test_Precomp4  = BitCountHelper<4>;
using Test_Precomp8  = BitCountHelper<8>;
using Test_Precomp12 = BitCountHelper<12>;
using Test_Precomp16 = BitCountHelper<16>;
using Test_Precomp22 = BitCountHelper<22>;

#define EnumerateTests(m) \
    m(Test_Naive) m(Test_DenseOnes) m(Test_SparseOnes) \
    m(Test_Builtin) \
    /*m(Test_Parallel) m(Test_NiftyParallel) m(Test_WP3)*/ m(Test_WP2) \
    m(Test_Precomp2) m(Test_Precomp4) m(Test_Precomp8) \
    m(Test_Precomp12) m(Test_Precomp16) m(Test_Precomp22)

static struct TestMethod
{
    unsigned(*const method)(TestType);
    const char* const name;

    double nocache_penalty;
} TestMethods[] =
{
    #define m(name) { name::bitcount, #name, 0. },
    EnumerateTests(m)
    #undef m
};

static std::vector<TestType> corpus_sparse, corpus_dense, corpus_random;

static volatile unsigned store_temp;

void GenerateCorpus(unsigned num_values)
{
    const auto count = Test_Builtin::bitcount;
    std::mt19937_64 random;

    for(auto* tests: { &corpus_sparse, &corpus_dense, &corpus_random })
        tests->reserve(num_values);

    unsigned bin_limit   = num_values / (TEST_BITS / 4);
    std::vector<TestType> bins[TEST_BITS];
    for(auto& b: bins) b.reserve(bin_limit);

    unsigned char bits[TEST_BITS];
    for(unsigned b=0; b<TEST_BITS; ++b) bits[b] = b;

    for(unsigned num_bits = 1; num_bits <= TEST_BITS; ++num_bits)
        for(unsigned n=0; n<bin_limit; ++n)
        {
            // Generate a random number that has exactly (num_bits) bits.
            std::random_shuffle(bits, bits+TEST_BITS);

            TestType value = 0;
            for(unsigned b=0; b<num_bits; ++b)
                value |= TestType(1) << bits[b];

            if( count(value) != num_bits ) abort();

            bins[num_bits-1].push_back(value);
        }

    // Populate the "random" corpus evenly from each bin
    for(unsigned n=0; n<num_values; ++n)
        corpus_random.push_back
            ( bins[n % TEST_BITS][n / TEST_BITS] );

    unsigned quarter = TEST_BITS / 4;
    unsigned threeq  = quarter * 3;
    unsigned half = num_values/2;

    // Populate the "sparse" one such that first 25% of bins
    // form 50% of content, and remaining 75% of bins form the other 50%.
    for(unsigned n=0; n<num_values; ++n)
        if(n < half)
            corpus_sparse.push_back(bins[n % quarter][n / quarter]);
        else
            corpus_sparse.push_back(bins[quarter + n % threeq][(n-half) / threeq]);

    // Populate the "dense" one such that first 75% of bins
    // form 50% of content, and remaining 25% of bins form the other 50%.
    for(unsigned n=0; n<num_values; ++n)
        if(n < half)
            corpus_dense.push_back(bins[n % threeq][n / threeq]);
        else
            corpus_dense.push_back(bins[threeq + n % quarter][(n-half) / quarter]);
}

void SanityCheck()
{
    for(const auto* tests: { &corpus_sparse, &corpus_dense, &corpus_random })
        for(auto i: *tests)
        {
            unsigned expect = ~0u;
            for(const auto& m: TestMethods)
            {
                unsigned result = m.method(i);
                if(expect == ~0u)
                    expect = result;
                else if(result != expect)
                    fprintf(stderr,
                        "For %08llX, %s gave %d -- should be %d\n",
                        (long long)i, m.name, (int) result, (int) expect );
            }
        }
}

struct cpu_clock
{
    static const unsigned MHz = 2394;
    //typedef std::chrono::duration<uint64_t, std::ratio<1,MHz*1000000>> duration;
    //typedef duration::rep rep;
    //typedef duration::period period;
    //typedef std::chrono::time_point<cpu_clock, duration> time_point;
    typedef uint64_t time_point;
    struct period { static const uint64_t num = 1;
                    static const uint64_t den = MHz * 1000000; };
    struct duration {
        uint64_t val;
        duration(uint64_t v) : val(v) { }
        uint64_t count() const { return val; }
    };

    static time_point now()
    {
        __asm__ __volatile__ ("cpuid" : : : "ebx","ecx");
        uint64_t t;
        //t = __rdtsc();
#ifdef __x86_64
        unsigned hi, lo;
        __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
        t = uint64_t(lo) | (uint64_t(hi) << 32);
#else
        __asm volatile ("rdtsc" : "=A"(t));
#endif
        //fprintf(stderr, "Cycles read: %llu\n", (unsigned long long)t);
        return t;//time_point(duration(t));
    }
};

template<typename T>
static T timediff(T begin, T end)
{
    if(end < begin) return end + T(1ull << 32) - begin;
    return end-begin;
}

static void BeginRT()
{
    int pid = getpid();
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(0, &mask);
    sched_setaffinity(pid, sizeof(mask), &mask);
    sched_param p;
    p.sched_priority = sched_get_priority_max(SCHED_FIFO);
    sched_setscheduler(pid, SCHED_FIFO, &p);
    sched_setparam(pid, &p);
}
static void EndRT()
{
    int pid = getpid();
    sched_param p;
    p.sched_priority = (sched_get_priority_max(SCHED_OTHER)
                      + sched_get_priority_min(SCHED_OTHER)) / 2;
    sched_setscheduler(pid, SCHED_OTHER, &p);
}

namespace
{
    //using timer_t = std::chrono::high_resolution_clock;
    //using timer_t = std::chrono::steady_clock;
    using timer_t = cpu_clock;

    timer_t timer;
    double calc_time(const timer_t::duration& dur, std::size_t scale)
    {
        typedef timer_t::period p;
        return dur.count() * (double)p::num / p::den / scale;
    }
}

template<typename T>
typename std::enable_if< NULL != T::bitcount_clflush,  void>::type
    MakeCacheProfile(TestMethod& method)
{
    printf("%s", method.name);
    fflush(stdout);

    // Profile the method as-is.
    ////
    double time_alone;
    {BeginRT();
     const unsigned loops = 20;
     auto begin = timer.now();
     for(unsigned a=0; a<loops; ++a)
      for(auto i: corpus_random)
       store_temp = T::bitcount(i);
     auto end   = timer.now();
     EndRT();
     time_alone = calc_time(timediff(begin,end), loops*corpus_random.size());
    }
    ////

    // Profile the method with no-op cache flushing.
    ////
    double time_dummyflush;
    {BeginRT();
     const unsigned loops = 8;
     auto begin = timer.now();
     for(unsigned a=0; a<loops; ++a)
      for(auto i: corpus_random)
       store_temp = T::bitcount_clflush_dummy(i);
     auto end   = timer.now();
     EndRT();
     time_dummyflush = calc_time(timediff(begin,end), loops*corpus_random.size());
    }
    ////

    // Profile the method with real cache flushing.
    ////
    double time_realflush;
    {BeginRT();
     const unsigned loops = 8;
     auto begin = timer.now();
     for(unsigned a=0; a<loops; ++a)
      for(auto i: corpus_random)
       store_temp = T::bitcount_clflush(i);
     auto end   = timer.now();
     EndRT();
     time_realflush = calc_time(timediff(begin,end), loops*corpus_random.size());
    }
    ////

    method.nocache_penalty = time_realflush - time_dummyflush;
    printf("-- Alone=%g, Dummy=%g, Flush=%g; Penalty=%g\n",
        time_alone, time_dummyflush, time_realflush,
        method.nocache_penalty);
}

template<typename T>
typename std::enable_if<T::No_clflush,  void>::type
    MakeCacheProfile(TestMethod& ) { }

static void MakeCacheProfiles()
{
    unsigned n=0;
    #define m(name) MakeCacheProfile<name>(TestMethods[n++]);
    EnumerateTests(m)
    #undef m
}

template<typename Functor>
void DoTimingTest(Functor& func, const char* name, double nocache_penalty)
{
    std::string n = name + 5;
    if(nocache_penalty)
        n += " @ Bad cache";
    else if(n.substr(0,7) == "Precomp")
        n += " @ Good cache";

    printf("%-22s", n.c_str());

    for(const auto* p: {&corpus_random, &corpus_dense, &corpus_sparse})
    {
        const unsigned loops = 50;
        BeginRT();
        auto begin = timer.now();
        for(unsigned n=0; n<loops; ++n)
        {
            for(auto i: *p) store_temp = func(i);
        }
        auto end   = timer.now();
        EndRT();
        double time = calc_time(timediff(begin,end), loops*corpus_random.size());
        time += nocache_penalty;
        // time = how many seconds per item
        printf("%10.2f Mcps", 1/time / 1e6);
    }

    printf("\n");
}

int main()
{
    printf("Generating random numbers...\n");
    GenerateCorpus(400000);
    for(unsigned a=0; a<sizeof(bits_in_dummy); ++a) bits_in_dummy[a] = ~a;

    printf("Sanity checking...\n");
    BeginRT();
    SanityCheck();
    EndRT();

    printf("Performing cache profile...\n");
    MakeCacheProfiles();

    printf("\n%d-bit integers, running on %d-bit platform\n",
        TEST_BITS,
        sizeof(long)==8 ? 64 : 32);

    printf("%-22s%15s%15s%15s\n", "---","---","---","---");
    printf("%-22s%15s%15s%15s\n",
        "Method name",
        "Random data",
        "Dense data",
        "Sparse data");
    printf("%-22s%15s%15s%15s\n", "---","---","---","---");

    #define m(name) DoTimingTest(name::bitcount, #name, 0.);
    EnumerateTests(m)
    #undef m

    unsigned n=0;
    #define m(name) { \
        if(TestMethods[n].nocache_penalty) \
            DoTimingTest(name::bitcount, #name, TestMethods[n].nocache_penalty); \
        ++n; }
    EnumerateTests(m)
    #undef m

    printf("---\n");
    printf("Mcps = million counts per second\n");
    printf("Random data has equal chance for any number of bits in it.\n");
    printf("Dense data has 75%% chance that more than half of its bits are set.\n");
    printf("Sparse data has 75%% chance that less than half of its bits are set.\n");
}
