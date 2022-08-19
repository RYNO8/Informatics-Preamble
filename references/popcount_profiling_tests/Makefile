CXX=g++
CXXFLAGS=-Wall -W -Ofast -std=c++11
CXXFLAGS.clang=-Wall -W -O3 -std=c++11
LDFLAGS=-static

FILES=\
	test.cc \
	precomp.inc \
	naive.hh \
	dense_ones.hh \
	sparse_ones.hh \
	parallel_count.hh \
	nifty_parallel.hh \
	precomp.hh \
	builtin.hh wp2.hh wp3.hh

all: \
	test_int32_on_32bit \
	test_int64_on_32bit \
	test_int32_on_64bit \
	test_int64_on_64bit \
	clang_int32_on_32bit \
	clang_int64_on_32bit \
	clang_int32_on_64bit \
	clang_int64_on_64bit \
	test_int16_on_32bit \
	test_int8_on_32bit \
	test_int16_on_64bit \
	test_int8_on_64bit \
	clang_int16_on_32bit \
	clang_int8_on_32bit \
	clang_int16_on_64bit \
	clang_int8_on_64bit

test_int32_on_32bit: $(FILES)
	$(CXX) -o $@ test.cc $(CXXFLAGS) $(LDFLAGS) -m32 -msse -msse2 -DTestType=uint32_t
test_int64_on_32bit: $(FILES)
	$(CXX) -o $@ test.cc $(CXXFLAGS) $(LDFLAGS) -m32 -msse -msse2 -DTestType=uint64_t
test_int32_on_64bit: $(FILES)
	$(CXX) -o $@ test.cc $(CXXFLAGS) $(LDFLAGS) -m64 -DTestType=uint32_t
test_int64_on_64bit: $(FILES)
	$(CXX) -o $@ test.cc $(CXXFLAGS) $(LDFLAGS) -m64 -DTestType=uint64_t

clang_int32_on_32bit: $(FILES)
	clang++ -o $@ test.cc $(CXXFLAGS.clang) $(LDFLAGS) -m32 -msse -msse2 -DTestType=uint32_t
clang_int64_on_32bit: $(FILES)
	clang++ -o $@ test.cc $(CXXFLAGS.clang) $(LDFLAGS) -m32 -msse -msse2 -DTestType=uint64_t
clang_int32_on_64bit: $(FILES)
	clang++ -o $@ test.cc $(CXXFLAGS.clang) $(LDFLAGS) -m64 -DTestType=uint32_t
clang_int64_on_64bit: $(FILES)
	clang++ -o $@ test.cc $(CXXFLAGS.clang) $(LDFLAGS) -m64 -DTestType=uint64_t

test_int8_on_32bit: $(FILES)
	$(CXX) -o $@ test.cc $(CXXFLAGS) $(LDFLAGS) -m32 -msse -msse2 -DTestType=uint8_t
test_int16_on_32bit: $(FILES)
	$(CXX) -o $@ test.cc $(CXXFLAGS) $(LDFLAGS) -m32 -msse -msse2 -DTestType=uint16_t
test_int8_on_64bit: $(FILES)
	$(CXX) -o $@ test.cc $(CXXFLAGS) $(LDFLAGS) -m64 -DTestType=uint8_t
test_int16_on_64bit: $(FILES)
	$(CXX) -o $@ test.cc $(CXXFLAGS) $(LDFLAGS) -m64 -DTestType=uint16_t

clang_int8_on_32bit: $(FILES)
	clang++ -o $@ test.cc $(CXXFLAGS.clang) $(LDFLAGS) -m32 -msse -msse2 -DTestType=uint8_t
clang_int16_on_32bit: $(FILES)
	clang++ -o $@ test.cc $(CXXFLAGS.clang) $(LDFLAGS) -m32 -msse -msse2 -DTestType=uint16_t
clang_int8_on_64bit: $(FILES)
	clang++ -o $@ test.cc $(CXXFLAGS.clang) $(LDFLAGS) -m64 -DTestType=uint8_t
clang_int16_on_64bit: $(FILES)
	clang++ -o $@ test.cc $(CXXFLAGS.clang) $(LDFLAGS) -m64 -DTestType=uint16_t


results:
	sudo ./test_int64_on_64bit |tail -n +10 > results-int64-on-amd64
	sudo ./test_int32_on_64bit |tail -n +10 > results-int32-on-amd64
	sudo ./test_int16_on_64bit |tail -n +10 > results-int16-on-amd64
	sudo ./test_int8_on_64bit  |tail -n +10 > results-int8-on-amd64
	sudo ./test_int64_on_32bit |tail -n +10 > results-int64-on-i386
	sudo ./test_int32_on_32bit |tail -n +10 > results-int32-on-i386
	sudo ./test_int16_on_32bit |tail -n +10 > results-int16-on-i386
	sudo ./test_int8_on_32bit  |tail -n +10 > results-int8-on-i386
	sudo ./clang_int64_on_64bit |tail -n +10 > cesults-int64-on-amd64
	sudo ./clang_int32_on_64bit |tail -n +10 > cesults-int32-on-amd64
	sudo ./clang_int16_on_64bit |tail -n +10 > cesults-int16-on-amd64
	sudo ./clang_int8_on_64bit  |tail -n +10 > cesults-int8-on-amd64
	sudo ./clang_int64_on_32bit |tail -n +10 > cesults-int64-on-i386
	sudo ./clang_int32_on_32bit |tail -n +10 > cesults-int32-on-i386
	sudo ./clang_int16_on_32bit |tail -n +10 > cesults-int16-on-i386
	sudo ./clang_int8_on_32bit  |tail -n +10 > cesults-int8-on-i386
