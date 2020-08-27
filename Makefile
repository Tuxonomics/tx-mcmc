CC = clang
CXX = clang++
PREFIX = /usr/local


CFLAGS = -std=c++11 -Wconversion -pedantic -Wno-c99-extensions

# Detect OS for MKL
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    MKL_OS = linux
endif
ifeq ($(UNAME_S),Darwin)
    MKL_OS = mac
endif

MKL_BASE = /opt/intel/compilers_and_libraries/$(MKL_OS)

ifeq ($(UNAME_S),Linux)
    MKL_PATH1 = $(MKL_BASE)/mkl/lib/intel64
    MKL_PATH2 = $(MKL_BASE)/lib/intel64
endif
ifeq ($(UNAME_S),Darwin)
    MKL_PATH1 = $(MKL_BASE)/mkl/lib
    MKL_PATH2 = $(MKL_BASE)/lib
endif

LFLAGS =  $(MKL_PATH1)/libmkl_intel_ilp64.a
LFLAGS =  $(MKL_PATH1)/libmkl_intel_lp64.a
LFLAGS += $(MKL_PATH1)/libmkl_intel_thread.a
LFLAGS += $(MKL_PATH1)/libmkl_core.a

ifeq ($(UNAME_S),Linux)
    LFLAGS += $(MKL_PATH1)/*.a $(MKL_PATH1)_lin/*.a
endif

LFLAGS += $(MKL_PATH2)/libiomp5.a


LFLAGS += -pthread -lm -ldl -lstdc++
DISABLED_WARNINGS = -Wno-switch


debug:		 		local_CFLAGS = -g -O0 -DDEBUG $(CFLAGS)
debug-asan:	 		local_CFLAGS = -g -O0 -DDEBUG $(CFLAGS)
release:     		local_CFLAGS = -O3 $(CFLAGS) -DRELEASE -march=native
tests-debug: 		local_CFLAGS = -O0 -fno-exceptions $(CFLAGS)
tests-asan:	 		local_CFLAGS = -O0 -fno-exceptions $(CFLAGS)
tests-debug/blas: 	local_CFLAGS = -O0 -fno-exceptions $(CFLAGS)
tests-asan/blas:	local_CFLAGS = -O0 -fno-exceptions $(CFLAGS)


TARGET = example
TEST = mcmc
TEST_TARGET = $(TEST)_tests
TEST_MAIN = $(TEST)_tests.cpp
TEST_LOG = $(TEST)_tests.log

all: 		debug
debug:   	clean $(TARGET)
debug-asan: clean $(TARGET)-asan
release: 	clean $(TARGET)


$(TARGET):
	$(CXX) src/$(TARGET).cpp -o $(TARGET) $(local_CFLAGS) $(LFLAGS) $(DISABLED_WARNINGS)

$(TARGET)-asan:
	$(CXX) src/$(TARGET).cpp -o $(TARGET) $(local_CFLAGS) $(LFLAGS) $(DISABLED_WARNINGS) -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./$(TARGET)



utilities/genTests:
	$(CXX) -o utilities/genTests utilities/gen_tests_main.cpp -std=c++11 -g

tests-debug: clean utilities/genTests
	./utilities/genTests $(shell find src -maxdepth 1 -type d) > $(TEST_MAIN)
	$(CXX) -o $(TEST_TARGET) $(TEST_MAIN) -DTEST $(LFLAGS) $(DISABLED_WARNINGS) $(local_CFLAGS) -g
	@./$(TEST_TARGET) 2> $(TEST_LOG)

tests-debug/blas: clean utilities/genTests
	./utilities/genTests $(shell find src -maxdepth 1 -type d) > $(TEST_MAIN)
	$(CXX) -o $(TEST_TARGET) $(TEST_MAIN) -DTEST -DUSE_BLAS $(LFLAGS) $(DISABLED_WARNINGS) $(local_CFLAGS) -g
	@./$(TEST_TARGET) 2> $(TEST_LOG)

tests-asan: clean utilities/genTests
	./utilities/genTests $(shell find src -maxdepth 1 -type d) > $(TEST_MAIN)
	$(CXX) -o $(TEST_TARGET) $(TEST_MAIN) -DTEST $(LFLAGS) $(DISABLED_WARNINGS) $(local_CFLAGS) -g -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./$(TEST_TARGET)

tests-asan/blas: clean utilities/genTests
	./utilities/genTests $(shell find src -maxdepth 1 -type d) > $(TEST_MAIN)
	$(CXX) -o $(TEST_TARGET) $(TEST_MAIN) -DTEST -DUSE_BLAS $(LFLAGS) $(DISABLED_WARNINGS) $(local_CFLAGS) -g -fsanitize=address
	ASAN_OPTIONS=symbolize=1 ASAN_SYMBOLIZER_PATH="$(shell which llvm-symbolizer)" ./$(TEST_TARGET)

tests: tests-debug
	@rm -f $(TEST_TARGET) $(TEST_MAIN)

tests/blas: tests-debug/blas
	@rm -f $(TEST_TARGET) $(TEST_MAIN)


clean:
	rm -f $(TARGET) $(TEST_TARGET) $(TEST_LOG) $(TEST_MAIN)


.PHONY: all clean debug release tests utilities/genTests
