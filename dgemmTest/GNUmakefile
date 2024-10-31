#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

include ../Make.include

naive.exe: benchmark.cpp dgemm-naive.cpp GNUmakefile
	$(CXX) $(CFLAGS) -o naive.exe benchmark.cpp dgemm-naive.cpp $(BLASFLAGS) 
blas.exe: benchmark.cpp dgemm-blas.cpp GNUmakefile
	$(CXX) $(CFLAGS) -o blas.exe benchmark.cpp dgemm-blas.cpp $(BLASFLAGS)

clean:
	rm -rf *.o *.exe *.d *.dSYM
