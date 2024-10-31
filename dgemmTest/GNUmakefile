#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

include ../Make.include

naive.exe: benchmark.cpp dgemm-naive.cpp GNUmakefile
	$(CXX) $(CFLAGS) $(BLASFLAGS) -o naive.exe benchmark.cpp dgemm-naive.cpp 
blas.exe: benchmark.cpp dgemm-blas.cpp GNUmakefile
	$(CXX) $(CFLAGS) $(BLASFLAGS) -o blas.exe benchmark.cpp dgemm-blas.cpp 

clean:
	rm -rf *.o *.exe *.d *.dSYM
