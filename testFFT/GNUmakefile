

DIM=3

include ../Make.include


CFLAGS += -I$(FFT_HOME) -I$(RMDA_HOME) -I$(FFTW_HOME)/include -DDIM=$(DIM) 


LIB_FLAGS:= -L$(LIBS_LOCAL) -L$(FFTW_HOME)/lib -lfftw3 -lfft1D


SRCFILES:= $(wildcard $(RMDA_HOME)/*.cpp)
OBJS:=$(patsubst %.cpp,%.o, $(SRCFILES))

EOBJS:= FFTTest1D.o $(addsuffix $(DIM)D.o, FFTMDTest)



$(LIBS_LOCAL)/libfft1D.a:$(wildcard $(FFT_HOME)/*.cpp) $(wildcard $(FFT_HOME)/*.H) GNUmakefile
	cd $(FFT_HOME);make libfft1D DIM=1 

$(LIBS_LOCAL)/libfft3D.a:$(wildcard $(FFT_HOME)/*.cpp) $(wildcard $(FFT_HOME)/*.H) GNUmakefile
	cd $(FFT_HOME);make libfft3D DIM=3 


test1d: FFT1DTest.cpp GNUmakefile $(LIBS_LOCAL)/libfft1D.a
	  $(CXX) $(CFLAGS) FFT1DTest.cpp $(LIB_FLAGS) -o test1d.exe

testMD: GNUmakefile $(OBJS) FFTMDTest$(osuffix) $(LIBS_LOCAL)/libfft3D.a $(LIBS_LOCAL)/libfft1D.a 
	$(CXX) $(CFLAGS) $(OBJS) FFTMDTest$(osuffix) $(LIB_FLAGS) -lfft3D  -o test$(DIM)D.exe


clean:
	rm  -rf *.o *.exe *.d *.dSYM $(LIBS_LOCAL)/*.a $(FFT_HOME)/*.d $(RMDA_HOME)/*.o $(RMDA_HOME)/*.d
