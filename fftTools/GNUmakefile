#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

include ../Make.include

CFLAGS += -I$(FFT_HOME) -I$(RMDA_HOME) -I$(FFTW_HOME)/include -DDIM=$(DIM)


LIB_FLAGS = -L$(LIBS_LOCAL) -L$(FFTW_HOME)/lib -lfftw3 -lfft1D


#SRCFILES:= $(wildcard $(RMDA_HOME)/*.cpp)
#OBJS:=$(patsubst %.cpp,%.o, $(SRCFILES))


# Archiving tool on MACs.
ifeq ($(system),Darwin)
  LIBTOOL = libtool -static -s -o
endif


1DFFTSRC = $(wildcard $(FFT_HOME)/FFT1D*.cpp)
1DFFTSRC += $(FFT_HOME)/PowerItoI.cpp 

1DFFTOBJS = $(1DFFTSRC:.cpp=.o)


libfft3D: GNUmakefile FFTMD.o
	$(LIBTOOL) libfft3D.a FFTMD.o 
	mkdir -p ../lib;mv libfft3D.a $(LIBS_LOCAL);rm *.o

libfft1D:GNUmakefile $(wildcard $(FFT_HOME)/*.cpp) $(wildcard $(FFT_HOME)/*.H)
	make $(1DFFTOBJS)
	$(LIBTOOL) libfft1D.a $(1DFFTOBJS)
	mkdir -p ../lib;mv libfft1D.a $(LIBS_LOCAL);rm *.o

clean:
	rm -rf *.o *.d ../lib/*.a

-include $(1DFFTOBJS:.o=.d)
-include $(MDFFTOBJS:.o=.d)
