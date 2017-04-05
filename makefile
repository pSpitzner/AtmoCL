CC=g++
platform=$(shell uname)
ifeq ($(platform),Darwin)
 CFLAGS = -std=c++11  -O3 -g -pedantic -I /opt/local/include -I ./classes/
 LFLAGS = -Wl,-framework,opencl -L/opt/local/lib -lgd
else ifeq ($(platform),Linux)
 CFLAGS = -std=c++11 -O3 -g -pedantic -I /net/nfs/opt/cuda-8.0/include -I /opt/local/include -I /opt/AMDAPPSDK-3.0/include -I ./classes/
 LFLAGS = -L /opt/AMDAPPSDK-3.0/lib/x86_64 -L /net/nfs/opt/cuda-8.0/lib64 -lOpenCL -L/opt/local/lib -lgd
endif

.SUFFIXES: .o .cpp
.cpp.o:  ; $(CC) -c $(CFLAGS) $*.cpp -o $*.o

OBJ= \
 ./classes/cllogger.o \
 ./classes/clcontext.o \
 ./classes/clbuffer.o \
 ./classes/clkernel.o \
 ./classes/clexport.o \
 ./classes/asystem.o \


all             : main.cpp $(OBJ)
	$(CC) $(CFLAGS) -o executable main.cpp ./classes/*.o $(LFLAGS)


clean           :
	(rm executable; rm ./classes/*.o; rm -r *.dSYM; rm ./output/img/*, rm ./output/logs/*)
