
# Makefile to build example MPI programs
# #
#
 CC=mpicc
#
 COMP=GNU
  ifeq ($(COMP), GNU)
   CFLAGS=-std=c99 -Wall
   endif

   EXE1=stencil.exe
#   #EXE2=derived-type.exe
EXES=$(EXE1) $(EXE2)

all: $(EXES)

$(EXES): %.exe : %.c
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean

clean:
	#\rm -f $(EXES)
	#\rm -f *.o

