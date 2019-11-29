CFLAGS = -g -fopenmp
LINKFLAGS = -fopenmp -lm
PROJECT = Trabalho2
CC = mpicc
CCSEQ = gcc
BINDIR = bin
SRCDIR = .
LIBDIR = .
LIBS := $(wildcard $(LIBDIR)/*h)
SRCS := $(wildcard $(SRCDIR)/*par.c)
OBJS := $(patsubst $(SRCDIR)/%.c, $(BINDIR)/%.o, $(SRCS))
SEQSRCS := $(wildcard $(SRCDIR)/*seq.c)
SEQOBJS := $(patsubst $(SRCDIR)/%.c, $(BINDIR)/%.o, $(SEQSRCS))

all : build 

run: build 
	mpirun --hostfile hosts -np 10 $(PROJECT) 10 10 10 10

runcompare: build buildseq
	mpirun --hostfile hosts -np 10 $(PROJECT) 10 10 10 10 > par.out
	./$(PROJECT)_seq 10 10 10 10 >seq.out
	diff seq.out par.out

runseq : buildseq
	./$(PROJECT)_seq 10 10 10 10

$(BINDIR)/%par.o : $(SRCDIR)/%par.c $(LIBS)
	$(CC) -c $< -I $(LIBDIR) $(CFLAGS) -o $@

build : $(BINDIR) $(OBJS)
	$(CC) $(LINKFLAGS) $(OBJS) -o $(PROJECT)

$(BINDIR)/%seq.o : $(SRCDIR)/%seq.c $(LIBS)
	$(CCSEQ) -c $< -I $(LIBDIR) $(CFLAGS) -o $@

buildseq : $(BINDIR) $(SEQOBJS)
	$(CCSEQ) $(LINKFLAGS) $(SEQOBJS) -o $(PROJECT)_seq

$(BINDIR) :
	mkdir -p $(BINDIR)

clean :
	rm -rf $(BINDIR)
	rm -f $(PROJECT)
	rm -f $(PROJECT)_seq
	clear

