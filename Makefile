CFLAGS = -g -Wall -O3 -fopenmp
LINKFLAGS = -fopenmp -lm
PROJECT = Trabalho2
GROUPNAME = grupo-10a
REMOTEADDRESS = andromeda.lasdpc.icmc.usp.br
REMOTEPORT = 2270
CC = mpicc
CCSEQ = gcc
BINDIR = bin
SRCDIR = .
LIBDIR = .
REMOTEDIR = remotefiles
REMOTEFILES := $(wildcard $(REMOTEDIR)/*)
LIBS := $(wildcard $(LIBDIR)/*.h)
SRCS := $(wildcard $(SRCDIR)/*par.c)
OBJS := $(patsubst $(SRCDIR)/%.c, $(BINDIR)/%.o, $(SRCS))
SEQSRCS := $(wildcard $(SRCDIR)/*seq.c)
SEQOBJS := $(patsubst $(SRCDIR)/%.c, $(BINDIR)/%.o, $(SEQSRCS))

all : build

run : build 
	./$(PROJECT) 10 10 10 10

runseq : buildseq
	./$(PROJECT)_seq 10 10 10 10

mpi_run: build 
	mpirun -np 2 $(PROJECT) 10 10 10 10

runcompare: build buildseq
	mpirun -np 2 $(PROJECT) 10 10 10 10 > par.out
	./$(PROJECT)_seq 10 10 10 10 >seq.out
	diff seq.out par.out

$(BINDIR)/%par.o : $(SRCDIR)/%par.c $(LIBS)
	$(CC) -c $< -I $(LIBDIR) $(CFLAGS) -o $@

build : $(BINDIR) $(OBJS)
	$(CC) $(LINKFLAGS) $(OBJS) -o $(PROJECT)

$(BINDIR)/%seq.o : $(SRCDIR)/%seq.c $(LIBS)
	$(CCSEQ) -c $< -I $(LIBDIR) $(CFLAGS) -o $@

buildseq : $(BINDIR) $(SEQOBJS)
	$(CCSEQ) $(LINKFLAGS) $(SEQOBJS) -o $(PROJECT)_seq

remote :
	ssh $(GROUPNAME)@$(REMOTEADDRESS) -p $(REMOTEPORT)

sendfiles :
	scp -P $(REMOTEPORT) $(SRCS) $(LIBS) $(REMOTEFILES) $(GROUPNAME)@$(REMOTEADDRESS):/home/$(GROUPNAME)

$(BINDIR) :
	mkdir -p $(BINDIR)

clean :
	rm -rf $(BINDIR)
	rm -f $(PROJECT).zip
	rm -f $(PROJECT)
	rm -f $(PROJECT)_seq
	clear

.zip : clean
	zip $(PROJECT).zip $(SRCS) $(LIBS) $(SEQSRCS) $(REMOTEFILES) *.pdf

