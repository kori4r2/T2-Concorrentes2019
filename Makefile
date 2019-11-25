CFLAGS = -g -Wall -O3 -fopenmp
LINKFLAGS = -fopenmp -lm
PROJECT = Trabalho2
REMOTEADDRESS = *ATUALIZEM ISSO AQUI*
REMOTEHOME = *ATUALIZEM ISSO AQUI*
REMOTEPORT = *ATUALIZEM ISSO AQUI*
CC = mpicc
DEBUGDIR = tests
BINDIR = bin
SRCDIR = src
LIBDIR = lib
REMOTEDIR = remotefiles
SERVERFILES := $(wildcard $(REMOTEDIR)/*)
LIBS := $(wildcard $(LIBDIR)/*h)
SRCS := $(wildcard $(SRCDIR)/*c)
OBJS := $(patsubst $(SRCDIR)/%.c, $(BINDIR)/%.o, $(SRCS))

all : build

build : $(BINDIR) $(OBJS)
	$(CC) $(LINKFLAGS) $(OBJS) -o $(PROJECT)

$(DEBUGDIR) :
	mkdir -p $(DEBUGDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)

$(BINDIR)/%.o : $(SRCDIR)/%.c $(LIBS)
	$(CC) -c $< -I $(LIBDIR) $(CFLAGS) -o $@

remote :
	ssh $(REMOTEADDRESS) -p $(REMOTEPORT)

sendfiles :
	scp -P $(REMOTEPORT) $(SRCS) $(LIBS) $(SERVERFILES) $(REMOTEADDRESS):$(REMOTEHOME)

clean :
	rm -rf $(BINDIR)
	rm -rf $(DEBUGDIR)
	rm -f $(PROJECT).zip
	rm -f $(PROJECT)
	#rm -f debug*.txt
	#rm -f resultado.txt
	clear

run : build
	./$(PROJECT)

mpi_run: build
	mpirun -np 2 $(PROJECT)

.zip : clean
	zip $(PROJECT).zip $(SRCS) $(LIBS) Makefile *.pdf

#debug: $(DEBUGDIR) all
#	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes mpirun -np 2 $(PROJECT) > $(DEBUGDIR)/output.txt 2> $(DEBUGDIR)/error.txt
#	diff resultadoCorreto.txt $(DEBUGDIR)/output.txt > $(DEBUGDIR)/diff.txt
