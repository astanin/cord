###############################################################################
#
#                dieses Makefile ben"otigt GNU make!!!
#
###############################################################################
#
# wenn auf der SGI kompiliert wird: Wahl zwischen GNU g++ und SGI-CC
#USE_COMPILER = GNU
USE_COMPILER = SGI
#
# optimiert kompilieren oder zum Debuggen
#FLAGTYPE = OPT
FLAGTYPE = DEBUG
#
###############################################################################

ifeq ($(HOSTTYPE),iris4d)
ifeq ($(USE_COMPILER),GNU)
# mit GCC auf SGI
#  GXXFLAGS =	-nostdinc						      \
#	-I/net/xtreme/u1/local/sources/libg++-2.7.2/libstdc++                 \
#	-I/net/xtreme/u1/local/mips-sgi-irix5.3/include                       \
#	-I/net/xtreme/u1/local/lib/gcc-lib/mips-sgi-irix5.3/2.7.2.2/include   \
#	-I/usr/include							      \
#	-I/net/xtreme/u1/local/include
  OPTFLAGS = -O2 -finline-functions -funroll-loops
  DBGFLAGS = -O0
  CXX      = g++ $(GXXFLAGS)
  FLIB     = -lftn
  BLAS	   = -lblas
  DEPFLAG  = -MM
else
# mit SGI-Compiler auf SGI
#  CXX	   = CC -n32 -experimental	# O2
  CXX	   = CC -n32			# Octane/Origin
  OPTFLAGS = -O3 -OPT:alias=any -IPA:addressing=on -INLINE
  DBGFLAGS = -g
  FLIB     = -lftn
  BLAS	   = -lblas
  DEPFLAG  = -M
endif
endif


ifeq ($(HOSTTYPE),i486-linux)
  CXX      = g++ -W -ansi -pedantic -Wall
  OPTFLAGS = -O3 -funroll-loops
  DBGFLAGS = -g
  FLIB     = -lf2c
  BLAS	   = -lblas
  DEPFLAG  = -MM
endif



ifeq ($(FLAGTYPE),OPT)
  CXXFLAGS = $(OPTFLAGS)
endif

ifeq ($(FLAGTYPE),DEBUG)
  CXXFLAGS = $(DBGFLAGS)
endif

MKDEP = $(CXX) $(DEPFLAG)

###############################################################################



# hier werden die Dependencies gespeichert
DEPFILE	 = Dependfile

# Name, der beim tarren verwendet wird (ls = linear solver)
PCKNAME  = ls

# heutiges Datum
DATE     = $(shell date +%D | sed 's/\// /g' | awk '{print $$2$$1$$3}')

# Files, die getarred werden sollen
TARFILES = README *.h www example

# Name des Tarfiles
TARNAME	 = $(PCKNAME)$(DATE)

# alle Sourcen
SRC	 = antisym.cc laplace2d.cc laplace3d.cc argument.cc

# make all
all: depend antisym laplace2d laplace3d


laplace3d: laplace3d.o argument.o
	$(CXX) $(CXXFLAGS) -o laplace3d laplace3d.o argument.o \
	-lm $(BLAS) $(FLIB)

laplace2d: laplace2d.o argument.o
	$(CXX) $(CXXFLAGS) -o laplace2d laplace2d.o argument.o \
	-lm $(BLAS) $(FLIB)

antisym: antisym.o argument.o
	$(CXX) $(CXXFLAGS) -o antisym antisym.o argument.o \
	-lm $(BLAS) $(FLIB)

# tgz-File erzeugen (tarren und gzippen)
tgz:
	cd ..; tar -cf - $(TARFILES) | gzip >$(TARNAME).tgz;

# tgz-File erzeugen und symbolic link vom www-Verzeichnis darauf setzen
dist: tgz
	cd ../www; rm -f lsolver.tgz; ln -s ../$(TARNAME).tgz lsolver.tgz

# Dependencies erstellen
depend:
	$(MKDEP) $(SRC) > $(DEPFILE)

# aufr"aumen
clean:
	rm -f *.o *~ core

# mehr aufr"aumen
proper: clean
	rm -rf ii_files $(DEPFILE) antisym laplace2d laplace3d

# Dependencies einf"ugen
-include $(DEPFILE)
