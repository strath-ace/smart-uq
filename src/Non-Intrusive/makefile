FC=g++
CC=gcc
CPP=g++
RM=rm -f
CPPFLAGS=-g -Wall -I../include/Math/Eigen/ -I../include/Non-Intrusive

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $<

all: polinterp_prb.x

PATH=. #./../src/Non-Intrusive

OBJ = ${PATH}/utils.o \
	${PATH}/sparse_grid_index.o \
	${PATH}/sparse_grid_cc_dataset.o \
	${PATH}/chebyshev_polynomial.o \
	${PATH}/multivariate_polynomials.o

polinterp_prb.x: $(OBJ) main.o
	$(CPP) $(CPPFLAGS) $(OBJ) main.o -o polinterp_prb.x

clean:
	$(RM) *~ *.o *.txt *.x \#*\#;