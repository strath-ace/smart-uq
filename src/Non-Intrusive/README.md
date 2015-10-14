
README

==================================

## Usage

To obtain the coefficients using a least-square method:

```
coe = full_multivariate_polynomials ( dim_num, deg_max, pts, trues );
```

where

- dim_num : spatial dimention

- deg_max : degree of expansion

- pts : interpolation points (of size dim_num x point_num )

- trues : true values at interpolation spoints (of size point_num)


## Example of *makefile*.

CPPFLAGS=-g -Wall -I../include/eigen/ -I../suit/alcor

%.o: %.cpp
        $(CPP) $(CPPFLAGS) -c $<

all: main.x

main.x: utils.o sparse_grid_index.o sparse_grid_cc_dataset.o \
        chebyshev_polynomial.o multivariate_polynomials.o main.o
        $(CPP) $(CPPFLAGS) utils.o sparse_grid_index.o sparse_grid_cc_dataset.o \
        chebyshev_polynomial.o multivariate_polynomials.o main.o -o main.x

clean:
        $(RM) *~ *.o *.txt *.x \#*\#;
