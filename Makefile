VPATH = ./
MPATH=${HOME}/ponyPlayground/cmatrix

DEBUG=NO

ifeq (${DEBUG},YES)
	CLANGFLAGS=-O0 -g -DDEBUG
	PONYCFLAGS=--debug
else
	CLANGFLAGS=-DAVX -march=native -O3
	PONYCFLAGS=
endif


default: matrix.c matrix.pony test.pony
	@make libs
	@make matrix

libs: matrix.c
	clang -o matrix.o -c matrix.c ${CLANGFLAGS}
	ar rcs libmatrix.a matrix.o

matrix: matrix.pony test.pony
	ponyc ${PONYCFLAGS} --path=${MPATH}
