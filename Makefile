VPATH = ./ ./optional ./matrix
LIBPATH=${shell pwd}/bin

DEBUG=NO

ifeq (${DEBUG},YES)
	CLANGFLAGS=-DAVX2 -march=native -O0 -g -DDEBUG
	PONYCFLAGS=--debug
else
	CLANGFLAGS=-DAVX2 -march=native -O3
	PONYCFLAGS=
endif

PONYSRC=./matrix/matrix.pony

default: matrix/*
	@cp ./optional/main.pony ${PONYSRC} ./
	@make libs -s
	@make ponymath -s
	@rm main.pony ${notdir ${PONYSRC}}

test: matrix/*
	@cp ./optional/_test.pony ${PONYSRC} ./
	@make libs -s
	@make ponymath -s
	@rm ./_test.pony ${notdir ${PONYSRC}}

libs: matrix/matrix.c
	@clang -o ${LIBPATH}/matrix.o -c matrix/matrix.c ${CLANGFLAGS}
	@ar rcs ${LIBPATH}/libmatrix.a ${LIBPATH}/matrix.o

ponymath: ${PONYSRC} optional/*
	@ponyc ${PONYCFLAGS} --path=${LIBPATH}

.PHONY: clean
clean:
	rm -rf ./bin/* ./ponymath
