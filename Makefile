VPATH = ./ ./test ./matrix
LIBPATH=${shell pwd}/bin

DEBUG=NO

ifeq (${DEBUG},YES)
	CLANGFLAGS=-DAVX -march=native -O0 -g -DDEBUG
	PONYCFLAGS=--debug
else
	CLANGFLAGS=-DAVX -march=native -O3
	PONYCFLAGS=
endif


default: matrix.c matrix.pony
	@cp ./optional/main.pony ./
	@make libs -s
	@make ponymath -s
	@rm ./main.pony

test: matrix.c matrix.pony
	@cp ./optional/_test.pony ./
	@make libs -s
	@make ponymath -s
	@rm ./_test.pony

libs: matrix.c
	@clang -o ${LIBPATH}/matrix.o -c matrix.c ${CLANGFLAGS}
	@ar rcs ${LIBPATH}/libmatrix.a ${LIBPATH}/matrix.o

ponymath: matrix.pony optional/main.pony optional/_test.pony
	@ponyc ${PONYCFLAGS} --path=${LIBPATH} --bin-name $@

.PHONY: clean
clean:
	rm -rf ./bin/* ./ponymath
