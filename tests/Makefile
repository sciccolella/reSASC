UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CC = /usr/bin/gcc
	CFLAGS = -std=gnu99 -DNDEBUG -O3 -fopenmp
endif
ifeq ($(UNAME), Darwin)
	CC = /usr/local/bin/gcc-8
	CFLAGS = -std=c99 -DNDEBUG -O3 -fopenmp
endif

.PHONY: all

all: test

test: test.o tap/basic.o ../mt19937ar.o ../sastep.o ../tree.o ../utils.o ../vector.o
	@echo "* Linking test"
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c
	@echo '* Compiling $<'
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY: clean
clean:
	rm test tap/*.o *.o ../*.o
