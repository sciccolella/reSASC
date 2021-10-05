test: test.o tap/basic.o
	gcc test.o tap/basic.o ../mt19937ar.o ../sastep.o ../tree.o ../vector.o ../utils.o -o test -lm -fopenmp

test.o: test.c
	gcc -c test.c -o test.o

tap/basic.o: tap/basic.c
	gcc -c tap/basic.c -o tap/basic.o

.PHONY: clean
clean:
	rm test tap/*.o *.o
