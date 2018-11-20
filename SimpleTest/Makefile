CFLAGS = -g
simple.tbin: makeIC
	makeIC 100 100 -1 > simple.tbin

makeIC: makeIC.o
	cc -o makeIC makeIC.o -lm

clean:
	rm -f makeIC makeIC.o
