CFLAGS = -g
simple.tbin: makeIC
	makeIC 100 100 -1 > simple.tbin

makeIC: makeIC.o
	cc -o makeIC makeIC.o -lm

direct_rad: direct_rad.o
	cc -o direct_rad direct_rad.o -lm

clean:
	rm -rf makeIC makeIC.o simple.chk? simple.000* simple.log simple.timings
