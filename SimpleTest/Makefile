CFLAGS = -g
simple.tbin: makeIC
	makeIC 100 100 -1 > simple.tbin

makeIC: makeIC.o
	cc -o makeIC makeIC.o -lm

direct_rad: direct_rad.o
	cc -o direct_rad direct_rad.o -lm

direct_thick: direct_thick.o
	cc -o direct_thick direct_thick.o -lm

clean:
	rm -rf makeIC makeIC.o simple.chk? simple.000* simple.log simple.timings
	rm -rf direct_rad.o direct_rad
