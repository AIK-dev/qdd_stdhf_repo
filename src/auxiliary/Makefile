### Makefile: Build QDD auxiliary programs.
# Simply invoke: 'make <program name WITHOUT the file extension>'
# E.g.: 'make spectr2' will build 'spectr2' out of 'spectr2.F90'
# and move it to '../../bin/'
COMPILER = ifort

%: %.f
	$(COMPILER) -o $@ $<
	mv -fv $@ ../../bin/

%: %.f90
	$(COMPILER) -o $@ $<
	mv -fv $@ ../../bin/

%: %.F90
	$(COMPILER) -o $@ $<
	mv -fv $@ ../../bin/
