! Lennard-Jones for Ar
PARAMETER(v0_lj=4*7.61E-4)/*4*120 kelvin IN units of[ry]*/
PARAMETER(sigma=6.54)/*0.34 nm IN units of[bohr]*/

DO i = 20, 140, 2
  x = i*0.1
  pot = v0_lj*((sigma/x)**12.0 - ((sigma/x)**6.0))
  WRITE (6, '(1x,2f6.1,1pg13.5)') x, x/2.0, pot
END DO
STOP
END
