&GLOBAL
kxbox=64, kybox=64, kzbox=72, dx=0.6,
nelect=12, nion=3, nspdw=6, numspin=2,
kstate=20, b2occ=0.3d0, gamocc=30.0d0, deocc=0.4d0,
init_ao=.F.,
ifsicp=2,ipsptyp=1,
&END

&STATIC
e0dmp=6d0,epswf=0.5d0,epsoro=1.0e-8,
ifhamdiag=10, isaves=000, tstat=.F.,
ismax=5000,
&END

&DYNAMIC
nabsorb=4,ispherabso=1,iangabso=0,powabso=0.0625D0,
irest=0,
dt1=0.0125,itmax=0,
jdip=10,
jinfo=10,jdiporb=0,
jenergy=10,jstinf=10,
jesc=10,jnorms=10,
&END

&PERIO    
 cc1(6)   =-0.4842988,
 cc2(6)   =-0.6363824,
 crloc(6) = 0.6,
 r0g(6)=0.6,
 r1g(6)=0.6,
 radiong(6)=2.6,
 h0_11g(6)= 1.215075,
&END



dt1=0.0125,itmax=165400,
