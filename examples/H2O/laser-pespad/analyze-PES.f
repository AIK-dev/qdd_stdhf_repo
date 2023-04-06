      program analyzePES
!     calculates an integrated PES of all directions given in pMP
!     and anisotropies

!     the angular weights for given arbritary orientations are 
!     calculated by searching for each point on the 
!     sphere the nearest orientation point
!     the fine grid of sphere points in initialized "irregular", 
!     i.e. for theta = 0 , one has only 1 phi 
!          for theta = 90, one has the maximum number of phis

      implicit real*8(a-h,o-z)
      parameter(kpoints=92)                        ! number of meas. points
      parameter(komega=10000)                       ! max. nr. of frequencies in PES
      character*(*) name
      parameter(name='H2O')                         ! set up the name

      
      parameter(ntheta =  361)                      ! fine grid
      parameter(itest  =   0)
      parameter(SMALL  = -100000.0D0)
      parameter(pi=3.141592653589793D0)
      parameter(dtheta = pi/(ntheta-1))
      dimension f(0:ntheta-1),nphi(0:ntheta-1)
      dimension w(kpoints),dir(3,kpoints),dire(3,kpoints),pes(kpoints)
      dimension dist(kpoints)
      dimension ksave(kpoints),thetz(kpoints),itheq(kpoints),
     &          pessave(kpoints)
      dimension energ(komega),yield(komega),beta(komega),
     &          beta4(komega),beta6(komega)

      character*6 asci6
      character*8 asci8

      test = 0.0D0
      wges = 0.0D0

      do k=1,kpoints
         w(k)     = 0.0D0
         ksave(k) = 0
      end do

!     init very fine grid in theta and phi on unit sphere
!     theta: equidistant
!     phi  : for theta=0,  1 phi
!            for theta=90, N phi's
!     every segment should have around the same surface area
      f(0)           = 2.0D0*pi*(1.0D0-cos(dtheta/2.0D0))
      f(ntheta-1)    = f(0)
      nphi(0)        = 1
      nphi(ntheta-1) = 1
      test           = f(0)+f(ntheta-1)

      do i=1,ntheta-2
        theta   = i*dtheta
        nphi(i) = DNINT(2.0D0*pi/f(0)*
     &        dabs(dcos(theta-dtheta/2.0D0)-dcos(theta+dtheta/2.0D0)))
        
        f(i)    = 2.0D0*pi/nphi(i)*
     &        dabs(dcos(theta-dtheta/2.0D0)-dcos(theta+dtheta/2.0D0))
        test    = test + f(i)*nphi(i)
      end do

!     open and prepare files
      open(unit=150,status='unknown',file='pPES.'//name)
      open(unit=151,status='unknown',file='iPES.'//name)
      open(unit=152,status='unknown',file='parametric.'//name)
      open(unit=153,status='unknown',file='distavphi.'//name)
      open(unit=155,status='unknown',file='distavphi-map.'//name)
      open(unit=154,status='unknown',file='velomap.'//name)

      write(151,'(a)') '# nr. of meas.points'
      write(151,'(a,1i6,a62)') '# ',kpoints,'weightage'
      write(152,'(a)') '# energy [Rydberg]   theta   distribution'
      write(153,'(a)') 
     &    '# energy [Rydberg]   theta   phi-averaged distribution'

!     init sampling, every direction point should lie on unit sphere
      read(150,*)
      read(150,*)
      do l=1,kpoints
         read(150,'(a,i6,a,3f12.1)') asci8,idum,asci6,
     &                                 dir(1,l),dir(2,l),dir(3,l)
         rad2      = dir(1,l)**2.0D0+dir(2,l)**2.0D0+dir(3,l)**2.0D0
         rad       = dsqrt(rad2) 
         dire(1,l) = dir(1,l)/rad
         dire(2,l) = dir(2,l)/rad
         dire(3,l) = dir(3,l)/rad
         thetz(l)  = dacos(dire(3,l))
         itheq(l)  = 1
      end do
 
      do m=1,kpoints
      do n=m+1,kpoints
        if(thetz(n) == thetz(m)) then
          itheq(m) = itheq(m) + 1
          thetz(n) = -1.0D0
        end if
      end do
      end do

!     searching for nearest orientation and sum up for weights
      do i=0,ntheta-1
         theta = i*dtheta

      do j=0,nphi(i)-1
         phi   = 2.0D0*pi*j/nphi(i)
         angsave = SMALL

      do k=1,kpoints
         ang = dire(1,k)*dsin(theta)*dcos(phi) + 
     &         dire(2,k)*dsin(theta)*dsin(phi) +
     &         dire(3,k)*dcos(theta)

         if(ang.eq.angsave) then       !  points with same distance share the weight surface
            ishare = ishare + 1
            ksave(ishare) = k
         end if
         if(ang.gt.angsave) then
            ishare     = 1
            angsave   = ang
            ksave(1)  = k
         end if
      end do
      if(ishare == 1) then
       w(ksave(1)) = w(ksave(1)) + f(i)
      else
       do l=1,ishare
       w(ksave(l)) = w(ksave(l)) + f(i)/ishare
       end do
      end if
      end do
      end do

      do k=1,kpoints
         wges = wges + w(k)
      end do

      write(*,*) '# point  x   y   z    w(point)'
      do ind=1,kpoints
      rad = dsqrt(dir(1,ind)**2.0D0+dir(2,ind)**2.0D0+dir(3,ind)**2.0D0)
         write(*,'(i4,10f8.4)') ind,dir(1,ind),dir(2,ind),dir(3,ind),
     &                       w(ind),wges,4.0D0*pi,rad
       end do

!     reading PES and sum up according weight factors
      do jnd=1,kpoints
        WRITE(151,'(a,i6,a,3f12.1,f12.5)') '# Point ',jnd , ' at : ',  
     &                dir(1,jnd),dir(2,jnd),dir(3,jnd),w(jnd)
      end do
      write(151,'(a8,4(a14))') '# omega ','int. PES','beta','beta4',
     &                          'beta6'

      pesc00tot = 0.0D0
      pesc20tot = 0.0D0
      pesc40tot = 0.0D0
      pesc60tot = 0.0D0
      dist      = 0.0D0
      do ind=1,komega
         read(150,*,err=98,end=98) energ(ind),pes
         kommax = ind
         if(ind == 1) delomega = energ(1)
         pesc00 = 0.0D0
         pesc20 = 0.0D0
         pesc40 = 0.0D0
         pesc60 = 0.0D0
         pessave     = 0.0D0
         dist        = dist + pes*delomega

         do k=1,kpoints
! +++ test
           if(itest == 1) then
            rad = dir(1,k)**2.0D0+dir(2,k)**2.0D0+dir(3,k)**2.0D0
            rad = 1.0D0   !dsqrt(rad)!            pes(k) = pes(k)*rad*rad
            pes(k) = 1.0D0/rad/rad*(1.0D0/dsqrt(4.0D0*pi) 
     &            + 1.0D0/dsqrt(5.0D0)*y20(dir(1,k),dir(2,k),dir(3,k))
     &            + 1.0D0/3.0D0*y40(dir(1,k),dir(2,k),dir(3,k))
     &            + 1.0D0/dsqrt(13.0D0)*y60(dir(1,k),dir(2,k),dir(3,k)))! dire**2.0  <=> beta = 2.0
           end if
! +++ end test
           pesc00 = pesc00 + w(k)*pes(k)/dsqrt(4.0D0*pi)
           pesc20 = pesc20 + w(k)*pes(k)*y20(dir(1,k),dir(2,k),dir(3,k))
           pesc40 = pesc40 + w(k)*pes(k)*y40(dir(1,k),dir(2,k),dir(3,k))
           pesc60 = pesc60 + w(k)*pes(k)*y60(dir(1,k),dir(2,k),dir(3,k))

           thetk = dacos(dire(3,k))
           do kl=1,kpoints
            if(thetk == thetz(kl)) pessave(kl) = pessave(kl) + pes(k)! averaging over phi 
           end do

         end do   ! k

         do km=1,kpoints
            if(thetz(km) >= 0D0) write(153,'(1f8.5,1pg16.5,1pg16.5)') 
     &             energ(ind),thetz(km)/pi*180.0D0,pessave(km)/itheq(km)
!           also velocity map possible
            if(thetz(km) >= 0D0) write(155,'(1f8.5,1pg16.5,1pg16.5)') 
     &        dsqrt(energ(ind))*dsin(thetz(km)),
     &        dsqrt(energ(ind))*dcos(thetz(km)),pessave(km)/itheq(km)
         end do
         write(153,*)
         write(155,*) 

         yield(ind) = dsqrt(4.0D0*pi)*pesc00                         ! the yield is the integrated PES
         beta(ind)  = dsqrt(5.0D0)*pesc20/pesc00
         beta4(ind) = 3.0D0*pesc40/pesc00
         beta6(ind) = dsqrt(13.0D0)*pesc60/pesc00
         write(151,'(1f8.5,1pg16.5,1pg16.5,2(1pg16.5))') energ(ind),
     &            yield(ind),beta(ind),beta4(ind),beta6(ind)
         pesc00tot = pesc00tot + pesc00*delomega
         pesc20tot = pesc20tot + pesc20*delomega
         pesc40tot = pesc40tot + pesc40*delomega
         pesc60tot = pesc60tot + pesc60*delomega
      end do
98    continue
      write(*,*) 'total yield = ',dsqrt(4.0D0*pi)*pesc00tot
      write(*,*) 'total beta  = ',dsqrt(5.0D0)*pesc20tot/pesc00tot
      write(*,*) 'total beta4 = ',3.0D0*pesc40tot/pesc00tot
      write(*,*) 'total beta6 = ',dsqrt(13.0D0)*pesc60tot/pesc00tot
      write(*,*) 'delomega    = ',delomega
      write(*,*) 'kpoints     = ',kpoints
      write(*,*) 'ntheta      = ',ntheta
      close(150)
      close(151)

      do k=1,kommax,10   !resolution
      do ith=0,180,9    !resolution
       theta = ith*pi/180.0D0
!      expansions end at beta2
       write(152,*) energ(k),theta/pi*180.0D0,
     &  yield(k)/(4D0*pi)*
     &  (1D0+beta(k)*(3D0*dcos(theta)**2D0-1D0)/2D0)
!    &  +beta4(k)*
!    &  (35D0*dcos(theta)**4D0-30D0*dcos(theta)**2D0+3D0)/8D0)
c       write(152,*) energ(k)*dsin(theta),energ(k)*dcos(theta),
c     &  yield(k)/(4D0*pi)*
c     &  (1D0+beta(k)*0.5D0*(3D0*dcos(theta)**2D0-1D0))
!     velocity maps
       write(154,'(5(1pg13.5))') 
     &  sqrt(energ(k))*dsin(theta),sqrt(energ(k))*dcos(theta),
     &  yield(k)/(4D0*pi)*sqrt(energ(k))*
     &  (1D0+beta(k)*(3D0*dcos(theta)**2D0-1D0)/2D0)
      end do
      write(152,*)
      write(154,*)
      end do 

      end program


      real*8 function y20(xx,yy,zz)
      implicit real*8(a-h,o-z)
      parameter(pi=3.141592653589793D0)
!      parameter(pi=3.141592654)
      if(xx.eq.0.0D0.and.yy.eq.0.0D0.and.zz.eq.0.0D0) then
         y20=0.0D0
      else
          y20=dsqrt(5.0D0/(16.0D0*pi))*
     &        (3.0D0*zz*zz/(xx*xx+yy*yy+zz*zz)-1.0D0)
      end if
      end


      real*8 function y40(xx,yy,zz)
      implicit real*8(a-h,o-z)
      parameter(pi=3.141592653589793D0)
!      parameter(pi=3.141592654)
      if(xx.eq.0.0D0.and.yy.eq.0.0D0.and.zz.eq.0.0D0) then
         y40=0.0D0
      else
         y40=3.0D0/16.0D0*dsqrt(1.0D0/pi)
     &       *(35.0D0*zz**4.0D0/(xx*xx+yy*yy+zz*zz)**2.0D0
     &         -30.0D0*zz*zz/(xx*xx+yy*yy+zz*zz)
     &                          +3.0D0)
      end if
      end


      real*8 function y60(xx,yy,zz)
      implicit real*8(a-h,o-z)
      parameter(pi=3.141592653589793D0)
!      parameter(pi=3.141592654)
      if(xx.eq.0.0D0.and.yy.eq.0.0D0.and.zz.eq.0.0D0) then
         y60=0.0D0
      else
         y60=dsqrt(13.0D0/pi)/32.0D0*
     &           (231.0D0*zz**6.0D0/(xx*xx+yy*yy+zz*zz)**3.0D0
     &           -315.0D0*zz**4.0D0/(xx*xx+yy*yy+zz*zz)**2.0D0
     &           +105.0D0*zz**2.0D0/(xx*xx+yy*yy+zz*zz)
     &            -5.0D0)
      end if
      end

