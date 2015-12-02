C NCLFORTSTART                                                        
      subroutine wrfcttcalc(prs,tk,qci,qcw,qvp,ght,ter,ctt,
     &                      haveqci,nz,ns,ew)

      implicit none
      integer nz,ns,ew,haveqci 
      double precision    ght(ew,ns,nz)
      double precision    prs(ew,ns,nz),tk(ew,ns,nz)
      double precision    qci(ew,ns,nz),qcw(ew,ns,nz)
      double precision    qvp(ew,ns,nz)
      double precision    ctt(ew,ns),ter(ew,ns)
c      double precision    znfac(nz)
C NCLEND
c
c
      integer i,j,k,mjx,miy,mkzh,ripk,wrfout
      double precision    vt,rgas,grav,opdepthu,opdepthd,dp
      double precision    ratmix,eps,arg1,arg2,agl_hgt,ussalr
      double precision    abscoefi,abscoef,fac,prsctt,celkel
c      double precision    ght(ew,ns,nz),stuff(ew,ns)
      double precision    pf(ns,ew,nz),p1,p2
c
c
       mjx      =   ew 
       miy      =   ns 
       mkzh     =   nz 
       eps      = 0.622d0
       ussalr   = .0065d0      ! deg C per m
       rgas     = 287.04d0     !J/K/kg
       grav     = 9.81d0
       abscoefi = .272d0  ! cloud ice absorption coefficient in m^2/g
       abscoef  =.145d0   ! cloud water absorption coefficient in m^2/g
       celkel   = 273.15d0 
       wrfout = 1


cCalculate the surface pressure 
       do j=1,ew
       do i=1,ns
           ratmix     = .001d0*qvp(j,i,1)
           arg1       = eps + ratmix
           arg2       = eps*(1.+ratmix)
           vt         =  tk(j,i,1) * arg1/arg2 !Virtual temperature
           agl_hgt    = ght(j,i,nz) - ter(j,i)
           arg1       = -grav/(rgas*ussalr)
           pf(i,j,nz) = prs(j,i,1)*
     &                        (vt/(vt+ussalr*(agl_hgt)))**(arg1)
       enddo
       enddo


c
       do j=1,ew
       do i=1,ns
          do k=1,nz-1
            ripk = nz-k+1
            pf(i,j,k)=.5d0*(prs(j,i,ripk)+prs(j,i,ripk-1))
          enddo
      enddo
      enddo

      do 190 j=1,ew
      do 190 i=1,ns
         opdepthd=0.d0
         k=0

c
c      Integrate downward from model top, calculating path at full
c      model vertical levels.
c
   20    opdepthu=opdepthd
         k=k+1
         ripk = nz-k+1

         if (k.eq.1) then
            dp=200.d0*(pf(i,j,1)-prs(j,i,nz))  ! should be in Pa
         else
            dp=100.d0*(pf(i,j,k)-pf(i,j,k-1))  ! should be in Pa
         endif 
         if (haveqci .eq. 0) then
            if (tk(i,j,k).lt.celkel) then
c             Note: abscoefi is m**2/g, qcw is g/kg,
c                   so no convrsion needed
               opdepthd=opdepthu+abscoefi*qcw(j,i,k)*dp/grav
            else
               opdepthd=opdepthu+abscoef*qcw(j,i,k)*dp/grav
            endif
         else
            opdepthd=opdepthd+(abscoef*qcw(j,i,ripk)+
     &                         abscoefi*qci(j,i,ripk))*dp/grav
         endif
         
          if (opdepthd.lt.1..and.k.lt.nz) then
            goto 20
         elseif (opdepthd.lt.1..and.k.eq.nz) then
            prsctt=prs(j,i,1)
         else
            fac=(1.-opdepthu)/(opdepthd-opdepthu)
            prsctt=pf(i,j,k-1)+fac*(pf(i,j,k)-pf(i,j,k-1))
            prsctt=min(prs(j,i,1),max(prs(j,i,nz),prsctt))
         endif

         do 30 k=2,nz
            ripk = nz-k+1
            p1   = prs(j,i,ripk+1)
            p2   = prs(j,i,ripk)
            if (prsctt .ge. p1 .and. prsctt .le .p2) then
               fac=(prsctt-p1)/(p2-p1)
               arg1 = fac*(tk(j,i,ripk)-tk(j,i,ripk+1))-celkel
               ctt(j,i) = tk(j,i,ripk+1)+ arg1
               goto 40
            endif
   30    continue 
   40    continue
 190  continue
      return
      end
