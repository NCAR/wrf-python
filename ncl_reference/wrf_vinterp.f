CThe subroutines in this file were taken directly from RIP code written
C by Dr. Mark Stoelinga.  They were modified by Sherrie
C Fredrick(NCAR/MMM) to work with NCL February 2015.
C NCLFORTSTART
      subroutine wrf_monotonic(out,in,lvprs,cor,idir,delta,
     &                         ew,ns,nz,icorsw)
       implicit none
       integer idir,ew,ns,nz,icorsw
       double precision delta
       double precision in(ew,ns,nz),out(ew,ns,nz)
       double precision lvprs(ew,ns,nz),cor(ew,ns)
C NCLEND

       integer i,j,k,ripk,k300

       do j=1,ns
       do i=1,ew
          if (icorsw.eq.1.and.cor(i,j).lt.0.) then
             do k=1,nz
                in(i,j,k)=-in(i,j,k)
              enddo
          endif   


c
c   First find k index that is at or below (height-wise) the 300 hPa
c      level.
c
      do k = 1,nz
         ripk =  nz-k+1
         if (lvprs(i,j,k) .le. 300.d0) then
            k300=k
            goto 40
         endif
      enddo
c
 40   continue

       do k = k300, 1,-1
         if (idir.eq.1) then
            out(i,j,k)=min(in(i,j,k),in(i,j,k+1)+delta)
         elseif (idir.eq.-1) then
            out(i,j,k)=max(in(i,j,k),in(i,j,k+1)-delta)
         endif
       enddo


       do k = k300+1, nz
         if (idir.eq.1) then
            out(i,j,k)=max(in(i,j,k),in(i,j,k-1)-delta)
         elseif (idir.eq.-1) then
            out(i,j,k)=min(in(i,j,k),in(i,j,k-1)+delta)
         endif
       enddo
      
      end do
      end do
      
      return
      end 

c--------------------------------------------------------------------

C NCLFORTSTART
      FUNCTION wrf_intrp_value (wvalp0,wvalp1,vlev,vcp0,vcp1,icase)
      implicit none

      integer icase
      double precision wvalp0,wvalp1,vlev,vcp0,vcp1
C NCLEND
      double precision valp0,valp1,rvalue,rgas,ussalr,sclht

      double precision wrf_intrp_value,chkdiff

      rgas    = 287.04d0     !J/K/kg
      ussalr  = 0.0065d0      ! deg C per m
      sclht   = rgas*256.d0/9.81d0

      valp0 = wvalp0
      valp1 = wvalp1
      if ( icase .eq. 2) then  !GHT
           valp0=exp(-wvalp0/sclht)
           valp1=exp(-wvalp1/sclht)
      end if

      chkdiff = vcp1 - vcp0
      if(chkdiff .eq. 0) then
         print *,"bad difference in vcp's"
         stop
      end if
 
      rvalue = (vlev-vcp0)*(valp1-valp0)/(vcp1-vcp0)+valp0
      if (icase .eq. 2) then  !GHT
          wrf_intrp_value = -sclht*log(rvalue)
      else
          wrf_intrp_value = rvalue
      endif

      return
      end
c------------------------------------------------------------
C NOTES:
c      vcarray is the array holding the values for the vertical 
c      coordinate.
c              It will always come in with the dimensions of 
c              the staggered U and V grid.  
C NCLFORTSTART

      subroutine  wrf_vintrp(datain,dataout,pres,tk,qvp,ght,terrain,
     &                       sfp,smsfp,vcarray,interp_levels,numlevels,
     &                       icase,ew,ns,nz,extrap,vcor,logp,rmsg)


      implicit none
      integer   ew,ns,nz,icase,extrap
      integer   vcor,numlevels,logp
      double precision   datain(ew,ns,nz),pres(ew,ns,nz),tk(ew,ns,nz)
      double precision   ght(ew,ns,nz)
      double precision   terrain(ew,ns),sfp(ew,ns),smsfp(ew,ns)
      double precision   dataout(ew,ns,numlevels),qvp(ew,ns,nz)
      double precision   vcarray(ew,ns,nz)
      double precision   interp_levels(numlevels),rmsg
C NCLEND     
       integer   njx,niy,nreqlvs,ripk
       integer   i,j,k,itriv,kupper
       integer   ifound,miy,mjx,isign
       double precision      rlevel,vlev,diff
       double precision      tempout(ew,ns),tmpvlev
       double precision      vcp1,vcp0,valp0,valp1
       double precision      rgas,rgasmd,sclht,ussalr,cvc,eps
       double precision      qvlhsl,ttlhsl,vclhsl,vctophsl
       double precision      wrf_intrp_value
       double precision      plhsl,zlhsl,ezlhsl,tlhsl,psurf,pratio,tlev
       double precision      ezsurf,psurfsm,zsurf,qvapor,vt
       double precision      rconst,expon,exponi
       double precision      ezlev,plev,zlev,ptarget,dpmin,dp
       double precision      pbot,zbot,tbotextrap,e
       double precision      tlclc1,tlclc2,tlclc3,tlclc4
       double precision      thtecon1,thtecon2,thtecon3
       double precision      tlcl,gamma,cp,cpmd,gammamd,gammam
       character cvcord*1

       rgas    = 287.04d0     !J/K/kg
       rgasmd  = .608d0
       ussalr  = .0065d0      ! deg C per m
       sclht   = rgas*256.d0/9.81d0
       eps     = 0.622d0
       rconst  = -9.81d0/(rgas * ussalr) 
       expon   =  rgas*ussalr/9.81d0
       exponi  =  1./expon
       tlclc1   = 2840.d0
       tlclc2   = 3.5d0
       tlclc3   = 4.805d0
       tlclc4   = 55.d0
       thtecon1 = 3376.d0 ! K
       thtecon2 = 2.54d0
       thtecon3 = 0.81d0
       cp       = 1004.d0
       cpmd     = 0.887d0
       gamma    = rgas/cp
       gammamd  = rgasmd-cpmd

       if(vcor .eq. 1) then
          cvcord = 'p'
       else if((vcor .eq. 2) .or. (vcor .eq. 3)) then
          cvcord = 'z'
       else if((vcor .eq. 4) .or. (vcor .eq. 5)) then
          cvcord = 't'
       end if


       miy = ns 
       mjx = ew
       njx = ew
       niy = ns


       do j = 1,mjx
          do i = 1,miy
               tempout(j,i) = rmsg
          end do
       end do



      do nreqlvs = 1,numlevels
         if(cvcord .eq. 'z') then
!Convert rlevel to meters from km

            rlevel = interp_levels(nreqlvs) * 1000.d0
            vlev = exp(-rlevel/sclht)              
         else if(cvcord .eq. 'p') then             
            vlev = interp_levels(nreqlvs)          
         else if(cvcord .eq. 't') then             
            vlev = interp_levels(nreqlvs)        
         end if                                    
 

         do j=1,mjx
         do i=1,miy
cGet the interpolated value that is within the model domain
            ifound = 0
            do k = 1,nz-1
               ripk  = nz-k+1
               vcp1  = vcarray(j,i,ripk-1)
               vcp0  = vcarray(j,i,ripk)
               valp0 = datain(j,i,ripk)
               valp1 = datain(j,i,ripk-1)
               if ((vlev.ge.vcp0.and.vlev.le.vcp1) .or. 
     &            (vlev.le.vcp0.and.vlev.ge.vcp1)) then
c                  print *,i,j,valp0,valp1
                  if((valp0 .eq. rmsg).or.(valp1 .eq. rmsg)) then
                     tempout(j,i) = rmsg
                     ifound=1
                  else
                     if(logp .eq. 1) then
                       vcp1  = log(vcp1)
                       vcp0  = log(vcp0)
                       if(vlev .eq. 0.0d0) then
                         print *,"Pressure value = 0"
                         print *,"Unable to take log of 0"
                         stop
                       end if
                       tmpvlev  = log(vlev)
                     else 
                       tmpvlev = vlev
                     end if
                     tempout(j,i) = wrf_intrp_value(valp0,valp1,
     &                                  tmpvlev,vcp0,vcp1,icase)
c                     print *,"one ",i,j,tempout(j,i)
                     ifound=1
                  end if
                  goto 115
               end if
             end do !end for the k loop
 115  continue


      if (ifound.eq.1) then !Grid point is in the model domain
          goto 333
      end if

cIf the user has requested no extrapolatin then just assign
call values above or below the model level to rmsg.
      if(extrap .eq. 0) then
         tempout(j,i) = rmsg
         goto 333  
      end if


c The grid point is either above or below the model domain    
c
c First we will check to see if the grid point is above the
c model domain.
      vclhsl   = vcarray(j,i,1) !lowest model level
      vctophsl = vcarray(j,i,nz)!highest model level
      diff     = vctophsl-vclhsl
      isign    = nint(diff/abs(diff))
C
      if(isign*vlev.ge.isign*vctophsl) then
C Assign the highest model level to the out array
         tempout(j,i)=datain(j,i,nz)  
C         print *,"at warn",j,i,tempout(j,i)
         goto 333
      endif
                           

c
c   Only remaining possibility is that the specified level is below
c   lowest model level.  If lowest model level value is missing,
c   set interpolated value to missing.
c
      if (datain(i,j,1) .eq. rmsg) then
          tempout(j,i) = rmsg
          goto 333
      endif

c
c   If the field comming in is not a pressure,temperature or height 
C   field we can set the output array to the value at the lowest 
c   model level.
c
      tempout(j,i) = datain(j,i,1)      
c
c   For the special cases of pressure on height levels or height on
c   pressure levels, or temperature-related variables on pressure or
c   height levels, perform a special extrapolation based on
c   US Standard Atmosphere.  Here we calcualate the surface pressure
c   with the altimeter equation.  This is how RIP calculates the 
c   surface pressure.
c
       if (icase.gt.0) then
           plhsl  = pres(j,i,1) * 0.01d0  !pressure at lowest model level
           zlhsl  = ght(j,i,1)            !grid point height a lowest model level
           ezlhsl = exp(-zlhsl/sclht)
           tlhsl  = tk(j,i,1)             !temperature in K at lowest model level
           zsurf  = terrain(j,i)
           qvapor = max((qvp(j,i,1)*.001d0),1.e-15)  
c virtual temperature
c          vt     = tlhsl * (eps + qvapor)/(eps*(1.0 + qvapor)) 
c           psurf  = plhsl * (vt/(vt+ussalr * (zlhsl-zsurf)))**rconst
           psurf    = sfp(j,i)
           psurfsm  = smsfp(j,i) 
           ezsurf   = exp(-zsurf/sclht)

cThe if for checking above ground
           if ((cvcord.eq.'z'.and.vlev.lt.ezsurf).or.
     &         (cvcord.eq.'p'.and.vlev.lt.psurf)) then
c
c      We are below the lowest data level but above the ground.
c      Use linear interpolation (linear in prs and exp-height).
c
               if (cvcord.eq.'p') then
                   plev=vlev
                   ezlev=((plev-plhsl)*ezsurf+(psurf-plev)*ezlhsl)/
     &                     (psurf-plhsl)
                   zlev=-sclht*log(ezlev)
                   if (icase .eq. 2) then
                       tempout(j,i)=zlev
                       goto 333
                   endif

                elseif (cvcord.eq.'z') then
                   ezlev=vlev
                   zlev=-sclht*log(ezlev)
                   plev=((ezlev-ezlhsl)*psurf+(ezsurf-ezlev)*plhsl)/
     &                    (ezsurf-ezlhsl)
                   if (icase .eq. 1) then
                      tempout(j,i)=plev
                      goto 333
                  endif
                endif

            else   !else for checking above ground
                ptarget=psurfsm-150.d0
                dpmin=1.e4
                do k=1,nz
                   ripk = nz-k+1
                   dp=abs((pres(j,i,ripk) * 0.01d0)-ptarget)
                   if (dp.gt.dpmin) goto 334
                   dpmin=min(dpmin,dp)
                enddo
 334            kupper=k-1

                ripk       = nz - kupper + 1
                pbot       = max(plhsl,psurf)
                zbot       = min(zlhsl,zsurf)
                pratio     = pbot/(pres(j,i,ripk) * 0.01d0)
                tbotextrap = tk(j,i,ripk)*(pratio)**expon
c virtual temperature
                vt = tbotextrap * (eps + qvapor)/(eps*(1.0d0+qvapor)) 
                if (cvcord.eq.'p') then
                   plev=vlev
                   zlev=zbot+vt/ussalr*(1.-(vlev/pbot)**expon)
                   if(icase .eq. 2) then
                      tempout(j,i)=zlev
                      goto 333
                   endif
                elseif (cvcord.eq.'z') then
                        zlev=-sclht*log(vlev)
                        plev=pbot*(1.+ussalr/vt*(zbot-zlev))**exponi
                        if (icase .eq. 1) then
                            tempout(j,i)=plev
                            goto 333
                        endif
                 endif                
            end if !end if for checking above ground
       end if !for icase gt 0
        

       if(icase .gt. 2) then !extrapolation for temperature
           tlev=tlhsl+(zlhsl-zlev)*ussalr
           qvapor = max(qvp(j,i,1),1.e-15)
           gammam = gamma*(1.+gammamd*qvapor)
           if(icase .eq. 3) then
              tempout(j,i) = tlev - 273.16d0
           else if(icase .eq. 4) then
              tempout(j,i) = tlev
C Potential temperature - theta
           else if (icase. eq. 5) then 
               tempout(j,i)=tlev*(1000.d0/plev)**gammam
C extraolation for equivalent potential temperature
           else if (icase .eq. 6) then 
               e    = qvapor*plev/(eps+qvapor)
               tlcl = tlclc1/(log(tlev**tlclc2/e)-tlclc3)+tlclc4
               tempout(j,i)=tlev*(1000.d0/plev)**(gammam)*
     &         exp((thtecon1/tlcl-thtecon2)*qvapor*
     &             (1.+thtecon3*qvapor))   
           end if            
       end if
       
 333  continue          

         end do
         end do
 !        print *,"----done----",interp_levels(nreqlvs)
          do i = 1,njx
             do j = 1,niy
                dataout(i,j,nreqlvs) = tempout(i,j)
           end do
          end do
       end do !end for the nreqlvs
       return
       end  !wrf_vinterp
