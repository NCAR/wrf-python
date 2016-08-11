C NCLFORTSTART                                                                    
      subroutine cloud_frac(pres,rh,lowc,midc,highc,nz,ns,ew)

      implicit none
      integer  nz,ns,ew
      real     pres(ew,ns,nz),rh(ew,ns,nz)
      real     lowc(ew,ns),midc(ew,ns),highc(ew,ns)
C NCLEND

      integer i,j,k
      integer kchi,kcmi,kclo 


      DO j = 1,ns
      DO i = 1,ew
         DO k = 1,nz-1

c          if((pres(i,j,k) .ge. 45000. ) .and.
c     &        (pres(i,j,k) .lt. 80000.))  then
c              kchi = k             

c          else if((pres(i,j,k) .ge. 80000.) .and.
c     &        (pres(i,j,k) .lt. 97000.)) then
c              kcmi = k

c         else if (pres(i,j,k) .ge. 97000.) then 
c              kclo = k
c         end if
          IF ( pres(i,j,k) .gt. 97000. ) kclo=k
          IF ( pres(i,j,k) .gt. 80000. ) kcmi=k
          IF ( pres(i,j,k) .gt. 45000. ) kchi=k
   
        end do

        DO k = 1,nz-1
          IF ( k .ge. kclo .AND. k .lt. kcmi ) then          
               lowc(i,j) = AMAX1(rh(i,j,k),lowc(i,j))
          else IF ( k .ge. kcmi .AND. k .lt. kchi ) then              !! mid cloud
              midc(i,j) = AMAX1(rh(i,j,k),midc(i,j))
          else if ( k .ge. kchi )  then                               !! high cloud
              highc(i,j) = AMAX1(rh(i,j,k),highc(i,j)) 
          end if
        END DO


        lowc(i,j)  = 4.0 * lowc(i,j)/100.-3.0
        midc(i,j)  = 4.0 * midc(i,j)/100.-3.0
        highc(i,j) = 2.5 * highc(i,j)/100.-1.5

       lowc(i,j)  = amin1(lowc(i,j),1.0)
       lowc(i,j)  = amax1(lowc(i,j),0.0)
       midc(i,j)  = amin1(midc(i,j),1.0)
       midc(i,j)  = amax1(midc(i,j),0.0)
       highc(i,j) = amin1(highc(i,j),1.0)
       highc(i,j) = amax1(highc(i,j),0.0)

       END DO
       END DO
       return
       end
