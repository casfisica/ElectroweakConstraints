*================================================================================* 
*                                                                                *
*     Sub rutinas para modificar el archivo fit.f y tener en cuenta el segundo   * 
*     orden en los coeficientes de Wilson                                        * 
*                                                                                * 
*================================================================================*



*--------------------------------------------------------------------------------* 
*                                                                                * 
*     Sub rutina para llamar Vmatrix, y organizar su salida en estructura        * 
*                                                                                * 
*--------------------------------------------------------------------------------* 


      subroutine RotMatrix(debug)
      implicit none
      logical debug             ! If is .true. use debuguin part
C      double complex VuL(3,3),VuR(3,3),VdL(3,3),VdR(3,3),VnL(3,3)
C      double complex VnR(3,3),VeL(3,3),VeR(3,3) !Deben de ir en comondelta      #
C      double complex Vmat(4,2,3,3)
      integer P,II,i,j
      include 'commondelta.f'   !the out Vmat is inside of commondelta.f

*     Import the values of the rotation matrix
      call Vmatrix(VuL,VuR,VdL,VdR,VnL,VnR,VeL,VeR)

*     Organize the values of rotation matrix as an estructure Vmat(P,I,i,j)
*     where P: Fermion type; 1:UP, 2:DOWN, 3:Neutrino, 4: Electron
*     I: Chirality; 1:Left, 2:Right
*     i,j: Flavor or Families; 1:electron, 2:muon, 3:Tau
      
      do 100, P=1, 4            !Fermion type
         do 200, II=1, 2        !Chirality
            do 300, i=1, 3
               do 400, j=1,3
                  if (P.eq.1) then
                     if (II.eq.1) then
                        Vmat(P,II,i,j)= VuL(i,j)
                     else if (II.eq.2) then
                        Vmat(P,II,i,j)= VuR(i,j)
                     else
                        STOP 'errVmat'
                     end if     !end I if
                  else if (p.eq.2) then
                     if (II.eq.1) then
                        Vmat(P,II,i,j)= VdL(i,j)
                     else if (II.eq.2) then
                        Vmat(P,II,i,j)= VdR(i,j)
                     else
                        STOP 'errVmat'
                     end if     !end I if
                  else if (P.eq.3) then
                     if (II.eq.1) then
                        Vmat(P,II,i,j)= VnL(i,j)
                     else if (II.eq.2) then
                        Vmat(P,II,i,j)= VnR(i,j)
                     else
                        STOP 'errVmat'
                     end if     !end I if
                  else if (P .eq. 4) then
                     if (II.eq.1) then
                        Vmat(P,II,i,j)= VeL(i,j)
                     else if (II.eq.2) then
                        Vmat(P,II,i,j)= VeR(i,j)
                     else
                        STOP 'errVmat'
                     end if     !end I if
                  else
                     STOP 'errVmat'
                  end if        !end P if
 400           continue
 300        continue
 200     continue
 100  continue
      
         
************************************DEBUGGING************************************   
         
         if (debug) then
            open (unit=10,file='RotMatrix.out',status='new')
 2000       format('',A5,4I2,A2,2E15.7)
          
            do 500, p=1, 4
               do 600, II=1, 2
                  do 700, i=1, 3
                     do 800, j=1, 3
                        if (Vmat(P,II,i,j).ne.0d0) then
                           write(10,2000) "Vmat(",P,II,i,j,")=",
     .                          Vmat(P,II,i,j)
                        endif
 800                 continue
 700              continue
 600           continue
 500        continue

         
         
      end if                    !End debugguing if
      close(unit=10)

      return
      end                       !End subroutine RotMatrix

**********************************************************************************


      
*--------------------------------------------------------------------------------* 
*                                                                                * 
*   Subrutine Epsilon(epsi,gcop,debug), call zprime and organize on a estructure * 
*                                                                                * 
*--------------------------------------------------------------------------------*


      subroutine Epsilon(debug)
      implicit none
      logical debug             ! If is .true. use debuguin part
C      double complex epsi(4,2,3), gcop(4,2,3)
      integer P,II,i,j
      double precision xval(27)
      
      include 'common.f'
      include 'commondelta.f'   !The outs epsi,gcop, are inside

C     Import the coupling values of Z' to fermions
      call zprimecoup(xval,27)  

      do 100, P=1, 4
         do 200, II=1, 2
            do 300, i=1, 3
               if (P.eq.1) then
                  if (II.eq.1) then
                     if (i.eq.1) then
                        epsi(P,II,i) = eps2_L(4) !UP Left
                     else if (i.eq.2) then
                        epsi(P,II,i) = eps2_L(5) !Charm Left
                     else if (i.eq.3) then
                        epsi(P,II,i) = eps2_L(6) !Top Left
                     else
                        STOP 'errEps'
                     end if
                  else if (II.eq.2) then
                     if (i.eq.1) then
                        epsi(P,II,i) = eps2_R(4) !UP Right
                     else if (i.eq.2) then
                        epsi(P,II,i) = eps2_R(5) !Charm Right
                     else if (i.eq.3) then
                        epsi(P,II,i) = eps2_R(6) !Top Right
                     else
                        STOP 'errEps'
                     end if
                  else
                     STOP 'errEps'
                  end if
               else if (P.eq.2) then
                  if (II.eq.1) then
                     if (i.eq.1) then
                        epsi(P,II,i) = eps2_L(7) !Down Left
                     else if (i.eq.2) then           
                        epsi(P,II,i) = eps2_L(8) !Strange Left
                     else if (i.eq.3) then           
                        epsi(P,II,i) = eps2_L(9) !Bottom Left
                     else                            
                        STOP 'errEps'      
                     end if                          
                  else if (II.eq.2) then             
                     if (i.eq.1) then                
                        epsi(P,II,i) = eps2_R(7) !Down Right
                     else if (i.eq.2) then           
                        epsi(P,II,i) = eps2_R(8) !Strange Right
                     else if (i.eq.3) then           
                        epsi(P,II,i) = eps2_R(9) !Bottom Right
                     else                            
                        STOP 'errEps'      
                     end if
                  else
                     STOP 'errEps'
                  end if
               else if (P.eq.3) then
                  if (II.eq.1) then
                     if (i.eq.1) then
                        epsi(P,II,i) = eps2_L(0)  !Neutrino e- Left
                     else if (i.eq.2) then           
                        epsi(P,II,i) = eps2_L(-1) !Neutrino M Left
                     else if (i.eq.3) then           
                        epsi(P,II,i) = eps2_L(-2) !Neutrino T Left
                     else                            
                        STOP 'errEps'      
                     end if                          
                  else if (II.eq.2) then             
                     if (i.eq.1) then                
                        epsi(P,II,i) = eps2_R(0)  !Neutrino e- Right
                     else if (i.eq.2) then           
                        epsi(P,II,i) = eps2_R(-1) !Neutrino M Right
                     else if (i.eq.3) then           
                        epsi(P,II,i) = eps2_R(-2) !Neutrino T Right
                     else                            
                        STOP 'errEps'      
                     end if
                  else
                     STOP 'errEps'
                  end if
               else if (P.eq.4) then
                  if (II.eq.1) then
                     if (i.eq.1) then
                        epsi(P,II,i) = eps2_L(1)  ! e- Left
                     else if (i.eq.2) then           
                        epsi(P,II,i) = eps2_L(2) ! Mu Left
                     else if (i.eq.3) then           
                        epsi(P,II,i) = eps2_L(3) ! Tau Left
                     else                            
                        STOP 'errEps'      
                     end if                          
                  else if (II.eq.2) then             
                     if (i.eq.1) then                
                        epsi(P,II,i) = eps2_R(1)  ! e- Right
                     else if (i.eq.2) then           
                        epsi(P,II,i) = eps2_R(2) ! Mu Right
                     else if (i.eq.3) then           
                        epsi(P,II,i) = eps2_R(3) !Tau Right
                     else                             
                        STOP 'errEps'       
                     end if
                  else
                     STOP 'errEps'
                  end if
               else
                  STOP 'errEps'   
               end if
 300        continue
 200     continue
 100  continue

        
         do 2400, P=1, 4
            do 2500, II=1, 2    !Para el caso de g 1:Axial, 2:Vectorial
               do 2600, i=1, 3
                  if (II.eq.1) then
                     gcop(P,II,i)= epsi(P,1,i)+
     .                    epsi(P,2,i)
                  else if (II.eq.2) then
                     gcop(P,II,i)= epsi(P,1,i)-
     .                    epsi(P,2,i)
                  else
                  end if
 2600          continue
 2500       continue
 2400    continue   
      


************************************DEBUGGING************************************   
         
      if (debug) then
         open (unit=11,file='Epsilon.out',status='new')
 2000    format('',A5,3I2,A2,2E15.7)

         do 400, P=1, 4
            do 500, II=1, 2
               do 600, i=1, 3
                  if (epsi(P,II,i).ne.0d0) then
                     write(11,2000) "epsi(",P,II,i,")=",epsi(P,II,i)
                  endif
 600           continue
 500        continue
 400     continue

         write(11,*) "-------------------------------------------"
         
         do 410, P=1, 4
            do 510, II=1, 2
               do 610, i=1, 3
                  if (gcop(P,II,i).ne.0d0) then
                     write(11,2000) "gcop(",P,II,i,")=",gcop(P,II,i)
                  endif
 610           continue
 510        continue
 410     continue
                  
      end if                    !End debugguing if
      close(unit=11)
      
      return
      end                       !End subroutine RotMatrix

**********************************************************************************


      
*--------------------------------------------------------------------------------* 
*                                                                                * 
*           Subrutine VPVt(VPT,debug) Calculate the product VPV^t                * 
*                                                                                * 
*--------------------------------------------------------------------------------*


      subroutine VPVt(debug)
      implicit none
C      double complex VPV(3,4,2,3,3),Vmat(4,2,3,3) !debe de ir en commondelta     #
      integer betha,i,j,P,II
      logical debug             ! If is .true. use debuguin part
      include 'commondelta.f'   ! The out VPV, is inside commondelta
      
*     Import the values of Vmatrix(rotation matrix)
C      call RotMatrix(Vmat,debug) !Only use in construction, the call is outside 

      do 500, betha=1, 3
         do 100, P=1, 4
            do 200, II=1, 2
               do 300, i=1, 3
                  do 400, j=1, 3
                     VPV(betha,P,II,i,j)= Vmat(P,II,i,betha)*
     .                    dconjg(Vmat(P,II,j,betha))
 400              continue
 300           continue
 200        continue
 100     continue
 500  continue
      

************************************DEBUGGING************************************   
         
      if (debug) then
         open (unit=12,file='VPV.out',status='new')
 2000    format('',A4,5I2,A2,2E15.7)
         
         do 1000, betha=1, 3
            do 600, P=1, 4
               do 700, II=1, 2
                  do 800, i=1, 3
                     do 900, j=1, 3
                        if (VPV(betha,P,II,i,j).ne.0d0) then
                           write(12,2000) "VPV(",betha,P,II,i,j,")=",
     .                          VPV(betha,P,II,i,j)
                        endif
 900                 continue
 800              continue
 700           continue
 600        continue
 1000    continue

         
 2001    format('',A5,4I2,A2,2E15.7)

         write(12,*) "-----------------------------------------"
         
            do 5001, p=1, 4
               do 6001, II=1, 2
                  do 7001, i=1, 3
                     do 8001, j=1, 3
                        if (Vmat(P,II,i,j).ne.0d0) then
                           write(12,2001) "Vmat(",P,II,i,j,")=",
     .                          Vmat(P,II,i,j)
                        endif
 8001                continue
 7001             continue
 6001          continue
 5001       continue

         
      end if                    !End debugguing if
      
      close(unit=12)
      
      return
      end                       !End subroutine VPVt

**********************************************************************************




      
*--------------------------------------------------------------------------------* 
*                                                                                * 
*                         Subrutine Delt(epsi,debug)                             * 
*                                                                                * 
*--------------------------------------------------------------------------------*
 
      subroutine Delt(debug)
      implicit none
C      double complex epsi(4,2,3), Del(3,3,4,2), gcop(4,2,3) !commondelta
      integer P,II,B,a,i
      logical debug             ! If is .true. use debuguin part
      include 'commondelta.f'   !The out Del, is inside
      
C      call Epsilon(epsi,gcop,debug) !Only use in construction, the call is outside

      do 100, B=1, 3
         do 200, a=1, 3
            do 300, P=1, 4
               do 400, II=1, 2
                  Del(B,a,P,II)=epsi(P,II,B)-epsi(P,II,a)
 400           continue
 300        continue
 200     continue
 100  continue
      
      
 
************************************DEBUGGING************************************  
         
      if (debug) then
         open (unit=20,file='Del.out',status='new')
*     Format $A#:character, $: number of entries of the type
*     #: number of spaces; I:integer; D: doubles; F
 1000    format('',A4,4I2,A2,2E15.7) 
         
         do 500, B=1, 3
            do 600, a=1, 3
               do 700, P=1, 4
                  do 800, II=1, 2
                     if (Del(B,a,P,II).ne.0d0) then
                        write(20,1000) "Del(",B,a,P,II,")=",
     .                       Del(B,a,P,II)
                     endif
 800              continue
 700           continue
 600        continue
 500     continue
                  
         write(20,*) "-------------------------------------------"
 
 2001    format('',A5,3I2,A2,2E15.7)

         do 4001, P=1, 4
            do 5001, II=1, 2
               do 6001, i=1, 3
                  if (epsi(P,II,i).ne.0d0) then
                     write(20,2001) "epsi(",P,II,i,")=",epsi(P,II,i)
                  endif
 6001           continue
 5001        continue
 4001     continue

         write(20,*) "-------------------------------------------"
         
         do 4101, P=1, 4
            do 5101, II=1, 2
               do 6101, i=1, 3
                  if (gcop(P,II,i).ne.0d0) then
                     write(20,2001) "gcop(",P,II,i,")=",gcop(P,II,i)
                  endif
 6101           continue
 5101        continue
 4101     continue


      end if                    !End debugguing if
      
      close(unit=20)
      
      return
      end                       !End subroutine Delt
 
**********************************************************************************

      
      
*--------------------------------------------------------------------------------* 
*                                                                                * 
*                           Subrutine DBDG                                       * 
*                                                                                * 
*--------------------------------------------------------------------------------*
 
      subroutine DBDG(debug)
      implicit none
C      double complex VPV(3,4,2,3,3), Del(3,3,4,2) !commndelta
C      double complex DB(3,3,4,2,3,3), DG(3,3,4,2,3,3) !Commondelta
      integer P,II,B,a,i,j,betha
      logical debug             ! If is .true. use debuguin part
      include 'commondelta.f'   !The outs DB,DG, are inside od commondelta.f 

C      call VPVt(VPV,debug)     !The call is made outside
C      call Delt(Del,debug)    


      do 100, B=1, 3
         do 200, a=1, 3
            do 300, P=1, 4
               do 400, II=1, 2
                  do 500, i=1, 3
                     do 600, j=1, 3
                        DB(B,a,P,II,i,j)= Del(B,a,P,II)*VPV(B,P,II,i,j)
 600                 continue
 500              continue
 400           continue
 300        continue
 200     continue
 100  continue

C     1: Vectorial, 2: Axial

      do 110, B=1, 3
         do 210, a=1, 3
            do 310, P=1, 4
               do 510, i=1, 3
                  do 610, j=1, 3
                     DG(B,a,P,1,i,j)=DB(B,a,P,1,i,j)+DB(B,a,P,2,i,j)
                     DG(B,a,P,2,i,j)=DB(B,a,P,1,i,j)-DB(B,a,P,2,i,j)
 610              continue
 510           continue
 310        continue
 210     continue
 110  continue
 

      
 
************************************DEBUGGING************************************  
         
      if (debug) then
         open (unit=14,file='DBDG.out',status='new')

 2000    format('',A3,6I2,A2,2E15.7)
         
      do 101, B=1, 3
         do 201, a=1, 3
            do 301, P=1, 4
               do 401, II=1, 2
                  do 501, i=1, 3
                     do 601, j=1, 3
                        if (DB(B,a,P,II,i,j).ne.0d0) then
                           write(14,2000) "DB(",B,a,P,II,i,j,")=",
     .                          DB(B,a,P,II,i,j)
                        endif
 601                 continue
 501              continue
 401           continue
 301        continue
 201     continue
 101  continue

      write(14,*)"-------------------------------------------"
      
      do 111, B=1, 3
         do 211, a=1, 3
            do 311, P=1, 4
               do 411, II=1, 2
                  do 511, i=1, 3
                     do 611, j=1, 3
                        if (DG(B,a,P,II,i,j).ne.0d0) then
                           write(14,2000) "DG(",B,a,P,II,i,j,")=",
     .                          DG(B,a,P,II,i,j)
                        endif
 611                 continue
 511              continue
 411           continue
 311        continue
 211     continue
 111  continue

      write(14,*) "-------------------------------------------"

      
 2001 format('',A4,5I2,A2,2E15.7)
         
         do 1001, betha=1, 3
            do 6001, P=1, 4
               do 7001, II=1, 2
                  do 8001, i=1, 3
                     do 9001, j=1, 3
                        if (VPV(betha,P,II,i,j).ne.0d0) then
                           write(14,2001) "VPV(",betha,P,II,i,j,")=",
     .                          VPV(betha,P,II,i,j)
                        endif
 9001                continue
 8001             continue
 7001          continue
 6001       continue
 1001    continue

         write(14,*) "-------------------------------------------"


 1000    format('',A4,4I2,A2,2E15.7) 
         
         do 5100, B=1, 3
            do 6100, a=1, 3
               do 7100, P=1, 4
                  do 8100, II=1, 2
                     if (Del(B,a,P,II).ne.0d0) then
                        write(14,1000) "Del(",B,a,P,II,")=",
     .                       Del(B,a,P,II)
                     endif
 8100             continue
 7100          continue
 6100       continue
 5100    continue
      
      end if                    !End debugguing if
      close(unit=14)
      
      return
      end                       !End subroutine RotMatrix
 
**********************************************************************************


*--------------------------------------------------------------------------------* 
*                                                                                * 
*                           Subrutine ZSM                                        * 
*                                                                                * 
*--------------------------------------------------------------------------------*

      subroutine ZSM(debug)
      dimension eps1_L(0:9),eps1_R(0:9), v(0:9),a(0:9) !Auxiliary until cuple
      integer P,II,i,j,f      
C      double complex epsi1(4,2,3), gcop1(4,2,3) !commondelta
      logical debug
      include 'commondelta.f'   !epsi1,gcop1, are inside
      
*******************************************************************************
*                                     Temporal                                *
*******************************************************************************     
      do 220 f = 0, 9                                                         ! 
         v(f) = (2*(f) -1)/2                                                  ! 
         a(f) =  f                                                            ! 
 220  continue                                                                ! 
*******************************************************************************
      
      do 197, j=0, 9
         eps1_L(j) = (v(j) + a(j))/2
         eps1_R(j) = (v(j) - a(j))/2
 197  continue

      do 100, P=1, 4
         do 200, II=1, 2
            do 300, i=1, 3
               if (P.eq.1) then
                  if (II.eq.1) then
                     if (i.eq.1) then
                        epsi1(P,II,i)= eps1_L(4) !UP Left
                     else if (i.eq.2) then
                        epsi1(P,II,i)= eps1_L(5) !Charm Left
                     else if (i.eq.3) then
                        epsi1(P,II,i)= eps1_L(6) !Top Left
                     else
                        STOP 'errEps'
                     end if
                  else if (II.eq.2) then
                     if (i.eq.1) then
                        epsi1(P,II,i)= eps1_R(4) !UP Right
                     else if (i.eq.2) then
                        epsi1(P,II,i)= eps1_R(5) !Charm Right
                     else if (i.eq.3) then
                        epsi1(P,II,i)= eps1_R(6) !Top Right
                     else
                        STOP 'errEps'
                     end if
                  else
                     STOP 'errEps'
                  end if
               else if (P.eq.2) then
                  if (II.eq.1) then
                     if (i.eq.1) then
                        epsi1(P,II,i)= eps1_L(7) !Down Left
                     else if (i.eq.2) then           
                        epsi1(P,II,i)= eps1_L(8) !Strange Left
                     else if (i.eq.3) then           
                        epsi1(P,II,i)= eps1_L(9) !Bottom Left
                     else                            
                        STOP 'errEps'      
                     end if                          
                  else if (II.eq.2) then             
                     if (i.eq.1) then                
                        epsi1(P,II,i)= eps1_R(7) !Down Right
                     else if (i.eq.2) then           
                        epsi1(P,II,i)= eps1_R(8) !Strange Right
                     else if (i.eq.3) then           
                        epsi1(P,II,i)= eps1_R(9) !Bottom Right
                     else                            
                        STOP 'errEps'      
                     end if
                  else
                     STOP 'errEps'
                  end if
               else if (P.eq.3) then
                  if (II.eq.1) then
                     if (i.eq.1) then
                        epsi1(P,II,i)= eps1_L(0) !Neutrino e- Left
                     else if (i.eq.2) then           
                        epsi1(P,II,i)= eps1_L(0) !Neutrino M Left
                     else if (i.eq.3) then           
                        epsi1(P,II,i)= eps1_L(0) !Neutrino T Left
                     else                            
                        STOP 'errEps'      
                     end if                          
                  else if (II.eq.2) then             
                     if (i.eq.1) then                
                        epsi1(P,II,i)= eps1_R(0)  !Neutrino e- Right
                     else if (i.eq.2) then           
                        epsi1(P,II,i)= eps1_R(0) !Neutrino M Right
                     else if (i.eq.3) then           
                        epsi1(P,II,i)= eps1_R(0) !Neutrino T Right
                     else                            
                        STOP 'errEps'      
                     end if
                  else
                     STOP 'errEps'
                  end if
               else if (P.eq.4) then
                  if (II.eq.1) then
                     if (i.eq.1) then
                        epsi1(P,II,i)= eps1_L(1)  ! e- Left
                     else if (i.eq.2) then           
                        epsi1(P,II,i)= eps1_L(2) ! Mu Left
                     else if (i.eq.3) then           
                        epsi1(P,II,i)= eps1_L(3) ! Tau Left
                     else                            
                        STOP 'errEps'      
                     end if                          
                  else if (II.eq.2) then             
                     if (i.eq.1) then                
                        epsi1(P,II,i)= eps1_R(1)  ! e- Right
                     else if (i.eq.2) then           
                        epsi1(P,II,i)= eps1_R(2) ! Mu Right
                     else if (i.eq.3) then           
                        epsi1(P,II,i)= eps1_R(3) !Tau Right
                     else                             
                        STOP 'errEps'       
                     end if
                  else
                     STOP 'errEps'
                  end if
               else
                  STOP 'errEps'   
               end if
 300        continue
 200     continue
 100  continue

      do 111, P=1, 4
         do 211, II=1, 2
            do 311, i=1, 3
               if (II.eq.1) then
                  gcop1(P,II,i)=epsi1(P,1,i)+epsi1(P,2,i)
               else if (II.eq.2) then
                  gcop1(P,II,i)=epsi1(P,1,i)-epsi1(P,2,i)
               else
                  STOP 'errEps'
               end if
 311        continue
 211     continue
 111  continue
      
C***********************************DEBUGGING***********************************C   
         
      if (debug) then
         open (unit=15,file='ZSM.out',status='new')

 2000    format('',A6,3I2,A2,2E15.7)
         
         do 400, P=1, 4
            do 500, II=1, 2
               do 600, i=1, 3
                  if (epsi1(P,II,i).ne.0d0) then
                     write(15,2000) "epsi1(",P,II,i,")=",epsi1(P,II,i)
                  endif
 600           continue
 500        continue
 400     continue


         do 410, P=1, 4
            do 510, II=1, 2
               do 610, i=1, 3
                  if (gcop1(P,II,i).ne.0d0) then
                     write(15,2000) "gcop1(",P,II,i,")=",gcop1(P,II,i)
                  endif
 610           continue
 510        continue
 410     continue

      end if                    !End debugguing if
      close(unit=15)
      
      return
      end                       !End subroutine ZSM

C********************************************************************************C

      
      
*--------------------------------------------------------------------------------*
*                                                                                *
*                                  Subrutine DCDCM                               *
*                                                                                *
*--------------------------------------------------------------------------------*
 
      subroutine DCDCM(w,y,debug)
      implicit none
C     double complex DB(3,3,4,2,3,3), DG(3,3,4,2,3,3) !Commondelta 
C     double complex DC(4,4,4,3,3,3,3), DCM(4,4,4,3,3,3,3),gcop1(4,2,3)
C     double complex epsi1(4,2,3),epsi(4,2,3),gcop(4,2,3)
      double complex w,y,deltaij(3,3)
      integer P,II,C,i,j,k,l,a,B
      logical debug             ! If is .true. use debuguin part
      include 'commondelta.f' !DC,DCM, are inside
      
*     II is a double indice, 1: LL(left, left), 2: LR, 3: RL, 4:  
   
      do 2200, i=1, 3
         do 2300, j=1, 3
            if (i.eq.j) then
               deltaij(i,j)= 1d0
            else
               deltaij(i,j)= 0d0
            end if
            
 2300    continue
 2200 continue

*      call DBDG(DB,DG,debug)    !The call it is outside
*      call ZSM(epsi1,gcop1,debug)    


      
      do 100, B=1, 3
         do 200, a=1, 3
            do 300, P=1, 4
               do 400, II=1, 4
                  do 500, i=1, 3
                     do 600, j=1, 3
                        do 700, C=1, 4
                           do 900, k=1, 3
                              do 1000, l=1, 3
                                 if (II.eq.1) then
                                    DC(P,II,C,i,j,k,l)=w*(deltaij(i,j)*
     .                                   epsi1(P,1,a)*DB(B,a,C,1,k,l)+
     .                                   deltaij(k,l)*
     .                                   epsi1(C,1,a)*DB(B,a,P,1,i,j))
     .                                   +y*(deltaij(i,j)*epsi(P,1,a)*
     .                                   DB(B,a,C,1,k,l)+
     .                                   deltaij(k,l)*epsi(C,1,a)*
     .                                   DB(B,a,P,1,i,j)+
     .                                   DB(B,a,P,1,i,j)*
     .                                   DB(B,a,C,1,k,l))
                                 else if (II.eq.2) then
                                    DC(P,II,C,i,j,k,l)=w*(deltaij(i,j)*
     .                                   epsi1(P,1,a)*DB(B,a,C,2,k,l)+
     .                                   deltaij(k,l)*
     .                                   epsi1(C,2,a)*DB(B,a,P,1,i,j))
     .                                   +y*(deltaij(i,j)*epsi(P,1,a)*
     .                                   DB(B,a,C,2,k,l)+
     .                                   deltaij(k,l)*epsi(C,2,a)*
     .                                   DB(B,a,P,1,i,j)+
     .                                   DB(B,a,P,1,i,j)*
     .                                   DB(B,a,C,2,k,l))
                                 else if (II.eq.3) then
                                    DC(P,II,C,i,j,k,l)=w*(deltaij(i,j)*
     .                                   epsi1(P,2,a)*DB(B,a,C,1,k,l)+
     .                                   deltaij(k,l)*
     .                                   epsi1(C,1,a)*DB(B,a,P,2,i,j))
     .                                   +y*(deltaij(i,j)*epsi(P,2,a)*
     .                                   DB(B,a,C,1,k,l)+
     .                                   deltaij(k,l)*epsi(C,1,a)*
     .                                   DB(B,a,P,2,i,j)+
     .                                   DB(B,a,P,2,i,j)*
     .                                   DB(B,a,C,1,k,l))
                                 else if (II.eq.4) then
                                    DC(P,II,C,i,j,k,l)=w*(deltaij(i,j)*
     .                                   epsi1(P,2,a)*DB(B,a,C,2,k,l)+
     .                                   deltaij(k,l)*
     .                                   epsi1(C,2,a)*DB(B,a,P,2,i,j))
     .                                   +y*(deltaij(i,j)*epsi(P,2,a)*
     .                                   DB(B,a,C,2,k,l)+
     .                                   deltaij(k,l)*epsi(C,2,a)*
     .                                   DB(B,a,P,2,i,j)+
     .                                   DB(B,a,P,2,i,j)*
     .                                   DB(B,a,C,2,k,l)) 
                                 else
                                    STOP 'errEps'
                                 end if
 1000                         continue
 900                       continue
 700                    continue
 600                 continue
 500              continue
 400           continue
 300        continue
 200     continue
 100  continue

C     de II indice is a double one, next it is the mining
C     II=1,2,3,4 1: Vectorial,Vectorial; 2: Vectorial,Axial; 3: A,V; 4: A,A
      
      do 110, B=1, 3
         do 210, a=1, 3
            do 310, P=1, 4
               do 410, II=1, 4
                  do 510, i=1, 3
                     do 610, j=1, 3
                        do 710, C=1, 4
                           do 910, k=1, 3
                              do 1100, l=1, 3
                                 if (II.eq.1) then
                                    DCM(P,II,C,i,j,k,l)=0d0
                                 else if (II.eq.2) then
                                    DCM(P,II,C,i,j,k,l)=0d0
                                 else if (II.eq.3) then
                                    DCM(P,II,C,i,j,k,l)=0d0
                                 else if (II.eq.4) then
                                    DCM(P,II,C,i,j,k,l)=0d0
                                 else
                                    STOP 'errEps'
                                 end if
 1100                         continue
 910                       continue
 710                    continue
 610                 continue
 510              continue
 410           continue
 310        continue
 210     continue
 110  continue
 
      
 
************************************DEBUGGING************************************  
         
      if (debug) then
         open (unit=15,file='DCDCM.out',status='new')
         
 2020    format('',A8,2I2,A2,2E15.7)
         
         do 2290, i=1, 3
            do 2390, j=1, 3
               write(15,2020) "deltaij(",i,j,")=",deltaij(i,j) 
 2390       continue
 2290    continue

         
 2000    format('',A3,7I2,A2,2E15.7) 
 2001    format('',A4,7I2,A2,2E15.7)
 
         write(15,*)"------------------DC-----------------------"
         
         do 101, B=1, 3
            do 201, a=1, 3
               do 301, P=1, 4
                  do 401, II=1, 4
                     do 501, i=1, 3
                        do 601, j=1, 3
                           do 701, C=1, 4
                              do 901, k=1, 3
                                 do 1001, l=1, 3
                                    if (DC(P,II,C,i,j,k,l).ne.0d0) then
                                       write(15,2000) "DC(",P,II,C,i,j,
     .                                      k,l,")=", DC(P,II,C,i,j,k,l)
                                    endif
 1001                            continue
 901                          continue
 701                       continue
 601                    continue
 501                 continue
 401              continue
 301           continue
 201        continue
 101     continue

         write(15,*)"--------------------DCM---------------------"
         
         do 111, B=1, 3
            do 211, a=1, 3
               do 311, P=1, 4
                  do 411, II=1, 4
                     do 511, i=1, 3
                        do 611, j=1, 3
                           do 711, C=1, 4
                              do 911, k=1, 3
                                 do 1101, l=1, 3
                                    if (DCM(P,II,C,i,j,k,l).ne.0d0) then 
                                       write(15,2001) "DCM(",P,II,C,i,j,
     .                                      k,l,")=",
     .                                      DCM(P,II,C,i,j,k,l)
                                    endif
 1101                            continue
 911                          continue
 711                       continue
 611                    continue
 511                 continue
 411              continue
 311           continue
 211        continue
 111     continue

         write(15,*) "------------------epsi1-----------------------"
         
 2002    format('',A6,3I2,A2,2E15.7)
         
         do 2400, P=1, 4
            do 2500, II=1, 2
               do 2600, i=1, 3
                  if (epsi1(P,II,i).ne.0d0) then
                     write(15,2002) "epsi1(",P,II,i,")=",epsi1(P,II,i)
                  endif
 2600          continue
 2500       continue
 2400    continue


         do 2410, P=1, 4
            do 2510, II=1, 2
               do 2610, i=1, 3
                  if (gcop1(P,II,i).ne.0d0) then
                     write(15,2002) "gcop1(",P,II,i,")=",gcop1(P,II,i)
                  endif
 2610          continue
 2510       continue
 2410    continue


         write(15,*) "-------------------DB------------------------"

 2003    format('',A3,6I2,A2,2E15.7)
         
      do 3101, B=1, 3
         do 3201, a=1, 3
            do 3301, P=1, 4
               do 3401, II=1, 2
                  do 3501, i=1, 3
                     do 3601, j=1, 3
                        if (DB(B,a,P,II,i,j).ne.0d0) then
                           write(15,2003) "DB(",B,a,P,II,i,j,")=",
     .                          DB(B,a,P,II,i,j)
                     endif
 3601                continue
 3501             continue
 3401          continue
 3301       continue
 3201    continue
 3101 continue

      write(15,*)"------------------DG-------------------------"
      
      do 3111, B=1, 3
         do 3211, a=1, 3
            do 3311, P=1, 4
               do 3411, II=1, 2
                  do 3511, i=1, 3
                     do 3611, j=1, 3
                        if (DG(B,a,P,II,i,j).ne.0d0) then
                           write(15,2003) "DG(",B,a,P,II,i,j,")=",
     .                          DG(B,a,P,II,i,j)
                        endif
 3611                continue
 3511             continue
 3411          continue
 3311       continue
 3211    continue
 3111 continue


      write(15,*)"------------------epsi-------------------------"
      

 2222 format('',A5,3I2,A2,2E15.7)

         do 4300, P=1, 4
            do 5300, II=1, 2
               do 6300, i=1, 3
                  if (epsi(P,II,i).ne.0d0) then
                     write(15,2222) "epsi(",P,II,i,")=",epsi(P,II,i)
                  endif
 6300          continue
 5300       continue
 4300    continue

         write(15,*) "-------------gcop---------------------------"
         
         do 4109, P=1, 4
            do 5109, II=1, 2
               do 6109, i=1, 3
                  if (gcop(P,II,i).ne.0d0) then
                     write(15,2222) "gcop(",P,II,i,")=",gcop(P,II,i)
                  endif
 6109          continue
 5109       continue
 4109    continue
      end if                    !End debugguing if
      close(unit=15)
      
      return
      end                       !End subroutine RotMatrix
 
**********************************************************************************
