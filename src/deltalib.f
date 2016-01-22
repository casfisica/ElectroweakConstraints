C================================================================================C 
C                                                                                C
C     Sub rutinas para modificar el archivo fit.f y tener en cuenta el segundo   C 
C     orden en los coeficientes de Wilson                                        C 
C                                                                                C 
C================================================================================C



C--------------------------------------------------------------------------------C 
C                                                                                C 
C     Sub rutina para llamar Vmatrix, y organizar su salida en estructura        C 
C                                                                                C 
C--------------------------------------------------------------------------------C 


      subroutine RotMatrix(Vmat,debug)
      implicit none
      logical debug             ! If is .true. use debuguin part
      double complex VuL(3,3),VuR(3,3),VdL(3,3),VdR(3,3),VnL(3,3)
      double complex VnR(3,3),VeL(3,3),VeR(3,3) !Deben de ir en comondelta      #
      double complex Vmat(4,2,3,3)
      integer P,II,i,j

C      include 'commondelta'

C     Import the values of the rotation matrix
      call Vmatrix(VuL,VuR,VdL,VdR,VnL,VnR,VeL,VeR)

C     Organize the values of rotation matrix as an estructure Vmat(P,I,i,j)
C     where P: Fermion type; 1:UP, 2:DOWN, 3:Neutrino, 4: Electron
C     I: Chirality; 1:Left, 2:Right
C     i,j: Flavor or Families; 1:electron, 2:muon, 3:Tau
      
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
      
         
C***********************************DEBUGGING***********************************C   
         
         if (debug) then
            open (unit=10,file='RotMatrix.out',status='new')
 2000       format('',A5,4I2,A2,2E15.7)
          
            do 500, p=1, 4
               do 600, II=1, 2
                  do 700, i=1, 3
                     do 800, j=1, 3
                        write(10,2000) "Vmat(",P,II,i,j,")=",
     .                       Vmat(P,II,i,j)            
 800                 continue
 700              continue
 600           continue
 500        continue

         
         
      end if                    !End debugguing if
      close(unit=10)

      return
      end                       !End subroutine RotMatrix

C********************************************************************************C


      
C--------------------------------------------------------------------------------C 
C                                                                                C 
C   Subrutine Epsilon(epsi,debug), call zprime and organize on a estructure      C 
C                                                                                C 
C--------------------------------------------------------------------------------C


      subroutine Epsilon(epsi,debug)
      implicit none
      logical debug             ! If is .true. use debuguin part
      double complex epsi(4,2,3)
      integer P,II,i,j
      double precision xval(27)
      
      include 'common.f'
C      include 'commondelta'

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




C***********************************DEBUGGING***********************************C   
         
      if (debug) then
         open (unit=11,file='Epsilon.out',status='new')
 2000    format('',A5,3I2,A2,2E15.7)

         do 400, P=1, 4
            do 500, II=1, 2
               do 600, i=1, 3
                  write(11,2000) "epsi(",P,II,i,")=",epsi(P,II,i)  
 600           continue
 500        continue
 400     continue
                  
      end if                    !End debugguing if
      close(unit=11)
      
      return
      end                       !End subroutine RotMatrix

C********************************************************************************C


      
C--------------------------------------------------------------------------------C 
C                                                                                C 
C           Subrutine VPVt(VPT,debug) Calculate the product VPV^t                C 
C                                                                                C 
C--------------------------------------------------------------------------------C


      subroutine VPVt(VPV,debug)
      implicit none
      double complex VPV(3,4,2,3,3),Vmat(4,2,3,3) !debe de ir en commondelta     #
      integer betha,i,j,P,II
      logical debug             ! If is .true. use debuguin part

C     Import the values of Vmatrix(rotation matrix)
      call RotMatrix(Vmat,debug)
      
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
      

C***********************************DEBUGGING***********************************C   
         
      if (debug) then
         open (unit=12,file='VPV.out',status='new')
 2000    format('',A4,5I2,A2,2E15.7)
         
         do 1000, betha=1, 3
            do 600, P=1, 4
               do 700, II=1, 2
                  do 800, i=1, 3
                     do 900, j=1, 3
                        write(12,2000) "VPV(",betha,P,II,i,j,")=",
     .                       VPV(betha,P,II,i,j)
 900                 continue
 800              continue
 700           continue
 600        continue
 1000    continue
         
      end if                    !End debugguing if
      
      close(unit=12)
      
      return
      end                       !End subroutine VPVt

C********************************************************************************C




      
C--------------------------------------------------------------------------------C 
C                                                                                C 
C   Subrutine Delt(epsi,debug)    C 
C                                                                                C 
C--------------------------------------------------------------------------------C
 
      subroutine Delt(Del,debug)
      implicit none
      double complex epsi(4,2,3), Del(3,3,4,2)
      integer P,II,B,a
      logical debug             ! If is .true. use debuguin part

      call Epsilon(epsi,debug)

      do 100, B=1, 3
         do 200, a=1, 3
            do 300, P=1, 4
               do 400, II=1, 2
                  Del(B,a,P,II)=epsi(P,II,B)-epsi(P,II,a)
 400           continue
 300        continue
 200     continue
 100  continue
      
      
 
C***********************************DEBUGGING***********************************C  
         
      if (debug) then
         open (unit=13,file='Del.out',status='new')
C     Format $A#:character, $: number of entries of the type
C     #: number of spaces; I:integer; D: doubles; F
 1000    format('',A4,4I2,A2,2E15.7) 


         
         do 500, B=1, 3
            do 600, a=1, 3
               do 700, P=1, 4
                  do 800, II=1, 2
                     write(13,1000) "Del(",B,a,P,II,")=",Del(B,a,P,II)
 800              continue
 700           continue
 600        continue
 500     continue
                  
      end if                    !End debugguing if

      close(unit=13)
      
      return
      end                       !End subroutine RotMatrix
 
C********************************************************************************C

      
      
CC--------------------------------------------------------------------------------C 
CC                                                                                C 
CC   Subrutine     C 
CC                                                                                C 
CC--------------------------------------------------------------------------------C
C 
C 
C 
C 
CC***********************************DEBUGGING***********************************C  
C         
C      if (debug) then
C         open (unit=12,file='.out',status='new')
C 
C         do 400, P=1, 4
C            do 500, II=1, 2
C               do 600, i=1, 3
C                  write(11,*) "epsi(",P,II,i,")=",epsi(P,II,i)  
C 600           continue
C 500        continue
C 400     continue
C                  
C      end if                    !End debugguing if
C      close(unit=11)
C      
C      return
C      end                       !End subroutine RotMatrix
C 
CC********************************************************************************C

      
