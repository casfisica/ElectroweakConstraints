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


      subroutine Epsilon(epsi,gcop,debug)
      implicit none
      logical debug             ! If is .true. use debuguin part
      double complex epsi(4,2,3), gcop(4,2,3)
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

        
         do 2400, P=1, 4
            do 2500, II=1, 2    !Para el caso de g 1:Axial, 2:Vectorial
               do 2600, i=1, 3
                  if (II.eq.1) then
                     gcop(P,II,i)= epsi(P,1,i)-
     .                    epsi(P,2,i)
                  else if (II.eq.2) then
                     gcop(P,II,i)= epsi(P,1,i)-
     .                    epsi(P,2,i)
                  else
                  end if
 2600          continue
 2500       continue
 2400    continue   
      


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

         do 410, P=1, 4
            do 510, II=1, 2
               do 610, i=1, 3
                  write(11,2000) "gcop(",P,II,i,")=",gcop(P,II,i)  
 610           continue
 510        continue
 410     continue
                  
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
C                         Subrutine Delt(epsi,debug)                             C 
C                                                                                C 
C--------------------------------------------------------------------------------C
 
      subroutine Delt(Del,debug)
      implicit none
      double complex epsi(4,2,3), Del(3,3,4,2), gcop(4,2,3)
      integer P,II,B,a
      logical debug             ! If is .true. use debuguin part

      call Epsilon(epsi,gcop,debug)

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
      end                       !End subroutine Delt
 
C********************************************************************************C

      
      
C--------------------------------------------------------------------------------C 
C                                                                                C 
C                           Subrutine DBDG                                       C 
C                                                                                C 
C--------------------------------------------------------------------------------C
 
      subroutine DBDG(DB,DG,debug)
      implicit none
      double complex VPV(3,4,2,3,3), Del(3,3,4,2)
      double complex DB(3,3,4,2,3,3), DG(3,3,4,2,3,3)
      integer P,II,B,a,i,j
      logical debug             ! If is .true. use debuguin part


      call VPVt(VPV,debug)
      call Delt(Del,debug)


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
 

      
 
C***********************************DEBUGGING***********************************C  
         
      if (debug) then
         open (unit=14,file='DBDG.out',status='new')

 2000    format('',A3,6I2,A2,2E15.7)
         

         
      do 101, B=1, 3
         do 201, a=1, 3
            do 301, P=1, 4
               do 401, II=1, 2
                  do 501, i=1, 3
                     do 601, j=1, 3
                        write(14,2000) "DB(",B,a,P,II,i,j,")=",
     .                       DB(B,a,P,II,i,j)
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
                        write(14,2000) "DG(",B,a,P,II,i,j,")=",
     .                       DG(B,a,P,II,i,j)
 611                 continue
 511              continue
 411           continue
 311        continue
 211     continue
 111  continue

      
      end if                    !End debugguing if
      close(unit=14)
      
      return
      end                       !End subroutine RotMatrix
 
C********************************************************************************C


C--------------------------------------------------------------------------------C 
C                                                                                C 
C                           Subrutine ZSM                                        C 
C                                                                                C 
C--------------------------------------------------------------------------------C

      subroutine ZSM(epsi1,gcop1,debug)
      dimension eps1_L(0:9),eps1_R(0:9), v(0:9),a(0:9)
      integer P,II,i,j,f      
      double complex epsi1(4,2,3), gcop1(4,2,3)
      logical debug


*******************************************************************************
C      Temporal
******************************************************************************      
      do 220 f = 0, 9
         v(f) = (2*(f) -1)/2
         a(f) =  f
 220  continue
******************************************************************************
      
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
                  write(15,2000) "epsi1(",P,II,i,")=",epsi1(P,II,i)  
 600           continue
 500        continue
 400     continue


         do 410, P=1, 4
            do 510, II=1, 2
               do 610, i=1, 3
                  write(15,2000) "gcop1(",P,II,i,")=",gcop1(P,II,i)  
 610           continue
 510        continue
 410     continue

      end if                    !End debugguing if
      close(unit=15)
      
      return
      end                       !End subroutine ZSM

C********************************************************************************C

      
      
C--------------------------------------------------------------------------------C 
C                                                                                C 
C   Subrutine DCDCM    C 
C                                                                                C 
C--------------------------------------------------------------------------------C
 
      subroutine DCDCM(w,y,DC,DCM,debug)
      implicit none
      double complex DB(3,3,4,2,3,3), DG(3,3,4,2,3,3), deltaij(3,3) !Commonlib 
      double complex DC(4,4,4,3,3,3,3), DCM(4,4,4,3,3,3,3),gcop1(4,2,3)
      double complex epsi1(4,2,3),w,y,epsi(4,2,3),gcop(4,2,3)
      integer P,II,C,i,j,k,l,a,B
      logical debug             ! If is .true. use debuguin part
 
C      include 'commondelta'
C     II is a double indice, 1: LL(left, left), 2: LR, 3: RL, 4:  
   
      do 2200, i=1, 3
         do 2300, j=1, 3
            if (i.eq.j) then
               deltaij(i,j)= 1d0
            else
               deltaij(i,j)= 0d0
            end if
            
 2300    continue
 2200 continue

      call DBDG(DB,DG,debug)    !Se llaman desde afuera para evitar calcularlos
      call ZSM(epsi1,gcop1,debug)     !muchas veces


      
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
                                    DC(P,II,C,i,j,k,l)= w
     .                                   *deltaij(i,j)*epsi1(P,1,a)
     .                                   *DB(a,B,C,1,k,l)
     .                                   +w*deltaij(k,l)*epsi1(C,1,a)
     .                                   *DB(a,B,P,1,i,j)
     .                                   +y*deltaij(i,j)*epsi(P,1,a)
     .                                   *DB(a,B,C,1,k,l)
     .                                   +y*deltaij(k,l)*epsi(C,1,a)
     .                                   *DB(a,B,P,1,i,j)
     .                                   +Y*DB(a,B,C,1,k,l)
     .                                   *DB(a,B,P,1,i,j)
                                 else if (II.eq.2) then
                                    DC(P,II,C,i,j,k,l)= w
     .                                   *deltaij(i,j)*epsi1(P,1,a)
     .                                   *DB(a,B,C,2,k,l)
     .                                   +w*deltaij(k,l)*epsi1(C,2,a)
     .                                   *DB(a,B,P,1,i,j)
     .                                   +y*deltaij(i,j)*epsi(P,1,a)
     .                                   *DB(a,B,C,2,k,l)
     .                                   +y*deltaij(k,l)*epsi(C,2,a)
     .                                   *DB(a,B,P,1,i,j)
     .                                   +Y*DB(a,B,C,2,k,l)
     .                                   *DB(a,B,P,1,i,j)
                                 else if (II.eq.3) then
                                    DC(P,II,C,i,j,k,l)= w
     .                                   *deltaij(i,j)*epsi1(P,2,a)
     .                                   *DB(a,B,C,1,k,l)
     .                                   +w*deltaij(k,l)*epsi1(C,1,a)
     .                                   *DB(a,B,P,2,i,j)
     .                                   +y*deltaij(i,j)*epsi(P,2,a)
     .                                   *DB(a,B,C,1,k,l)
     .                                   +y*deltaij(k,l)*epsi(C,1,a)
     .                                   *DB(a,B,P,2,i,j)
     .                                   +Y*DB(a,B,C,1,k,l)
     .                                   *DB(a,B,P,2,i,j)
                                 else if (II.eq.4) then
                                    DC(P,II,C,i,j,k,l)= w
     .                                   *deltaij(i,j)*epsi1(P,2,a)
     .                                   *DB(a,B,C,2,k,l)
     .                                   +w*deltaij(k,l)*epsi1(C,2,a)
     .                                   *DB(a,B,P,2,i,j)
     .                                   +y*deltaij(i,j)*epsi(P,2,a)
     .                                   *DB(a,B,C,2,k,l)
     .                                   +y*deltaij(k,l)*epsi(C,2,a)
     .                                   *DB(a,B,P,2,i,j)
     .                                   +Y*DB(a,B,C,2,k,l)
     .                                   *DB(a,B,P,2,i,j)
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
                                    DCM(P,II,C,i,j,k,l)=(1d0/4d0)*
     .                                   (w*deltaij(i,j)*gcop1(P,1,a)*
     .                                   DG(a,B,C,1,k,l)+
     .                                   w*deltaij(k,l)*gcop1(C,1,a)*
     .                                   DG(a,B,P,1,i,j)+
     .                                   y*deltaij(i,j)*gcop(P,1,a)*
     .                                   DG(a,B,C,1,k,l)+
     .                                   y*deltaij(k,l)*gcop(C,1,a)*
     .                                   DG(a,B,P,1,i,j)+
     .                                   y*DG(a,B,P,1,i,j)*
     .                                   DG(a,B,C,1,k,l))
                                 else if (II.eq.2) then
                                    DCM(P,II,C,i,j,k,l)=(-1d0/4d0)*
     .                                   (w*deltaij(i,j)*gcop1(P,1,a)*
     .                                   DG(a,B,C,2,k,l)+
     .                                   w*deltaij(k,l)*gcop1(C,2,a)*
     .                                   DG(a,B,P,1,i,j)+
     .                                   y*deltaij(i,j)*gcop(P,1,a)*
     .                                   DG(a,B,C,2,k,l)+
     .                                   y*deltaij(k,l)*gcop(C,2,a)*
     .                                   DG(a,B,P,1,i,j)+
     .                                   y*DG(a,B,P,1,i,j)*
     .                                   DG(a,B,C,2,k,l))
                                 else if (II.eq.3) then
                                    DCM(P,II,C,i,j,k,l)= (-1d0/4d0)*
     .                                   (w*deltaij(i,j)*gcop1(P,2,a)*
     .                                   DG(a,B,C,1,k,l)+
     .                                   w*deltaij(k,l)*gcop1(C,1,a)*
     .                                   DG(a,B,P,2,i,j)+
     .                                   y*deltaij(i,j)*gcop(P,2,a)*
     .                                   DG(a,B,C,1,k,l)+
     .                                   y*deltaij(k,l)*gcop(C,1,a)*
     .                                   DG(a,B,P,2,i,j)+
     .                                   y*DG(a,B,P,2,i,j)*
     .                                   DG(a,B,C,1,k,l))
                                 else if (II.eq.4) then
                                    DCM(P,II,C,i,j,k,l)= (1d0/4d0)*
     .                                   (w*deltaij(i,j)*gcop1(P,2,a)*
     .                                   DG(a,B,C,2,k,l)+
     .                                   w*deltaij(k,l)*gcop1(C,2,a)*
     .                                   DG(a,B,P,2,i,j)+
     .                                   y*deltaij(i,j)*gcop(P,2,a)*
     .                                   DG(a,B,C,2,k,l)+
     .                                   y*deltaij(k,l)*gcop(C,2,a)*
     .                                   DG(a,B,P,2,i,j)+
     .                                   y*DG(a,B,P,2,i,j)*
     .                                   DG(a,B,C,2,k,l))
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
 
      
 
C***********************************DEBUGGING***********************************C  
         
      if (debug) then
         open (unit=15,file='DCDCM.out',status='new')
 
 2000    format('',A3,7I2,A2,2E15.7) 
 2001    format('',A4,7I2,A2,2E15.7)
         
         do 101, B=1, 3
            do 201, a=1, 3
               do 301, P=1, 4
                  do 401, II=1, 4
                     do 501, i=1, 3
                        do 601, j=1, 3
                           do 701, C=1, 4
                              do 901, k=1, 3
                                 do 1001, l=1, 3
                                    write(15,2000) "DC(",P,II,C,i,j,k,l,
     .                                   ")=", DC(P,II,C,i,j,k,l)
 1001                            continue
 901                          continue
 701                       continue
 601                    continue
 501                 continue
 401              continue
 301           continue
 201        continue
 101     continue

         write(15,*)"-------------------------------------------"
         
         do 111, B=1, 3
            do 211, a=1, 3
               do 311, P=1, 4
                  do 411, II=1, 4
                     do 511, i=1, 3
                        do 611, j=1, 3
                           do 711, C=1, 4
                              do 911, k=1, 3
                                 do 1101, l=1, 3
                                    write(15,2001) "DCM(",P,II,C,i,j,k,l
     .                                   ,")=", DCM(P,II,C,i,j,k,l)
 1101                            continue
 911                          continue
 711                       continue
 611                    continue
 511                 continue
 411              continue
 311           continue
 211        continue
 111     continue
         




      end if                    !End debugguing if
      close(unit=15)
      
      return
      end                       !End subroutine RotMatrix
 
C********************************************************************************C
