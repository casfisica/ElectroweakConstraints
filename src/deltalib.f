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
      double complex VnR(3,3),VeL(3,3),VeR(3,3)
      double complex Vmat(4,2,3,3)
      integer P,II,i,j

C      include 'common.f'
C      include 'commondelta'

C     Import the values of the rotation matrix
      call Vmatrix(VuL,VuR,VdL,VdR,VnL,VnR,VeL,VeR)

C     Organize the values of rotation matrix as an estructure Vmat(P,I,i,j)
C     where P: Fermion type; 1:UP, 2:DOWN, 3:Neutrino, 4: Electron
C     I: Chirality; 1:Left, 2:Right
C     i,j: Flavor or Families; 1:electron, 2:muon, 3:Tau
      
      do 100, P=1, 4            !Fermion type
         do 200, II=1, 2        !Chirality
            if (P.eq.1) then
               if (II.eq.1) then
                  Vmat(P,II,i,j)= VuL(i,j)
               else if (II.eq.2) then
                  Vmat(P,II,i,j)= VuR(i,j)
               else
                  STOP 'errVmat'
               end if           !end I if
            else if (p.eq.2) then
               if (II.eq.1) then
                  Vmat(P,II,i,j)= VdL(i,j)
               else if (II.eq.2) then
                  Vmat(P,II,i,j)= VdR(i,j)
               else
                  STOP 'errVmat'
               end if           !end I if
            else if (P.eq.3) then
               if (II.eq.1) then
                  Vmat(P,II,i,j)= VnL(i,j)
               else if (II.eq.2) then
                  Vmat(P,II,i,j)= VnR(i,j)
               else
                  STOP 'errVmat'
               end if           !end I if
            else if (P .eq. 4) then
               if (II.eq.1) then
                  Vmat(P,II,i,j)= VeL(i,j)
               else if (II.eq.2) then
                  Vmat(P,II,i,j)= VeR(i,j)
               else
                  STOP 'errVmat'
               end if           !end I if
            else
               STOP 'errVmat'
            end if              !end P if
 200     continue
 100  continue
      
         
C***********************************DEBUGGING***********************************C   
         
         if (debug) then
            open (unit=10,file='RotMatrix.out',status='new')
            
            do 400, p=1, 4
               do 500, II=1, 2
                  do 600, i=1, 3
                     do 700, j=1, 3
                        write(10,*) "Vmat(",P,II,i,j,")=",Vmat(P,II,i,j)            
 700                 continue
 600              continue
 500           continue
 400        continue

         
         
      end if                    !End debugguing if
      close(unit=10)

      return
      end                       !End subroutine RotMatrix
