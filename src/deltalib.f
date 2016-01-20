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
      double complex V(4,2,3,3) 
      integer P,I,i,j

      include 'common.f'
C      include 'commondelta'

C     Import the values of the rotation matrix
      call Vmatrix(VuL,VuR,VdL,VdR,VnL,VnR,VeL,VeR)

C     Organize the values of rotation matrix as an estructure Vmat(P,I,i,j)
C     where P: Fermion type; 1:UP, 2:DOWN, 3:Neutrino, 4: Electron
C     I: Chirality; 1:Left, 2:Right
C     i,j: Flavor or Families; 1:electron, 2:muon, 3:Tau
      
      do 100, P=1, 4            !Fermion type
         do 200, I=1, 2         !Chirality
            do 300, i=1, 3      !Family
               do 400, j=1, 3   !Family

                  if (P .eq. 1) then

                     if (I .eq. 1) then

                        Vmat(P,I,i,j)= VuL(i,j)
                        
                     else if (I .eq. 2) then

                        Vmat(P,I,i,j)= VuR(i,j)
                        
                     else

                        STOP 'errVmat'
                        
                     end if     !end I if

                     
                  else if (p .eq. 2) then

                     if (I .eq. 1) then

                        Vmat(P,I,i,j)= VdL(i,j)
                        
                     else if (I .eq. 2) then

                        Vmat(P,I,i,j)= VdR(i,j)
                        
                     else

                        STOP 'errVmat'
                        
                     end if     !end I if
                     
                  else if (P .eq. 3) then
                     
                     if (I .eq. 1) then

                        Vmat(P,I,i,j)= VnL(i,j)
                        
                     else if (I .eq. 2) then

                        Vmat(P,I,i,j)= VnR(i,j)
                        
                     else

                        STOP 'errVmat'
                        
                     end if     !end I if
                     
                  else if (P .eq. 4) then

                     if (I .eq. 1) then

                        Vmat(P,I,i,j)= VeL(i,j)
                        
                     else if (I .eq. 2) then

                        Vmat(P,I,i,j)= VeR(i,j)
                        
                     else

                        STOP 'errVmat'
                        
                     end if     !end I if
                     
                  else

                     STOP 'errVmat'
                     
                  end if        !end P if
                     

               
 300        continue
            
 200     continue
         
 100  continue
