      program main
      implicit none
      logical debug
      character(len=32) arg
      double complex Vmat(4,2,3,3)
      double complex VuL(3,3),VuR(3,3),VdL(3,3),VdR(3,3),VnL(3,3)
      double complex VnR(3,3),VeL(3,3),VeR(3,3)
      integer d

      
C*************************************Arguments**********************************C
      debug=.false.

      do d=1, iargc()
      call getarg(d, arg)
      end do

      if (arg.eq.'d') then
         debug=.true.
      end if
      
C*******************************************************************************C

      
      
      write(*,*) arg     
      call RotMatrix(Vmat,debug)
      
      
      stop
      end