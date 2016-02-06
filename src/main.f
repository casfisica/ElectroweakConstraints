      program main
      implicit none
      logical debug, flinit
      character(len=32) arg
      integer d, f, w,y

      include 'common.f'
      include 'commondelta.f'
      
C*************************************Arguments**********************************C
      debug=.false.

      do d=1, iargc()
      call getarg(d, arg)
      end do

      if (arg.eq.'-d') then
         debug=.true.
      end if
      
C*******************************************************************************C
      flagzp=22
      flinit=.true.



      
C      call sin2thetaw(flinit)
C      write(*,*) arg     
      call RotMatrix(debug)
      call Epsilon(debug)
      call VPVt(debug)
      call Delt(debug)
      call DBDG(debug)
      call ZSM(debug)
      call DCDCM(21,200,debug)
      
      stop
      end
