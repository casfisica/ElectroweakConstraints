      program main
      implicit none
      logical debug, flinit
      character(len=32) arg
      double complex Vmat(4,2,3,3),VPV(3,4,2,3,3), Del(3,3,4,2)
      double complex VuL(3,3),VuR(3,3),VdL(3,3),VdR(3,3),VnL(3,3)
      double complex VnR(3,3),VeL(3,3),VeR(3,3),epsi(4,2,3)
      double complex DB(3,3,4,2,3,3), DG(3,3,4,2,3,3),epsi1(4,2,3)
      double complex DC(4,4,4,3,3,3,3), DCM(4,4,4,3,3,3,3)
      integer d, f, w,y

      include 'common.f'
      
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
C      call RotMatrix(Vmat,debug)
C      call Epsilon(epsi,gcop,debug)
C      call VPVt(VPV,debug)
C      call Delt(Del,debug)
C      call DBDG(DB,DG,debug)
C      call ZSM(epsi1,debug)
      call DCDCM(2,2,DC,DCM,debug)
      
      stop
      end
