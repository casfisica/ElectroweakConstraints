      double complex Vmat(4,2,3,3),VPV(3,4,2,3,3), Del(3,3,4,2)
      double complex VuL(3,3),VuR(3,3),VdL(3,3),VdR(3,3),VnL(3,3)
      double complex VnR(3,3),VeL(3,3),VeR(3,3),epsi(4,2,3),gcop(4,2,3)
      double complex DB(3,3,4,2,3,3), DG(3,3,4,2,3,3),epsi1(4,2,3)
      double complex DC(4,4,4,3,3,3,3), DCM(4,4,4,3,3,3,3)
      double complex gcop1(4,2,3)

      common /deltalib/ Vmat,VPV,Del,VuL,VuR,VdL,VdR,VnL,        
     .                  VnR,VeL,VeR,epsi,gcop,DB,DG,epsi1,     
     .                  DC,DCM,gcop1
