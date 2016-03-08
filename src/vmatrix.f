      subroutine Vmatrix(VuL,VuR,VdL,VdR,VnL,VnR,VeL,VeR)     
      implicit none
      double complex VuL(3,3),VuR(3,3),VdL(3,3),VdR(3,3),
     . VnL(3,3),VnR(3,3),VeL(3,3),VeR(3,3)
      double precision mc,mt,du,mu,md,ms,mb
      integer i,j

      do 100 i=1,3
      do 101 j=1,3
      VuL(3,3)=0d0
      VuR(3,3)=0d0
      VdL(3,3)=0d0
      VdR(3,3)=0d0
      VnL(3,3)=0d0
      VnR(3,3)=0d0
      VeL(3,3)=0d0
      VeR(3,3)=0d0
101   continue
100   continue
 



        mt = 172d0 
        mc = 0.560d0
        mu = 0.0024d0 
        mb = 2.89d0 
        ms = 0.06d0 
        md = 0.0029d0
        du = 171.721d0

      
       VuL(1,1)= (dsqrt(mc)*dsqrt(mt)*dsqrt(du - mu))/
     .           (dsqrt(du)*dsqrt(mt - mu)*dsqrt(mc + mu))
       VuL(1,2)= (dsqrt(du - mu)*dsqrt(mu))/
     .           (dsqrt(mt - mu)*dsqrt(mc + mu))
       VuL(1,3)=- dsqrt(du + mc)*dsqrt(-du + mt)*dsqrt(mu)/
     .           (dsqrt(du)*dsqrt(mt - mu)*dsqrt(mc + mu)) 
       VuL(2,1)= dsqrt(du + mc)*dsqrt(mt)*dsqrt(mu)/
     .           (dsqrt(du)*dsqrt(mc + mt)*dsqrt(mc + mu))
       VuL(2,2)=- dsqrt(mc)*dsqrt(du + mc)/
     .           (dsqrt(mc + mt)*dsqrt(mc + mu))
       VuL(2,3)= (dsqrt(mc)*dsqrt(-du + mt)*dsqrt(du - mu))/
     .           (dsqrt(du)*dsqrt(mc + mt)*dsqrt(mc + mu))
       VuL(3,1)= (dsqrt(mc)*dsqrt(-du + mt)*dsqrt(mu))/
     .           (dsqrt(du)*dsqrt(mc + mt)*dsqrt(mt - mu))
       VuL(3,2)= (dsqrt(mt)*dsqrt(-du + mt))/
     .           (dsqrt(mc + mt)*dsqrt(mt - mu))
       VuL(3,3)= dsqrt(du + mc)*dsqrt(mt)*dsqrt(du - mu)/
     .           (dsqrt(du)*dsqrt(mc + mt)*dsqrt(mt - mu))  



      VuR(1,1)=1d0
      VuR(2,2)=1d0
      VuR(3,3)=1d0
  


     
      VdL(1,1)=  dSqrt(ms/(md + ms))
      VdL(1,2)=  dSqrt(md/(md + ms))
      VdL(1,3)=  0d0
      VdL(2,1)= -dSqrt(md/(md + ms))
      VdL(2,2)=  dSqrt(ms/(md + ms))
      VdL(2,3)=  0d0
      VdL(3,1)=  0d0
      VdL(3,2)=  0d0
      VdL(3,3)=  1d0
   
      VdR(1,1)=1d0
      VdR(2,2)=1d0
      VdR(3,3)=1d0




      VnL(1,1)=1d0
      VnL(2,2)=1d0
      VnL(3,3)=1d0

      VnR(1,1)=1d0
      VnR(2,2)=1d0
      VnR(3,3)=1d0


      VeL(1,1)=1d0
      VeL(2,2)=1d0
      VeL(3,3)=1d0

      VeR(1,1)=1d0
      VeR(2,2)=1d0
      VeR(3,3)=1d0



      

      
       end
