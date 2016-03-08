      subroutine zprimecoup(xval,npar)

      implicit none
      integer f,npar,imin
      double precision norm,alpzLR,xval(27) ! npar =27
      double precision aparam,bparam,cotalp,ratbsa,pm
      double precision cosprima2, cosprima,yu,ye,gbar,sinprima,gprima,xl
      double precision cosT, sinT, sinW, cosW, MZ1, MW1, sinW2, cte4, cte1
      double precision cte2, THETHA, tanW, cte3,cte6,cte7,c2w,equix
      include 'common.f'

C for the Z_chi defined through SO(10) -->  SU(5) X U(1)_chi:       flagzp = 1;
C for the Z_psi defined through    E_6 --> SO(10) X U(1)_psi:       flagzp = 2;
C for the Z_eta defined through sqrt(3/8) Z_chi - sqrt(5/8) Z_psi:  flagzp = 3;
C for the Z_LR defined through SU(2)_R X U(1) --> U(1)_Y X U(1)_LR: flagzp = 4;
C for the Z^prime defined to have the same couplings as the Z:      flagzp = 5;
C for the Z^string:                                                 flagzp = 6;
C for the general E_6 boson:                                        flagzp = 7;
C for the general family-universal Z^prime boson:                   flagzp = 8;
C for the model independent Z^prime boson:                          flagzp = 9;
C for the model independent Z^prime boson in v,a parametrization:   flagzp =10;
C for the general E_6 boson in alpha, beta parametrization:         flagzp =11;
C for the twisted E_6 boson:                                        flagzp =12;
C for the seesaw neutrino mass related one (Adhikari, Erler, Ma):   flagzp =13;
C for litlest Higgs                                                 flagzp =15;
C for the Z_N defined through 0.249*Z_chi + 0.969*Z_psi:            flagzp =14;
C for the Z_S  sequantial                                           flagzp =16;
c for  Z_{I}                                                        flagzp =17; 
c for  10+x5*                                                       flagzp =18;
c for  B-xL                                                         flagzp =19;
c for  q+xu                                                         flagzp =20;
c for  d-ux                                                         flagzp =21;
c for model A                                                       flagzp =22;   
c for model B                                                       flagzp =23;
c for model C                                                       flagzp =24;
c for model D                                                       flagzp =25;
c for model E                                                       flagzp =26;
c for model F                                                       flagzp =27;
c for model G                                                       flagzp =28;
c for model H                                                       flagzp =29;
c for model minimal                                                 flagzp =30;


c         if(universal.eqv..true.)then
c         imin=0
c         else
c         imin=-2
c         endif
          imin = 0

         if (flagzp.eq.1) then

          norm = 2*dsqrt(10.d0)
         
         eps2_L(0) =  3/norm
         eps2_L(1) =  3/norm
         eps2_L(2) =  3/norm
         eps2_L(3) =  3/norm
         eps2_L(4) = -1/norm
         eps2_L(5) = -1/norm
         eps2_L(6) = -1/norm
         eps2_L(7) = -1/norm
         eps2_L(8) = -1/norm
         eps2_L(9) = -1/norm
         
         eps2_R(0) =  5/norm
         eps2_R(1) =  1/norm
         eps2_R(2) =  1/norm
         eps2_R(3) =  1/norm
         eps2_R(4) =  1/norm
         eps2_R(5) =  1/norm
         eps2_R(6) =  1/norm
         eps2_R(7) = -3/norm
         eps2_R(8) = -3/norm
         eps2_R(9) = -3/norm

      else if (flagzp.eq.2) then
         
         eps2_L(0) =  1/dsqrt(24.d0)
         eps2_R(0) = -1/dsqrt(24.d0)
         do 10 f = 1, 9
            eps2_L(f) = eps2_L(0)
            eps2_R(f) = eps2_R(0)
 10      continue

      else if (flagzp.eq.3) then

         norm = 2*dsqrt(15.d0)
         
         eps2_L(0) =  1/norm
         eps2_L(1) =  1/norm
         eps2_L(2) =  1/norm
         eps2_L(3) =  1/norm
         eps2_L(4) = -2/norm
         eps2_L(5) = -2/norm
         eps2_L(6) = -2/norm
         eps2_L(7) = -2/norm
         eps2_L(8) = -2/norm
         eps2_L(9) = -2/norm
         
         eps2_R(0) =  5/norm
         eps2_R(1) =  2/norm
         eps2_R(2) =  2/norm
         eps2_R(3) =  2/norm
         eps2_R(4) =  2/norm
         eps2_R(5) =  2/norm
         eps2_R(6) =  2/norm
         eps2_R(7) = -1/norm
         eps2_R(8) = -1/norm
         eps2_R(9) = -1/norm

      else if (flagzp.eq.4) then

C         norm = 2*dsqrt(5/3.d0)
c         alpzLR = dsqrt(ratgRL**2*coshat/sinhat - 1)
c
c         eps2_L(0) =  1/norm/alpzLR
c         eps2_L(1) =  1/norm/alpzLR
c         eps2_L(2) =  1/norm/alpzLR
c         eps2_L(3) =  1/norm/alpzLR
c         eps2_L(4) = -1/norm/alpzLR/3
c         eps2_L(5) = -1/norm/alpzLR/3
c         eps2_L(6) = -1/norm/alpzLR/3
c         eps2_L(7) = -1/norm/alpzLR/3
c         eps2_L(8) = -1/norm/alpzLR/3
c         eps2_L(9) = -1/norm/alpzLR/3
c         
c         do 20 f = 0, 9
c            eps2_R(f) = eps2_L(f) + 2*i3(f)*alpzLR/norm
c 20      continue

         eps2_L(0) =   0.253996607
         eps2_L(1) =   0.253996607
         eps2_L(2) =   0.253996607 
         eps2_L(3) =   0.253996607
         eps2_L(4) =  -0.0846655357
         eps2_L(5) =  -0.0846655357
         eps2_L(6) =  -0.0846655357 
         eps2_L(7) =  -0.0846655357
         eps2_L(8) =  -0.0846655357
         eps2_L(9) =  -0.0846655357

         eps2_R(0) =   0
         eps2_R(1) =  -0.336562463
         eps2_R(2) =  -0.336562463 
         eps2_R(3) =  -0.336562463 
         eps2_R(4) =   0.505893534 
         eps2_R(5) =   0.505893534  
         eps2_R(6) =   0.505893534
         eps2_R(7) =  -0.675224605
         eps2_R(8) =  -0.675224605
         eps2_R(9) =  -0.675224605



      else if (flagzp.eq.5) then

C         do 30 f = 0, 9
C            eps2_L(f) =  i3(f) - q(f)*sinhat
C            eps2_R(f) =        - q(f)*sinhat
C 30      continue

c          eps2_L(0) = 0.5   
c          eps2_L(1) = -0.268735409
c          eps2_L(2) = -0.268735409
c          eps2_L(3) = -0.268735409
c          eps2_L(4) =  0.345823606
c          eps2_L(5) =  0.345823606
c          eps2_L(6) =  0.345823606
c          eps2_L(7) =  0.345823606
c          eps2_L(8) =  0.345823606 
c          eps2_L(9) =  0.345823606

c          eps2_r(0)=  0
c          eps2_r(1)=  0.231264591 
c          eps2_r(2)=  0.231264591 
c          eps2_r(3)=  0.231264591
c          eps2_r(4)= -0.154176394
c          eps2_r(5)= -0.154176394
c          eps2_r(6)= -0.154176394 
c          eps2_r(7)=  0.0770881972
c          eps2_r(8)=  0.0770881972
c          eps2_r(9)=  0.0770881972

         q(0) =  0.d0
         q(1) = -1.d0
         q(2) = -1.d0
         q(3) = -1.d0
         q(4) =  2.d0/3
         q(5) =  2.d0/3
         q(6) =  2.d0/3
         q(7) = -1.d0/3
         q(8) = -1.d0/3
         q(9) = -1.d0/3



         i3(0) =  1.d0/2
         i3(1) = -1.d0/2
         i3(2) = -1.d0/2
         i3(3) = -1.d0/2
         i3(4) =  1.d0/2
         i3(5) =  1.d0/2
         i3(6) =  1.d0/2
         i3(7) = -1.d0/2
         i3(8) = -1.d0/2
         i3(9) = -1.d0/2

         sinhat =   0.231285882d0
          
            do 37 f = 0, 9
            eps2_L(f) =  i3(f) - q(f)*sinhat
            eps2_R(f) =        - q(f)*sinhat
 37       continue





      else if (flagzp.eq.6) then

         norm = 100.d0
         
         eps2_L(0) = - 65/norm
         eps2_L(1) = -204/norm
         eps2_L(2) = - 65/norm
         eps2_L(3) =   74/norm
         eps2_L(4) =   68/norm
         eps2_L(5) =   68/norm
         eps2_L(6) =  -71/norm
         eps2_L(7) =   68/norm
         eps2_L(8) =   68/norm
         eps2_L(9) =  -71/norm
         
         eps2_R(0) =    0/norm
         eps2_R(1) =    9/norm
         eps2_R(2) =    9/norm
         eps2_R(3) = -130/norm
         eps2_R(4) = -  6/norm
         eps2_R(5) = -  6/norm
         eps2_R(6) =  133/norm
         eps2_R(7) =    3/norm
         eps2_R(8) =    3/norm
         eps2_R(9) = -136/norm

      else if (flagzp.eq.7) then

C   estos son los valores de xval para Z_R        
        xval(19) = 0d0
        xval(20) = -1d0
        xval(21) = 0d0


         norm = dsqrt(16*xval(19)*xval(20)/3 + 20*xval(19)*xval(21)/3
     .        +                                 4*xval(20)*xval(21)
     .        + 40*xval(19)**2/9 + 4*xval(20)**2 + 16*xval(21)**2)

         eps2_L(1) = (-   xval(19)                          )/norm
         eps2_L(2) = (-   xval(19)                          )/norm
         eps2_L(3) = (-   xval(19)                          )/norm
         eps2_L(4) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(5) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(6) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(7) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(8) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(9) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(0) = (-   xval(19)                          )/norm
         
         eps2_R(1) = (-   xval(19)/3 + xval(20) -   xval(21))/norm
         eps2_R(2) = (-   xval(19)/3 + xval(20) -   xval(21))/norm
         eps2_R(3) = (-   xval(19)/3 + xval(20) -   xval(21))/norm
         eps2_R(4) = (-   xval(19)/3 - xval(20) -   xval(21))/norm
         eps2_R(5) = (-   xval(19)/3 - xval(20) -   xval(21))/norm
         eps2_R(6) = (-   xval(19)/3 - xval(20) -   xval(21))/norm
         eps2_R(7) = (    xval(19)   + xval(20)             )/norm
         eps2_R(8) = (    xval(19)   + xval(20)             )/norm
         eps2_R(9) = (    xval(19)   + xval(20)             )/norm
         eps2_R(0) = (- 5*xval(19)/3 - xval(20) - 2*xval(21))/norm

C         eps2_L(3) = (+ 2*xval(19)/3 + xval(20) -   xval(21))/norm
C         eps2_R(9) = (- 2*xval(19)/3            +   xval(21))/norm

      else if ((flagzp.eq.11).or.(flagzp.eq.12)) then
c         xval(19)=aalpha
c         xval(20)=bbeta
         cotalp = 1d0/dtan(xval(19))
         ratbsa = dtan(xval(20))/dsin(xval(19))

c        pm = 1/+1  for  erler 1999 and  erler 2009 respectively
         
        pm = 1d0

         if (xval(19).eq.0.d0) then
            aparam = 0.d0
            bparam = -4*pm*tan(xval(20))/(tan(xval(20)) + 
     .              dsqrt(27/5.d0))/3
         else
            aparam = 10/(pm*dsqrt(10.d0)*ratbsa 
     .       + pm*dsqrt(54.d0)*cotalp - 6)
            bparam = -pm*dsqrt(8/45.d0)*aparam*ratbsa
         endif

         norm = dsqrt(16*aparam/3 + 20*bparam/3 + 4*aparam*bparam
     .        + 40/9.d0 + 4*aparam**2 + 16*bparam**2)

         eps2_L(1) = (- 1                         )/norm
         eps2_L(2) = (- 1                         )/norm
         eps2_L(3) = (- 1                         )/norm
         eps2_L(4) = (  1/3.d0          +   bparam)/norm
         eps2_L(5) = (  1/3.d0          +   bparam)/norm
         eps2_L(6) = (  1/3.d0          +   bparam)/norm
         eps2_L(7) = (  1/3.d0          +   bparam)/norm
         eps2_L(8) = (  1/3.d0          +   bparam)/norm
         eps2_L(9) = (  1/3.d0          +   bparam)/norm
         eps2_L(0) = (- 1                         )/norm
         
         eps2_R(1) = (- 1/3.d0 + aparam -   bparam)/norm
         eps2_R(2) = (- 1/3.d0 + aparam -   bparam)/norm
         eps2_R(3) = (- 1/3.d0 + aparam -   bparam)/norm
         eps2_R(4) = (- 1/3.d0 - aparam -   bparam)/norm
         eps2_R(5) = (- 1/3.d0 - aparam -   bparam)/norm
         eps2_R(6) = (- 1/3.d0 - aparam -   bparam)/norm
         eps2_R(7) = (  1      + aparam           )/norm
         eps2_R(8) = (  1      + aparam           )/norm
         eps2_R(9) = (  1      + aparam           )/norm
         eps2_R(0) = (- 5/3.d0 - aparam - 2*bparam)/norm

         if (flagzp.eq.12) then
            eps2_R(7) = (- 2/3.d0            + bparam)/norm
            eps2_R(8) = (- 2/3.d0            + bparam)/norm
            eps2_R(9) = (- 2/3.d0            + bparam)/norm
         endif

      else if (flagzp.eq.8) then

         eps2_L(1) = xval(19)
         eps2_L(2) = xval(19)
         eps2_L(3) = xval(19)
         eps2_L(4) = xval(21)
         eps2_L(5) = xval(21)
         eps2_L(6) = xval(21)
         eps2_L(7) = xval(21)
         eps2_L(8) = xval(21)
         eps2_L(9) = xval(21)
         eps2_L(0) = xval(19)
         
         eps2_R(1) = xval(20)
         eps2_R(2) = xval(20)
         eps2_R(3) = xval(20)
         eps2_R(4) = xval(22)
         eps2_R(5) = xval(22)
         eps2_R(6) = xval(22)
         eps2_R(7) = xval(23)
         eps2_R(8) = xval(23)
         eps2_R(9) = xval(23)
         eps2_R(0) = 0.d0

      else if (flagzp.eq.9) then

         eps2_L(1) = xval(19)
         eps2_L(2) = xval(19)
         eps2_L(3) = xval(24)
         eps2_L(4) = xval(21)
         eps2_L(5) = xval(21)
         eps2_L(6) = xval(26)
         eps2_L(7) = xval(21)
         eps2_L(8) = xval(21)
         eps2_L(9) = xval(26)
         eps2_L(0) = (2*xval(19) + xval(24))/3
         
         eps2_R(1) = xval(20)
         eps2_R(2) = xval(20)
         eps2_R(3) = xval(25)
         eps2_R(4) = xval(22)
         eps2_R(5) = xval(22)
         eps2_R(6) = 0.d0
         eps2_R(7) = xval(23)
         eps2_R(8) = xval(23)
         eps2_R(9) = xval(27)
         eps2_R(0) = 0.d0

      endif

      do 40 f = 0, 9 
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 40   continue

      if (flagzp.eq.10) then

         v2(1) = xval(19)
         v2(2) = xval(19)
         v2(3) = xval(24)
         v2(4) = xval(21)
         v2(5) = xval(21)
         v2(6) = 0.d0
         v2(7) = xval(22)
         v2(8) = xval(22)
         v2(9) = xval(26)
         v2(0) = (xval(19) + xval(20) + (xval(24) + xval(25))/2)/3

         a2(1) = xval(20)
         a2(2) = xval(20)
         a2(3) = xval(25)
         a2(4) = xval(23)
         a2(5) = xval(23)
         a2(6) = 0.d0
         a2(7) = xval(21) - xval(22) + xval(23)
         a2(8) = xval(21) - xval(22) + xval(23)
         a2(9) = xval(27)
         a2(0) = (xval(19) + xval(20) + (xval(24) + xval(25))/2)/3

         do 50 f = 0, 9 
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2
 50      continue

      endif

      if (flagzp.eq.13) then

         norm = 1.d0
         
         eps2_L(0) =  xval(19)/norm
         eps2_L(1) =  xval(19)/norm
         eps2_L(2) =  xval(19)/norm
         eps2_L(3) =  xval(19)/norm
         eps2_L(4) =  xval(20)/norm
         eps2_L(5) =  xval(20)/norm
         eps2_L(6) =  xval(20)/norm
         eps2_L(7) =  xval(20)/norm
         eps2_L(8) =  xval(20)/norm
         eps2_L(9) =  xval(20)/norm
         
C        eps2_R(0) =  (3*xval(19) + 9*xval(20))/4/norm  ! model (B)
         eps2_R(0) =  (3*xval(19) + 9*xval(20))/8/norm  ! model (C)
         eps2_R(1) =  (5*xval(19) - 9*xval(20))/4/norm
         eps2_R(2) =  (5*xval(19) - 9*xval(20))/4/norm
         eps2_R(3) =  (5*xval(19) - 9*xval(20))/4/norm
         eps2_R(4) = -(3*xval(19) - 7*xval(20))/4/norm
         eps2_R(5) = -(3*xval(19) - 7*xval(20))/4/norm
         eps2_R(6) = -(3*xval(19) - 7*xval(20))/4/norm
         eps2_R(7) =  (3*xval(19) +   xval(20))/4/norm
         eps2_R(8) =  (3*xval(19) +   xval(20))/4/norm
         eps2_R(9) =  (3*xval(19) +   xval(20))/4/norm

      do 60 f = 0, 9 
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 60   continue

         endif

         if (flagzp.eq.15) then

         norm = 0.5d0
         cosprima = 0.73
         cosprima2 = cosprima*cosprima
         sinprima = dsqrt(1-cosprima2)
         gprima = 0.346410162
         gbar = gprima/(2*sinprima*cosprima)
         xl = 0.2
         ye = 3/5
         yu = -2/5


        
         eps2_L(0) = 0
         eps2_L(1) =  gbar*(ye-1+cosprima2)/(norm) 
         eps2_L(2) =  gbar*(ye-1+cosprima2)/(norm)
         eps2_L(3) =  gbar*(ye-1+cosprima2)/(norm)
         eps2_L(4) =  gbar*(yu+2/3-(2*cosprima2)/3)/(norm)
         eps2_L(5) =  gbar*(yu+2/3-(2*cosprima2)/3)/(norm)
         eps2_L(6) =  gbar*(yu+7/5-(1*cosprima2)/3)/(norm)
         eps2_L(7) =  gbar*(yu+4/15+cosprima2/3)/(norm)
         eps2_L(8) =  gbar*(yu+4/15+cosprima2/3)/(norm)
         eps2_L(9) =  gbar*(yu+2/3-(2*cosprima2)/3-xl/5)/(norm)

        
         eps2_R(0) =  gbar*(ye-4/5+cosprima2/2)/(norm)  
         eps2_R(1) =  gbar*(ye-4/5+cosprima2/2)/(norm)
         eps2_R(2) =  gbar*(ye-4/5+cosprima2/2)/(norm)
         eps2_R(3) =  gbar*(ye-4/5+cosprima2/2)/(norm)
         eps2_R(4) =  gbar*(yu+7/15-(cosprima2)/6)/(norm)
         eps2_R(5) =  gbar*(yu+7/15-(cosprima2)/6)/(norm)
         eps2_R(6) =  gbar*(yu+7/15-(cosprima2)/3)/(norm)
         eps2_R(7) =  gbar*(yu+2/5-cosprima2/6)/(norm)
         eps2_R(8) =  gbar*(yu+2/5-cosprima2/6)/(norm)
         eps2_R(9) =  gbar*(yu+2/5-cosprima2/6)/(norm)
         do 70 f = 0, 9
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 70   continue

         endif


C introduzo ZN aqui
         if (flagzp.eq.14) then

         norm = dsqrt(40.d0)

         eps2_L(0) =  2/norm
         eps2_L(1) =  2/norm
         eps2_L(2) =  2/norm
         eps2_L(3) =  2/norm
         eps2_L(4) =  1/norm
         eps2_L(5) =  1/norm
         eps2_L(6) =  1/norm
         eps2_L(7) =  1/norm
         eps2_L(8) =  1/norm
         eps2_L(9) =  1/norm

         eps2_R(0) =  0.d0
         eps2_R(1) =  -1/norm
         eps2_R(2) =  -1/norm
         eps2_R(3) =  -1/norm
         eps2_R(4) =  -1/norm
         eps2_R(5) =  -1/norm
         eps2_R(6) =  -1/norm
         eps2_R(7) =  -2/norm
         eps2_R(8) =  -2/norm
         eps2_R(9) =  -2/norm
         
         endif

         if (flagzp.eq.16) then

         norm =  2*dsqrt(15.d0)

         eps2_L(0) =  4/norm
         eps2_L(1) =  4/norm
         eps2_L(2) =  4/norm
         eps2_L(3) =  4/norm
         eps2_L(4) = -0.5/norm
         eps2_L(5) = -0.5/norm
         eps2_L(6) = -0.5/norm
         eps2_L(7) = -0.5/norm
         eps2_L(8) = -0.5/norm
         eps2_L(9) = -0.5/norm


         eps2_R(0) = 5/norm
         eps2_R(1) = 0.5/norm
         eps2_R(2) = 0.5 /norm
         eps2_R(3) = 0.5/norm
         eps2_R(4) = 0.5/norm
         eps2_R(5) = 0.5/norm
         eps2_R(6) = 0.5/norm
         eps2_R(7) = -4/norm
         eps2_R(8) = -4/norm
         eps2_R(9) = -4/norm

       do 80 f = 0, 9
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 80    continue

        endif
        
         if (flagzp.eq.17) then

         norm =  2d0

         eps2_L(0) = -1/norm
         eps2_L(1) = -1/norm
         eps2_L(2) = -1/norm
         eps2_L(3) = -1/norm
         eps2_L(4) = 0/norm
         eps2_L(5) = 0/norm
         eps2_L(6) = 0/norm
         eps2_L(7) = 0/norm
         eps2_L(8) = 0/norm
         eps2_L(9) = 0/norm


         eps2_R(0) = -1/norm
         eps2_R(1) = 0/norm
         eps2_R(2) = 0/norm
         eps2_R(3) = 0/norm
         eps2_R(4) = 0/norm
         eps2_R(5) = 0/norm
         eps2_R(6) = 0/norm
         eps2_R(7) = 1/norm
         eps2_R(8) = 1/norm
         eps2_R(9) = 1/norm

       do 90 f = 0, 9
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 90     continue
         endif

         if (flagzp.eq.18) then ! 10+x5*

         norm =  2d0*dsqrt(4+equix+equix**2)/3d0

         eps2_L(0) = equix/(3*norm)
         eps2_L(1) = equix/(3*norm)
         eps2_L(2) = equix/(3*norm)
         eps2_L(3) = equix/(3*norm)
         eps2_L(4) = 1d0/(3d0*norm)
         eps2_L(5) = 1d0/(3d0*norm)
         eps2_L(6) = 1d0/(3d0*norm)
         eps2_L(7) = 1d0/(3d0*norm)
         eps2_L(8) = 1d0/(3d0*norm)
         eps2_L(9) = 1d0/(3d0*norm)
         eps2_L(10)= 0d0
         eps2_L(11)=-(1d0+equix)/(3*norm)
         eps2_L(12)=-2d0/(3d0*norm)

         eps2_R(0) =  (-2d0+equix)/(3d0*norm)
         eps2_R(1) =  -1d0/(3d0*norm)
         eps2_R(2) =  -1d0/(3d0*norm)
         eps2_R(3) =  -1d0/(3d0*norm)
         eps2_R(4) =  -1d0/(3d0*norm)
         eps2_R(5) =  -1d0/(3d0*norm)
         eps2_R(6) =  -1d0/(3d0*norm)
         eps2_R(7) =  -equix/(3d0*norm)
         eps2_R(8) =  -equix/(3d0*norm)
         eps2_R(9) =  -equix/(3d0*norm)
         eps2_R(10)= (-1d0-equix/3d0)/norm
         eps2_R(11)=  2d0/(3d0*norm)
         eps2_R(12)= (1d0+equix)/(3*norm)

         do 100 f = 0, 9
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 100     continue
         endif

          if (flagzp.eq.19) then 

          equix = 1d0  ! for B-xL

          norm =  -dsqrt(4d0*(1d0+equix**2)/(3d0))

         eps2_L(0) = -equix/norm
         eps2_L(1) = -equix/norm
         eps2_L(2) = -equix/norm
         eps2_L(3) = -equix/(norm)
         eps2_L(4) = 1d0/(3d0*norm)
         eps2_L(5) = 1d0/(3d0*norm)
         eps2_L(6) = 1d0/(3d0*norm)
         eps2_L(7) = 1d0/(3d0*norm)
         eps2_L(8) = 1d0/(3d0*norm)
         eps2_L(9) = 1d0/(3d0*norm)


         eps2_R(0) =  -1d0/norm
         eps2_R(1) = -equix/norm
         eps2_R(2) = -equix/norm
         eps2_R(3) = -equix/norm
         eps2_R(4) =  1d0/(3d0*norm)
         eps2_R(5) =  1d0/(3d0*norm)
         eps2_R(6) =  1d0/(3d0*norm)
         eps2_R(7) =  1d0/(3d0*norm)
         eps2_R(8) =  1d0/(3d0*norm)
         eps2_R(9) =  1d0/(3d0*norm)

         do 110 f = 0, 9
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 110     continue



         elseif (flagzp.eq.20) then ! q+xu

          norm =  (2d0/3d0)*dsqrt(7-2*equix+equix**2)

         eps2_L(0) = -1d0/norm
         eps2_L(1) = -1d0/norm
         eps2_L(2) = -1d0/norm
         eps2_L(3) = -1d0/norm
         eps2_L(4) = 1d0/(3d0*norm)
         eps2_L(5) = 1d0/(3d0*norm)
         eps2_L(6) = 1d0/(3d0*norm)
         eps2_L(7) = 1d0/(3d0*norm)
         eps2_L(8) = 1d0/(3d0*norm)
         eps2_L(9) = 1d0/(3d0*norm)
         eps2_L(10)= 0d0
         eps2_L(11)= 0d0
         eps2_L(12)= 0d0


         eps2_R(0) = (-4d0+equix)/(3d0*norm)
         eps2_R(1) = -(2d0+equix)/(3d0*norm)
         eps2_R(2) = -(2d0+equix)/(3d0*norm)
         eps2_R(3) = -(2d0+equix)/(3d0*norm)
         eps2_R(4) =  equix/(3d0*norm)
         eps2_R(5) =  equix/(3d0*norm)
         eps2_R(6) =  equix/(3d0*norm)
         eps2_R(7) =  (2d0-equix)/(3d0*norm)
         eps2_R(8) =  (2d0-equix)/(3d0*norm)
         eps2_R(9) =  (2d0-equix)/(3d0*norm)
         eps2_R(10)= 0d0
         eps2_R(11)= 0d0
         eps2_R(12)= 0d0

         do 120 f = 0, 9
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 120     continue


         elseif (flagzp.eq.21) then ! !  d-ux

         equix = 1d10 ! for down phobic

          norm =  (2d0/3d0)*dsqrt(1-equix+equix**2)

         eps2_L(0) = (-1d0+equix)/(3*norm)
         eps2_L(1) = (-1d0+equix)/(3*norm)
         eps2_L(2) = (-1d0+equix)/(3*norm)
         eps2_L(3) = (-1d0+equix)/(3*norm)
         eps2_L(4) =  0d0
         eps2_L(5) =  0d0
         eps2_L(6) =  0d0
         eps2_L(7) =  0d0
         eps2_L(8) =  0d0
         eps2_L(9) =  0d0
         eps2_L(10)= 0d0
         eps2_L(11)= -2*equix/(5*norm)
         eps2_L(12)= (1-4*equix/5d0)/(3d0*norm)  

         eps2_R(0) =  -equix/(3d0*norm)
         eps2_R(1) =   equix/(3d0*norm)
         eps2_R(2) =   equix/(3d0*norm)
         eps2_R(3) =   equix/(3d0*norm)
         eps2_R(4) =  -equix/(3d0*norm)
         eps2_R(5) =  -equix/(3d0*norm)
         eps2_R(6) =  -equix/(3d0*norm)
         eps2_R(7) =   1d0/(3d0*norm)
         eps2_R(8) =   1d0/(3d0*norm)
         eps2_R(9) =   1d0/(3d0*norm)
         eps2_R(10)=   0d0
         eps2_R(11)=  (-1d0+equix/5d0)/(3*norm)
         eps2_R(12)=  equix/(15d0*norm)

         do 130 f = 0, 9
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 130     continue

         elseif (flagzp.eq.0) then

          norm =  1d0

         eps2_L(0) =  0d0
         eps2_L(1) =  0d0
         eps2_L(2) =  0d0
         eps2_L(3) =  0d0
         eps2_L(4) =  0d0
         eps2_L(5) =  0d0
         eps2_L(6) =  0d0
         eps2_L(7) =  0d0
         eps2_L(8) =  0d0
         eps2_L(9) =  1d0


         eps2_R(0) =  0d0
         eps2_R(1) =  0d0
         eps2_R(2) =  0d0
         eps2_R(3) =  0d0
         eps2_R(4) =  0d0
         eps2_R(5) =  0d0
         eps2_R(6) =  0d0
         eps2_R(7) =  0d0
         eps2_R(8) =  0d0
         eps2_R(9) =  1d0

         do 140 f = 0, 9
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
140     continue
c model A
         elseif (flagzp.eq.22) then

         MW1 = 80.401d0
         MZ1 = 91.1876d0
         cosT = 1d0
         sinT = 0d0
         cosW = MW1/MZ1
         sinW = dsqrt(1d0-cosW*cosW)
         sinW2 = sinW*sinW
         cte4 = dsqrt(4d0*cosW*cosW-1d0)   
       
         v2(-2)= (0.5d0-sinW2)*cosT/cte4
         v2(-1)= (0.5d0-sinW2)*cosT/cte4 
         v2(0) = (0.5d0-sinW2)*cosT/cte4
         v2(1) = -(-0.5d0+2d0*sinW2)*cosT/cte4
         v2(2) = -(-0.5d0+2d0*sinW2)*cosT/cte4
         v2(3) = -(-0.5d0+2d0*sinW2)*cosT/cte4
         v2(4) = -(0.5d0-4d0*sinW2/3d0)*cosT/cte4
         v2(5) = -(0.5d0-4d0*sinW2/3d0)*cosT/cte4
         v2(6) = (0.5d0+sinW2/3d0)*cosT/cte4
         v2(7) = -(0.5d0-sinW2/3d0)*cosT/cte4
         v2(8) = -(0.5d0-sinW2/3d0)*cosT/cte4
         v2(9) = -(-0.5d0+2d0*sinW2/3d0)*cosT/cte4

         a2(-2)= (0.5d0-sinW2)*cosT/cte4
         a2(-1)= (0.5d0-sinW2)*cosT/cte4
         a2(0) = (0.5d0-sinW2)*cosT/cte4
         a2(1) = 0.5d0*cosT/cte4
         a2(2) = 0.5d0*cosT/cte4
         a2(3) = 0.5d0*cosT/cte4
         a2(4) = -0.5d0*cosT/cte4
         a2(5) = -0.5d0*cosT/cte4
         a2(6) = -(-0.5d0+sinW2)*cosT/cte4
         a2(7) = -(0.5d0-sinW2)*cosT/cte4
         a2(8) = -(0.5d0-sinW2)*cosT/cte4
         a2(9) = 0.5d0*cosT/cte4       
  
         do 150 f = imin, 9
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2
 150     continue

C model B
         elseif (flagzp.eq.23) then

         MW1 = 80.401d0
         MZ1 = 91.1876d0
         cosT = 1d0
         sinT = 0d0
         cosW = MW1/MZ1
         sinW = dsqrt(1d0-cosW*cosW)
         sinW2 = sinW*sinW
         cte4 = dsqrt(4d0*cosW*cosW-1d0)



         v2(0) = -0.5d0*cosT/cte4
         v2(1) = -(0.5d0+sinw2)*cosT/cte4
         v2(2) = -(0.5d0+sinw2)*cosT/cte4
         v2(3) = -(0.5d0+sinw2)*cosT/cte4
         v2(4) = (0.5d0+sinW2/3d0)*cosT/cte4
         v2(5) = (0.5d0+sinW2/3d0)*cosT/cte4
         v2(6) = (-0.5d0+4d0*sinW2/3d0)*cosT/cte4
         v2(7) = (0.5d0-2d0*sinW2/3d0)*cosT/cte4
         v2(8) = (0.5d0-2d0*sinW2/3d0)*cosT/cte4
         v2(9) = (-0.5d0+sinW2/3d0)*cosT/cte4


         a2(0) =  -0.5d0*cosT/cte4
         a2(1) =  -(0.5d0-sinw2)*cosT/cte4
         a2(2) =  -(0.5d0-sinw2)*cosT/cte4
         a2(3) =  -(0.5d0-sinw2)*cosT/cte4 
         a2(4) =  -(-0.5d0+sinw2)*cosT/cte4
         a2(5) =  -(-0.5d0+sinw2)*cosT/cte4
         a2(6) =  -0.5d0*cosT/cte4
         a2(7) =  0.5d0*cosT/cte4
         a2(8) =  0.5d0*cosT/cte4
         a2(9) =   -(0.5d0-sinw2)*cosT/cte4

         do 160 f = 0, 9
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2
         

160 	continue

c model C
         elseif (flagzp.eq.24) then

         MW1 = 80.401d0
         MZ1 = 91.1876d0
         cosT = 1d0
         sinT = 0d0
         cosW = MW1/MZ1
         sinW = dsqrt(1d0-cosW*cosW)
         sinW2 = sinW*sinW
         cte4 = dsqrt(4d0*cosW*cosW-1d0)


c con estas dos entradas iguales 
         v2(-2) = -(-0.5d0+sinW2)*cosT/cte4  !nu electronico 
         v2(-1) = -(-0.5d0+sinW2)*cosT/cte4  !nu muonico
	 v2(0) =  -0.5d0*cosT/cte4            !nu tau
         v2(1) =  -(0.5d0+sinW2)*cosT/cte4
         v2(2) =  -(-0.5d0+2d0*sinW2)*cosT/cte4
         v2(3) =  -3d0*(-0.5d0+sinW2)*cosT/cte4
         v2(4) =  -(0.5d0-4d0*sinW2/3d0)*cosT/cte4
         v2(5) =  -(0.5d0-4d0*sinW2/3d0)*cosT/cte4
         v2(6) =  -(-0.5d0-sinW2/3d0)*cosT/cte4
         v2(7) =  -(0.5d0-sinW2/3d0)*cosT/cte4
         v2(8) =  -(0.5d0-sinW2/3d0)*cosT/cte4
         v2(9) =  -(-0.5d0+2d0*sinW2/3d0)*cosT/cte4

         a2(-2) =   -(-0.5d0+sinW2)*cosT/cte4
         a2(-1) =   -(-0.5d0+sinW2)*cosT/cte4 
         a2(0) =    -0.5d0*cosT/cte4
         a2(1) =   -(0.5d0-sinW2)*cosT/cte4
         a2(2) =   0.5d0*cosT/cte4
         a2(3) =   (0.5d0-cosW*cosW)*cosT/cte4
         a2(4) =   -0.5d0*cosT/cte4
         a2(5) =   -0.5d0*cosT/cte4
         a2(6) =   -(-0.5d0+sinW2)*cosT/cte4
         a2(7) =   -(0.5d0-sinW2)*cosT/cte4
         a2(8) =   -(0.5d0-sinW2)*cosT/cte4
         a2(9) =   0.5d0*cosT/cte4
         imin = -2 
         do 170 f = imin, 9
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2


170     continue

C modelo D

	elseif (flagzp.eq.25) then

         MW1 = 80.401d0
         MZ1 = 91.1876d0
         cosT = 1d0
         sinT = 0d0
         cosW = MW1/MZ1
         sinW = dsqrt(1d0-cosW*cosW)
         sinW2 = sinW*sinW
         cte4 = dsqrt(4d0*cosW*cosW-1d0)


c tenemos dos posibilidades
         v2(-2) = -0.5d0*cosT/cte4
	 v2(-1) = -0.5d0*cosT/cte4
         v2(0) =  (0.5d0-sinW2)*cosT/cte4
         v2(1) =  -cosT*(0.5d0+sinW2)/cte4
         v2(2) =  -cosT*(0.5d0+sinW2)/cte4
         v2(3) =  -cosT*(-0.5d0+2d0*sinW2)/cte4
         v2(4) =  cosT*(0.5d0+sinW2/3d0)/cte4
         v2(5) =  cosT*(0.5d0+sinW2/3d0)/cte4
         v2(6) =  -cosT*(0.5d0-4d0*sinW2/3d0)/cte4
         v2(7) =  -cosT*(-0.5d0+2d0*sinW2/3d0)/cte4
         v2(8) =  -cosT*(-0.5d0+2d0*sinW2/3d0)/cte4
         v2(9) = -cosT*(0.5d0-sinW2/3d0)/cte4

         a2(-2) = -0.5d0*cosT/cte4
         a2(-1) = -0.5d0*cosT/cte4
         a2(0) =   (0.5d0-sinW2)*cosT/cte4 
         a2(1) =   -cosT*(0.5d0-sinW2)/cte4
         a2(2) =   -cosT*(0.5d0-sinW2)/cte4
         a2(3) =   0.5d0*cosT/cte4
         a2(4) =   -cosT*(-0.5d0+sinW2)/cte4
         a2(5) =   -cosT*(-0.5d0+sinW2)/cte4
         a2(6) =   -0.5d0*cosT/cte4
         a2(7) =   0.5d0*cosT/cte4
         a2(8) =   0.5d0*cosT/cte4
         a2(9) =   -cosT*(0.5d0-sinW2)/cte4
        imin = -2 
	do 180 f = imin, 9
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2


180     continue


c   Modelo E
	elseif (flagzp.eq.26) then
	

	 MW1 = 80.401d0
         MZ1 = 91.1876d0
         cosT = 1d0
         sinT = 0d0
         cosW = MW1/MZ1
         sinW = dsqrt(1d0-cosW*cosW)
         sinW2 = sinW*sinW
         cte4 = dsqrt(4d0*cosW*cosW-1d0)

c se tiene dos posibilidades
         v2(-2)= (0.5d0-sinW2)*cosT/cte4
         v2(-1)= (0.5d0-sinW2)*cosT/cte4
         v2(0) =  -0.5d0*cosT/cte4
         v2(1) =  cosT*(1.5d0-3d0*sinW2)/cte4  ! changed 26 marzo 5 -->5
         v2(2) =  cosT*(1.5d0-3d0*sinW2)/cte4  ! changed 26 marzo  "
         v2(3) = -cosT*(0.5d0+sinW2)/cte4      ! changed "  "        
         v2(4) =  cosT*(-0.5d0+4d0*sinW2/3d0)/cte4
         v2(5) =  cosT*(-0.5d0+4d0*sinW2/3d0)/cte4
         v2(6) =  cosT*(0.5d0+sinW2/3d0)/cte4
         v2(7) =  cosT*(-0.5d0+sinW2/3d0)/cte4
         v2(8) =  cosT*(-0.5d0+sinW2/3d0)/cte4
         v2(9) =  cosT*(0.5d0-2d0*sinW2/3d0)/cte4

        


         a2(-2)=   (0.5d0-sinW2)*cosT/cte4
         a2(-1)=   (0.5d0-sinW2)*cosT/cte4
         a2(0) =   -0.5d0*cosT/cte4
         a2(1) =   cosT*(0.5d0-cosW*cosW)/cte4
         a2(2) =   cosT*(0.5d0-cosW*cosW)/cte4
         a2(3) =   cosT*(-0.5d0+sinW2)/cte4
         a2(4) =   -0.5d0*cosT/cte4
         a2(5) =   -0.5d0*cosT/cte4
         a2(6) =   cosT*(0.5d0-sinW2)/cte4
         a2(7) =   cosT*(-0.5d0+sinW2)/cte4
         a2(8) =   cosT*(-0.5d0+sinW2)/cte4
         a2(9) =   0.5d0*cosT/cte4
	

         imin = -2
        do 190 f = imin, 9
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2


190     continue

c    Modelo F
	elseif (flagzp.eq.27) then


         MW1 = 80.401d0
         MZ1 = 91.1876d0
         cosT = 1d0
         sinT = 0d0
         cosW = MW1/MZ1
         sinW = dsqrt(1d0-cosW*cosW)
         sinW2 = sinW*sinW
         cte4 = dsqrt(4d0*cosW*cosW-1d0)

c se tiene dos posibilidades
         v2(-2) =  -0.5d0*cosT/cte4
         v2(-1) =  -0.5d0*cosT/cte4
         v2(0) =  cosT*(0.5d0-sinW2)/cte4
         v2(1) =  cosT*(-0.5d0-sinW2)/cte4
         v2(2) =  cosT*(-0.5d0-sinW2)/cte4
         v2(3) =  cosT*(1.5d0-3d0*sinW2)/cte4  ! changed 26 marzo 5->3
         v2(4) =  cosT*(0.5d0+sinW2/3d0)/cte4
         v2(5) =  cosT*(0.5d0+sinW2/3d0)/cte4
         v2(6) =  cosT*(-0.5d0+4d0*sinW2/3d0)/cte4
         v2(7) =  cosT*(0.5d0-2d0*sinW2/3d0)/cte4
         v2(8) =  cosT*(0.5d0-2d0*sinW2/3d0)/cte4
         v2(9) =  cosT*(-0.5d0+sinW2/3d0)/cte4

	 a2(-2) =  -0.5d0*cosT/cte4
	 a2(-1) =  -0.5d0*cosT/cte4
         a2(0) =   cosT*(0.5d0-sinW2)/cte4
         a2(1) =   cosT*(-0.5d0+sinW2)/cte4
         a2(2) =   cosT*(-0.5d0+sinW2)/cte4
         a2(3) =   cosT*(0.5d0-cosW*cosW)/cte4
         a2(4) =   cosT*(0.5d0-sinW2)/cte4
         a2(5) =   cosT*(0.5d0-sinW2)/cte4
         a2(6) =  -0.5d0*cosT/cte4
         a2(7) =   0.5d0*cosT/cte4
         a2(8) =   0.5d0*cosT/cte4
         a2(9) =   cosT*(-0.5d0+sinW2)/cte4
         imin = -2
        do 200 f = imin, 9
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2


200     continue

c model G
	elseif (flagzp.eq.28) then


         MW1 = 80.401d0
         MZ1 = 91.1876d0
         cosT = 1d0
         sinT = 0d0
         cosW = MW1/MZ1
         sinW = dsqrt(1d0-cosW*cosW)
         sinW2 = sinW*sinW
         cte4 = dsqrt(4d0*cosW*cosW-1d0)
         c2w = cosW*cosW-sinW*sinW

	 v2(0) =   0.5d0*(-sinT+cosT*(1d0-2*sinW2)/cte4)
         v2(1) =   -sinT*(-0.5d0+2d0*sinW2)+(0.5d0-sinW2)*3d0*cosT/cte4
         v2(2) =   -sinT*(-0.5d0+2d0*sinW2)+(0.5d0-sinW2)*3d0*cosT/cte4
         v2(3) =   -sinT*(-0.5d0+2d0*sinW2)+(0.5d0-sinW2)*3d0*cosT/cte4
         v2(4) =   (0.5d0 - 4d0*sinW2/3d0)*(-sinT-cosT/cte4)
         v2(5) =   (0.5d0 - 4d0*sinW2/3d0)*(-sinT-cosT/cte4)
         v2(6) =   (0.5d0 - 4d0*sinW2/3d0)*(-sinT-cosT/cte4)
         v2(7) = (-0.5d0 + 2d0*sinW2/3d0)*(-sinT)-(0.5d0-sinW2/3d0)*cosT/cte4
         v2(8) = (-0.5d0 + 2d0*sinW2/3d0)*(-sinT)-(0.5d0-sinW2/3d0)*cosT/cte4
         v2(9) = (-0.5d0 + 2d0*sinW2/3d0)*(-sinT)-(0.5d0-sinW2/3d0)*cosT/cte4


         a2(0) =   0.5d0*(-sinT+cosT*(1d0-2d0*sinW2)/cte4)
         a2(1) =   0.5d0*(sinT)+cosT*(0.5d0-cosW*cosW)/cte4
         a2(2) =   0.5d0*(sinT)+cosT*(0.5d0-cosW*cosW)/cte4
         a2(3) =   0.5d0*(sinT)+cosT*(0.5d0-cosW*cosW)/cte4
         a2(4) =   (-0.5d0*sinT-0.5d0*cosT/cte4)
         a2(5) =   (-0.5d0*sinT-0.5d0*cosT/cte4)
         a2(6) =   (-0.5d0*sinT-0.5d0*cosT/cte4)
         a2(7) =   -(cosT*0.5d0*c2w/cte4)
         a2(8) =   -(cosT*0.5d0*c2w/cte4)
         a2(9) =   -(cosT*0.5d0*c2w/cte4)

        do 210 f = 0, 9
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2


210     continue

c modelo H
	elseif (flagzp.eq.29) then


         MW1 = 80.401d0
         MZ1 = 91.1876d0
         cosT = 1d0
         sinT = 0d0
         cosW = MW1/MZ1
         sinW = dsqrt(1d0-cosW*cosW)
         sinW2 = sinW*sinW
         cte4 = dsqrt(4d0*cosW*cosW-1d0)

	 v2(0) = -0.5d0*cosT/cte4
         v2(1) = (-0.5d0-sinw2)*cosT/cte4
         v2(2) = (-0.5d0-sinw2)*cosT/cte4
         v2(3) = (-0.5d0-sinw2)*cosT/cte4
         v2(4) = (0.5d0+sinW2/3d0)*cosT/cte4
         v2(5) = (0.5d0+sinW2/3d0)*cosT/cte4
         v2(6) = (0.5d0+sinW2/3d0)*cosT/cte4
         v2(7) = (0.5d0-2d0*sinW2/3d0)*cosT/cte4
         v2(8) = (0.5d0-2d0*sinW2/3d0)*cosT/cte4
         v2(9) = (0.5d0-2d0*sinW2/3d0)*cosT/cte4


         a2(0) =  -0.5d0*cosT/cte4
         a2(1) =  (-0.5d0+sinw2)*cosT/cte4
         a2(2) =  (-0.5d0+sinw2)*cosT/cte4
         a2(3) =  (-0.5d0+sinw2)*cosT/cte4
         a2(4) =  (0.5d0-sinw2)*cosT/cte4
         a2(5) =  (0.5d0-sinw2)*cosT/cte4
         a2(6) =  (0.5d0-sinw2)*cosT/cte4
         a2(7) =  0.5d0*cosT/cte4
         a2(8) =  0.5d0*cosT/cte4
         a2(9) =  0.5d0*cosT/cte4

        do 220 f = 0, 9
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2

220     continue

c modelo minimal
        elseif (flagzp.eq.30) then

         MW1   = 80.401d0
         MZ1   = 91.1876d0
         cosT  = 1d0
         sinT  = dsqrt(1d0-cosT*cosT)
         cosW  = MW1/MZ1
         sinW  = dsqrt(1d0-cosW*cosW)
         sinW2 = sinW*sinW
         cte6  =  dsqrt(1d0-4d0*sinW2)
         cte7  = 2d0*dsqrt(3d0)  ! I remove the 2d0*
         
         
c     neutrino   
         v2(0) = cte6/cte7
c     e-
         v2(1) = 3d0*cte6/cte7
c     mu-
         v2(2) = 3d0*cte6/cte7
c     tau-       
         v2(3) = 3d0*cte6/cte7
c     up
         v2(4) = (-1d0+6d0*sinW2)/(cte6*cte7)
c     ch         
         v2(5) = (-1d0+6d0*sinW2)/(cte6*cte7)
c     top        
         v2(6) = (1d0+4d0*sinW2)/(cte6*cte7)
c     down
         v2(7) = (-1d0)/(cte6*cte7)
c     strange    
         v2(8) = (-1d0)/(cte6*cte7)
c     bottom     
         v2(9) = (1d0-2d0*sinW2)/(cte6*cte7)
                 
         a2(0) = cte6/cte7
         a2(1) = -cte6/cte7
         a2(2) = -cte6/cte7
         a2(3) = -cte6/cte7
         a2(4) = -(1d0+2d0*sinW2)/(cte6*cte7)
         a2(5) = -(1d0+2d0*sinW2)/(cte6*cte7)
         a2(6) = 2d0*(1d0-4d0*sinW2)/(cte6*cte7) ! pongo un factor 2 . t quark 
         a2(7) = (-1d0+4d0*sinW2)/(cte6*cte7)
         a2(8) = (-1d0+4d0*sinW2)/(cte6*cte7)
         a2(9) = (1d0+2d0*sinW2)/(cte6*cte7)
 
         do 230 f = 0, 9 
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2
 
 
230      continue



c         MW1 = 80.401d0
c         MZ1 = 91.1876d0
c         cosT = 1d0
c         sinT = 0d0
c         cosW = MW1/MZ1
c         sinW = dsqrt(1d0-cosW*cosW)
c         sinW2 = sinW*sinW
c         cte4 = dsqrt(4d0*cosW*cosW-1d0)

c         v2(0) = 0.288675134
c         v2(1) = -0.288675134
c         v2(2) = -0.288675134
c         v2(3) = -0.288675134
c         v2(4) = -0.288675134
c         v2(5) = 0.288675134
c         v2(6) = 0.288675134
c         v2(7) = -0.288675134
c         v2(8) = 0.288675134
c         v2(9) = 0.288675134


c         a2(0) =  0.288675134
c         a2(1) =  -0.288675134
c         a2(2) =  -0.288675134
c         a2(3) =  -0.288675134
c         a2(4) =  -0.288675134
c         a2(5) =  0.288675134
c         a2(6) =  0.288675134
c         a2(7) =  -0.288675134
c         a2(8) =  0.288675134
c         a2(9) =  0.288675134

c        do 230 f = 0, 9
c            eps2_L(f) = (v2(f) + a2(f))/2
c            eps2_R(f) = (v2(f) - a2(f))/2
c
c230	continue

         endif   

        return
        end

