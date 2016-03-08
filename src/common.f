      dimension q(0:9),i3(0:9),nc(0:9),v(0:9),a(0:9),v2(-2:12),
     .          a2(-2:12),eps2_L(-2:12),eps2_R(-2:12),
     .          asgrid(0:2000),an(4,0:2000),bn(4,0:2000)

      double precision alpha,gf,alfas0,asgrid,mu0,dahad3,pol4hf,an,bn,
     .                 pi1,pi2,zeta2,zeta3,zeta4,zeta5,zeta6,zeta7
      double precision q,i3,nc,v,a,alphat,sinhat,coshat,rhonc
      double precision mz,mh,mt,mb,mc,ms,md,mu,mtau,mmu,me,mp,mw
      double precision mw2,mz2,mt2,mh2
      double precision ratzw2,rathw2,rathz2,rattw2,rattz2,ratth2
      double precision SWpar,Spar,Tpar,Upar,Brho,Bkappa,Rpar,Zpar,
     .                 delrwh,delrzh,delr,rhohat,rho2,rho3,rho2a
      double precision mzp,mzp2,mz02,sinth,sinth2,costh,costh2,ratg21,
     .                 rhoeff,rhoezp,rhozzp,v2,a2,eps2_L,eps2_R,ratgRL
      double precision prob,mtp,sigma

      logical          flagmr,ffermi,fa2mt4,fa2mt2,fa2mt0,fla2im,
     .                 fasmt2,fasmt0,fas2mt,flagmf,fpolew,flgech,fobliq,
     .                 f4lqcd,falas2,fbayes,flagmh,flagmt,flagmc,flagS,
     .                 flagT,flgrho,fkappa,fzprim,fsinth,fwrite,flprob,
     .                 fhiggs,flagal,fsplot,flhdec

      integer          flgblm,flagzp,zprime

      common /inputs/  alpha,gf,alfas0,asgrid,mu0,dahad3,pol4hf,an,bn,
     .                 pi1,pi2,zeta2,zeta3,zeta4,zeta5,zeta6,zeta7
      common /coupls/  q,i3,nc,v,a,alphat,sinhat,coshat,rhonc
      common /masses/  mz,mh,mt,mb,mc,ms,md,mu,mtau,mmu,me,mp,mw
      common /mass2/   mw2,mz2,mt2,mh2
      common /ratios/  ratzw2,rathw2,rathz2,rattw2,rattz2,ratth2
      common /obliqe/  SWpar,Spar,Tpar,Upar,Brho,Bkappa,Rpar,Zpar,
     .                 delrwh,delrzh,delr,rhohat,rho2,rho3,rho2a
      common /zprime/  mzp,mzp2,mz02,sinth,sinth2,costh,costh2,ratg21,
     .                 rhoeff,rhoezp,rhozzp,v2,a2,eps2_L,eps2_R,ratgRL
      common /limits/  prob,mtp,sigma,zprime
      common /flags/   flagmr,flgblm,ffermi,fa2mt4,fa2mt2,fa2mt0,fla2im,
     .                 fasmt2,fasmt0,fas2mt,flagmf,fpolew,flgech,fobliq,
     .                 f4lqcd,falas2,fbayes,flagzp,flagmh,flagmt,flagmc,
     .                 flgrho,fkappa,flagS,flagT,fzprim,fsinth,fwrite,
     .                 flprob,fhiggs,flagal,fsplot,flhdec


