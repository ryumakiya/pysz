      SUBROUTINE qgaus2_n14(func,a,b,ss,q,r) ! one parameter
      REAL*8 a,b,ss,func,q,r
      EXTERNAL func
      INTEGER j
      REAL*8 dx,xm,xr,w(7),x(7)
      SAVE w,x
      DATA w/.2152638534631578d0,.2051984637212956d0,
     *.1855383974779378d0,.1572031671581935d0,
     *.1215185706879032d0,.0801580871597602d0,
     *.0351194603317519d0/
      DATA x/.1080549487073437d0,.3191123689278897d0,
     *.5152486363581541d0,.6872929048116855d0,
     *.8272013150697650d0,.9284348836635735d0,
     *.9862838086968123d0/
      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      ss=0d0
      do 11 j=1,7
        dx=xr*x(j)
        ss=ss+w(j)*(func(xm+dx,q,r)+func(xm-dx,q,r))
11    continue
      ss=xr*ss

      return
      END

C     (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
