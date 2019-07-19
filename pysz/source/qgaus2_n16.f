      SUBROUTINE qgaus2_n16(func,a,b,ss,q,r) ! one parameter
      REAL*8 a,b,ss,func,q,r
      EXTERNAL func
      INTEGER j
      REAL*8 dx,xm,xr,w(8),x(8)
      SAVE w,x
      DATA w/.1894506104550685d0,.1826034150449236d0,
     *.1691565193950025d0,.1495959888165767d0,
     *.1246289712555339d0,.0951585116824928d0,
     *.0622535239386479d0,.0271524594117541d0/
      DATA x/.0950125098376374d0,.2816035507792589d0,
     *.4580167776572274d0,.6178762444026438d0,
     *.7554044083550030d0,.8656312023878318d0,
     *.9445750230732326d0,.9894009349916499d0/
      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      ss=0d0
      do 11 j=1,8
        dx=xr*x(j)
        ss=ss+w(j)*(func(xm+dx,q,r)+func(xm-dx,q,r))
11    continue
      ss=xr*ss

      return
      END

C     (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
