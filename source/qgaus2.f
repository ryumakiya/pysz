      SUBROUTINE qgaus2(func,a,b,ss,q,r) ! one parameter
      REAL*8 a,b,ss,func,q,r
      EXTERNAL func
      INTEGER j
      REAL*8 dx,xm,xr,w(5),x(5)
      SAVE w,x
      DATA w/.2955242247d0,.2692667193d0,.2190863625d0,.1494513491d0,
     *.0666713443d0/
      DATA x/.1488743389d0,.4333953941d0,.6794095682d0,.8650633666d0,
     *.9739065285d0/
      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      ss=0d0
      do 11 j=1,5
        dx=xr*x(j)
        ss=ss+w(j)*(func(xm+dx,q,r)+func(xm-dx,q,r))
11    continue
      ss=xr*ss

      return
      END

C     (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
