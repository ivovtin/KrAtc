      SUBROUTINE ATC_SPECN(n,cont,cont1,cont2)
c*******************************************************************
c                                                                  *
c  This subroutine makes convolution of arrays cont and cont1.     *
c  cont2 - output array.                                           *
c                                                                  *
c*******************************************************************
c
      dimension cont(1024),cont1(1024),cont2(1024)
c
      rn = n
c
      s  = 0.
      s1 = 0.
c
      do i = 1,n
      c  = cont(i)
      s  = s + c
      c1 = cont1(i)
      s1 = s1 + c1
      end do
c
      do i = 1,n
      cont(i)  = cont(i)/s
      cont1(i) = cont1(i)/s1
      end do
c
      do i = 0,n-1
      c = 0
c
      do j = 0,i
      c = c + cont(j+1) * cont1(i-j+1)
      end do
c
      cont2(i+1) = c
c       print *,' c =',c
      end do
c
      RETURN
      END

