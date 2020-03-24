C Program to compute Pi using monte carlo methods
       program test_rand
       implicit none

       integer niter,i,j
       integer seed
       real*4 count
       real *8 x,y,pi(100),z
       real*8   rand

C initialize random numbers
       seed = 35791246
       call srand (seed)
       do j= 1,10
          niter = niter+1
          count =0
          do i=1,niter
c              call grndm(gamx,1)
c              call grndm(gamy,1)
c              x=gamx
c              y=gamy
c               x=0
c               y=0
             x=rand()
             y=rand()
             print*,'x=',x,' y=',y
             z= x*x +y*y
             if (z .le. 1) count =count+1
          end do
          pi(j)= count/niter*4.
          write(*,10) niter,pi(j)
c          print*,'test'
c10        format('Number of trials is:   estimate of pi is:',f8.5)
10        format('Number of trials is: 'I5,'  estimate of pi is:',f8.5)
       end do
       end
