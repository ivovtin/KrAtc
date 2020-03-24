      SUBROUTINE ATC_HIST
c************************************************************************
c                                                                       *
c  Calculation of amplitude distribution for integer Npe = 1,2,3,..     *
c                                                                       *
c  Used ROOT file for one-single photoelectron spectrum from MCP PMT    *
c                                                                       *
c  17.10.2002. V.A.Tayursky:                                            *
c  20.12.2018  I.V.Ovtin                                                *
c                                                                       *
c************************************************************************
c
c  SPEC(1024,0:I,cnt) - Array of amplitude distribution from MCP PMT
c                       for integer I photoelectrons.
c
      COMMON /ATCPMT/ RES, SPEC(1024,0:40,160)
      REAL MEAN
      integer*4 IC
c
      DATA SPEC/6717440*0./
c  Read out the single-photoelectron distribution from ROOT file
c-----------------------------------------------------------------
c  Transfotmation one-single photoelectron histogram
c  have been read out from ROOT file to SPEC array
c-----------------------------------------------------------------
c........
c   X axis scale parameters
c
c      RES  = ((XMA - XMI) / NX) / MEAN   ! resolution in ph.el./channel
c                                         ! (maximum 0.1 ph.el./channel)
      RES=0.04
      NPE_MAX = 40                        ! maximum number of photoelectrons
      NCHATC = NPE_MAX / RES              ! number of channels
c
c  Put one-single photoelectron histogram to SPEC array
c
      DO ICNT=1,160
c         print *,' ICNT=',ICNT
      DO I=1,1024
         CALL atc_read_1phe(%val(ICNT),%val(I),bins)      !read histogram with single-ph.e. spectrum
         SPEC(I,1,ICNT)=bins                              !fill channels for single-ph.e. spectrum
c         print *,' I = ',I,' BIN = ',bins,' SPEC(I,1,ICNT)=',SPEC(I,1,ICNT)
c
      END DO
c
c  Fill SPEC array              !насчитываем амплитудные спектры для целого числа ф.э. и заносим их в обший блок SPEC(1024,0:30)
c
      DO I=0,NPE_MAX-1
        IF(I.EQ.0) THEN
          SPEC(1,I,ICNT) = 1.
        ELSE
c
          CALL ATC_SPECN(NCHATC,SPEC(1,1,ICNT),SPEC(1,I,ICNT),SPEC(1,I+1,ICNT)) ! convolution
c          print *,' ICNT = ',ICNT,' I=',I,' SPEC(1,I+1,ICNT) = ',SPEC(1,I+1,ICNT)
c
        END IF
      END DO
      END DO
c
      RETURN
      END
