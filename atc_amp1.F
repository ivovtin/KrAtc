#include "geant321/pilot.h"
#include "KcSys/fortran.h"
c
        SUBROUTINE ATC_AMP1(IC,NPE,AMP)
c
c*******************************************************************
c                                                                  *
c   This subroutine generates the amplitude from ATC counter       *
c   Other subroutines in this module: ATC_HIST, ATC_SPECN          *
c                                                                  *
c   Input parameter: NPE - avarage number of photoelectrons        *
c                                                                  *
c   Output parameter: AMP - amplitude from ATC counter (in p-e)    *
c                                                                  *
c   Author: S.F.Ganzhur, I.V.Ovtin                                 *
c                                                                  *
c*******************************************************************
c  SPEC(500,I) - ampl. distribution from MCP PMT for integer I pe
c
      COMMON /ATCPMT/ RES, SPEC(1024,0:40,160)
c
c  PMT(500) - array of amplitude distribution from MCP PMT for Npe.
c
      DIMENSION PMT(1024), POI(0:40)
      integer*4 IC
      REAL NPE
      REAL AMP
c........
c-> X axis scale parameters:
c     RES     - Resolution in phe/chennel (maximum 0.1 phe/ch)
c     NPE_MAX - maximum number of photoelectron
c     NCHATC     - number of channels
c
      DATA ICALL/0/
      integer*4 ev
c
c      print *,' IC=',IC,' NPE=',NPE
      RES=0.04
c
      IF(ICALL.GT.0) GOTO 1
      CALL ATC_HIST                                     ! preparatiion of arrays for ampl.generation
c
c**********************
      ICALL=1
c
1     NPE_MAX = 30
      NCHATC = NPE_MAX / RES                            !number of channels×
c.......
c
c  POI(I) - Poisson distribution for NPE                ! P(n)=((mu)^n/n!)*e^(-mu)
c
      DO J = 0,NPE_MAX
      IF (J .EQ. 0) THEN
      POI(0)  = EXP(-NPE)                               !probably of zero
      ELSE
      POI(J)=(NPE**J)*EXP(-NPE)
c
      DO I=1,J
      POI(J)=POI(J)/I
      END DO
c
      ENDIF
c      CALL atc_poisson(%val(POI(J)),%val(J))            !test - draw Poisson distribution
      END DO
c      print *,'IC=',IC,' NPE=',NPE
c      CALL atc_poisson_fill(%val(IC),%val(NPE))
c----------------------------------------------------
c---> PMT(I) calculation
c
      DO I =1,1024
      PMT(I) = 0.
      END DO
c
      DO I = 0,NPE_MAX
c
      DO J=1,NCHATC
      PMT(J)=PMT(J)+POI(I)*SPEC(J,I,IC)                 !convolution Poisson distribution with single ph.e. spectrum
c
c      CALL atc_pmt_signal(%val(PMT(J)),%val(J))
c      print *,'I=',I,' J=',J,' PMT(J)=',PMT(J),' SPEC(J,I,IC)=',SPEC(J,I,IC)
      END DO
c
      END DO
c      CALL atc_pmt_signal_fill(%val(IC),%val(NPE))
c
c----------------------------------------------------------
c---> Transformation average number of photoelectrons from
c     the ATC counter to the real one taking into account
c     electronics and statistic distribution.
c----------------------------------------------------------
c-> Y axis scale parameter:
c
      RMAX = PMT(1)
c
      DO I = 1,NCHATC-1
      IF(RMAX .LE. PMT(I+1)) RMAX = PMT(I+1)
      END DO
ccc........
c-> X generation
999   GAMX=rand()                                       !random number generator in the interval (0,1)
c      print *,' GAMX=',GAMX
      GAMX = GAMX * NPE_MAX                             !get random number ph.e. from 0 to NPE_MAX
      I_HI = (GAMX / RES) + 1                           !get number of channel
c
c      print *,' RES=',RES
c
      FUN  = PMT(I_HI)                                  !probably in this channel
c      print *,' GAMX=',GAMX,' RES=',RES,' I_HI=',I_HI,' FUN=',FUN
c........
c-> Y generation
c
      GAMY=rand()
c      print *,' GAMY=',GAMY
c
      GAMY = GAMY * RMAX
c........
c---> Checking of condition
C
      IF(FUN.LT.GAMY) GOTO 999
      AMP = GAMX                                        !Amplitude in number ph.e. - random from spectrum PMT(J)
c      print *,' AMP=',AMP
      RETURN
c
      END
