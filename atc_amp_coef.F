      SUBROUTINE ATC_AMP_COEF(IC,AMPL,Iaer_CH,Itef_CH,Ish_CH,Iaer_SC,Itef_SC,Ish_SC,I_area,I_NLEVEL)
c -----------------------------------------------------------------
c                                                                 |
c  Determines average number of photo-electrons in Ch. counters   |
c  Inhomogeneity of light collection is taken into account.       |
c                                                                 |
c  Output parameters:                                             |
c                                                                 |
c    IC  - number of aerogel counter (1-160), particle in aerogel |
c    AMPL - Amplitude in number of photo-electrons                |
c                                                                 |
c   Last change 20.12.2017    I.V.Ovtin                           |
c------------------------------------------------------------------
c
      double precision Narea(68), Katc(68), Err_Kligt(68), Kcntf(160), cnt(160)
      REAL AMPL
      REAL Iaer_CH,Itef_CH,Ish_CH,Iaer_SC,Itef_SC,Ish_SC
      REAL Klight, Kcnt
c
      character fnameAS*128
      character fnameAQ*128
      character fnameAR*128
      character fnameAY*128
      character fnameEX*128
      character fnameKcnt*128
      character Dir*80
      character PrefixAS*20
      character PrefixAQ*20
      character PrefixAR*20
      character PrefixAY*20
      character PrefixEX*20
      character PrefixKcnt*20
      character ext*20
c
      AMPL=0.
      DNFOT=0.
      Aaer_CH=0.
      Ash_CH=0.
      Atef_CH=0.
      Aaer_SC=0.
      Ash_SC=0.
      Atef_SC=0.
c
      NLEVEL=I_NLEVEL
c      print *,' ****atc_amp_coef.F****** '
c      print *,'IC=',IC,' I_area = ',I_area,' NLEVEL = ',NLEVEL                       !print name of geometry volume
c      print *,'   Iaer_CH = ',Iaer_CH,' Ish_CH = ',Ish_CH,' Itef_CH = ',Itef_CH
c      print *,'   Iaer_SC = ',Iaer_SC,' Ish_SC = ',Ish_SC,' Itef_SC = ',Itef_SC
c
        DATA ICALL2/0/
c
        IF(ICALL2.GT.0) GOTO 1
	CALL get_atccalib_range(IREG)
c	print *,' IREG = ',IREG
c
	ICALL2=1
c
c        Dir='/home/ovtin/simulation/'
        Dir='/spool/users/atc/SimCalibCoef/'              !directory with calibration light of inhomogeneity
        PrefixAS='data_AS_'
        PrefixAQ='data_AQ_'
        PrefixAR='data_AR_'
        PrefixAY='data_AY_'
        PrefixEX='data_EX_'
        PrefixKcnt='data_Kcnt_all_'
        ext = 'range.dat'
c
        write(fnameAS,'(a,a,i2.2,a)') Dir(:lenocc(Dir)),
     &          PrefixAS(:lenocc(PrefixAS)),IREG,ext
        print *,'fnameAS:',fnameAS
c
        write(fnameAQ,'(a,a,i2.2,a)') Dir(:lenocc(Dir)),
     &          PrefixAQ(:lenocc(PrefixAQ)),IREG,ext
        print *,'fnameAQ:',fnameAQ
c
        write(fnameAR,'(a,a,i2.2,a)') Dir(:lenocc(Dir)),
     &          PrefixAR(:lenocc(PrefixAR)),IREG,ext
        print *,'fnameAR:',fnameAR
c
        write(fnameAY,'(a,a,i2.2,a)') Dir(:lenocc(Dir)),
     &          PrefixAY(:lenocc(PrefixAY)),IREG,ext
        print *,'fnameAY:',fnameAY
c
        write(fnameEX,'(a,a,i2.2,a)') Dir(:lenocc(Dir)),
     &          PrefixEX(:lenocc(PrefixEX)),IREG,ext
        print *,'fnameEX:',fnameEX
c
        write(fnameKcnt,'(a,a,i2.2,a)') Dir(:lenocc(Dir)),
     &          PrefixKcnt(:lenocc(PrefixKcnt)),IREG,ext
        print *,'fnameKcnt:',fnameKcnt
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  !calculate signal in Teflon - scintillation and cherenkov
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1	IF ( NLEVEL.EQ.3) THEN
c
      	 Atef_SC=0.25*Itef_SC
         Atef_CH=0.05*Itef_CH
c
c	 print *,'TEFLON NLEVEL = ',NLEVEL,' IC = ',IC

	ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c !calculate signal in Aerogel - scintillation and cherenkov
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ADEL=0
c
c  !define value
	Klight=4.0;
	Kcnt=1.5;
ccc
cc   !Take Klight for region of aerogel - long counters
ccc
        IF(Iaer_CH.GT.0) THEN
	DO i=1,68
	 	IF ( I_area.EQ.i .AND. ((40.LT.IC) .AND. (IC.LE. 60)) ) THEN
         		ADEL=1
         		open(98,file=fnameAS, status='old')
         		DO k=1,68
    	    		read(98,*) Narea(k), Katc(k), Err_Kligt(k)
               		IF(k.EQ.i) THEN
               			Klight=Katc(k)
               		ENDIF
			END DO
	    		close(98)
c	    		print *,Narea, Klight, Err_Kligt
c         		print *,'   i=',i,' Klight = ',Klight
         	ENDIF
c
	 	IF ( I_area.EQ.i .AND. ((100.LT.IC) .AND. (IC.LE. 120))) THEN
         		ADEL=1
         		open(98,file=fnameAR, status='old')
         		DO k=1,68
    	    		read(98,*) Narea(k), Katc(k), Err_Kligt(k)
               		IF(k.EQ.i) THEN
               			Klight=Katc(k)
	              	ENDIF
	    		END DO
	    		close(98)
c	    		print *,Narea, Klight, Err_Kligt
c         		print *,'   i=',i,' Klight = ',Klight
         	ENDIF
	END DO
ccc
cc   !Take Klight for region of aerogel - short counters
ccc
	DO i=1,58
	 	IF ( I_area.EQ.i .AND. ((20.LT.IC) .AND. (IC.LE. 40)) ) THEN
         		ADEL=1
         		open(98,file=fnameAY, status='old')
         		DO k=1,58
    	    		read(98,*) Narea(k), Katc(k), Err_Kligt(k)
               		IF(k.EQ.i) THEN
               			Klight=Katc(k)
               		ENDIF
	    		END DO
	    		close(98)
c	    		print *,Narea, Klight, Err_Kligt
c         		print *,'   i=',i,' Klight = ',Klight
         	ENDIF
c
	 	IF ( I_area.EQ.i .AND. ((120.LT.IC) .AND. (IC.LE. 140)) ) THEN
         		ADEL=1
         		open(98,file=fnameAQ, status='old')
         		DO k=1,58
    	    		read(98,*) Narea(k), Katc(k), Err_Kligt(k)
               		IF(k.EQ.i) THEN
               			Klight=Katc(k)
               		ENDIF
	    		END DO
	    		close(98)
c	    		print *,Narea, Klight, Err_Kligt
c         		print *,'   i=',i,' Klight = ',Klight
         	ENDIF
	END DO
ccc
cc   !Take Klight for region of aerogel - endcap counters
ccc
	DO i=1,28
	 	IF ( I_area.EQ.i .AND. (((1.LT.IC) .AND. (IC.LE. 20)) .OR. ((60.LT.IC) .AND. (IC.LE. 100)) .OR. ((140.LT.IC) .AND. (IC.LE. 160))) ) THEN
         		ADEL=1
ccc         		 open(98,file='/spool/users/atc/SimCalibCoef/data_EX_20092017_10range.dat', status='old')
         		open(98,file=fnameEX, status='old')
         		DO k=1,28
    	    		read(98,*) Narea(k), Katc(k), Err_Kligt(k)
               		IF(k.EQ.i) THEN
               			Klight=Katc(k)
               		ENDIF
	    		END DO
	    		close(98)
c	    		print *,Narea, Klight, Err_Kligt
c         		print *,'   i=',i,' Klight = ',Klight
         	ENDIF
	END DO
       	ENDIF
c
c Refrectiv coefficient for aerogel
c
         open(98,file=fnameKcnt, status='old')
         DO k=1,160
    	  read(98,*) cnt(k), Kcntf(k)
           IF(k.EQ.IC) THEN
              Kcnt=Kcntf(k)
c	      print *,cnt(k), Kcntf(k)
           ENDIF
	 END DO
	close(98)
c
      IF (NLEVEL.EQ.4 .OR. ADEL.EQ.1) THEN                 !eсли заряд не равен 0 и частица попала в аэрогель
      	   Aaer_SC=0.*Iaer_SC                              !Aaer_SC=0.00*DESTEP  Scintillation in aerogel is zero
  	   Aaer_CH=Klight*Kcnt*Iaer_CH                     !Aaer_CH=STEP*CHARGE*CHARGE*DNFOT
c      	   print *,'AEROGEL NLEVEL = ',NLEVEL,' IC = ',IC
      ENDIF
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c !calculate signal in Shifter - scintillation and cherenkov
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	Klight=20.00
c
        IF(Ish_CH.GT.0) THEN
	DO i=1,46
	 	IF ( I_area.EQ.i .AND. (((1.LT.IC) .AND. (IC.LE. 20)) .OR. ((60.LT.IC) .AND. (IC.LE. 100)) .OR. ((140.LT.IC) .AND. (IC.LE. 160))) ) THEN     !endcaps
 		open(98,file='/spool/users/atc/SimCalibCoef/data_SE.dat', status='old')
		DO k=1,46
		 read(98,*) Narea(k), Katc(k), Err_Kligt(k)
		 IF(k.EQ.i) THEN
		  Klight=Katc(k)
		 ENDIF
		END DO
		close(98)
c	         print *,Narea, Klight, Err_Kligt
c        	print *,'   i=',i,' Klight = ',Klight
        	ENDIF
	END DO
c
	DO i=1,66                                                                      !inner long barrel counter
	 	IF ( I_area.EQ.i .AND. ((40.LT.IC) .AND. (IC.LE. 60)) ) THEN
		open(98,file='/spool/users/atc/SimCalibCoef/data_SBIL.dat', status='old')
		DO k=1,66
		 read(98,*) Narea(k), Katc(k), Err_Kligt(k)
		 IF(k.EQ.i) THEN
		  Klight=Katc(k)
		 ENDIF
		END DO
		close(98)
c        	print *,'   i=',i,' Klight = ',Klight
        	ENDIF
	END DO
c
	DO i=1,58                                                                      !inner short barrel counter
	 	IF ( I_area.EQ.i .AND. ((20.LT.IC) .AND. (IC.LE. 40)) ) THEN
		open(98,file='/spool/users/atc/SimCalibCoef/data_SBIS.dat', status='old')
		DO k=1,58
		 read(98,*) Narea(k), Katc(k), Err_Kligt(k)
		 IF(k.EQ.i) THEN
		  Klight=Katc(k)
		 ENDIF
		END DO
		close(98)
c        	print *,'   i=',i,' Klight = ',Klight
        	ENDIF
	END DO
c
	DO i=1,66                                                                      !outer long barrel counterк
	 	IF ( I_area.EQ.i .AND. ((100.LT.IC) .AND. (IC.LE. 120)) ) THEN
		open(98,file='/spool/users/atc/SimCalibCoef/data_SBOL.dat', status='old')
		DO k=1,66
		 read(98,*) Narea(k), Katc(k), Err_Kligt(k)
		 IF(k.EQ.i) THEN
		  Klight=Katc(k)
		 ENDIF
		END DO
		close(98)
c        	print *,'   i=',i,' Klight = ',Klight
        	ENDIF
	END DO
c
	DO i=1,58                                                                      !outer short barrel counter
	 	IF ( I_area.EQ.i .AND. ((120.LT.IC) .AND. (IC.LE. 140)) ) THEN
		open(98,file='/spool/users/atc/SimCalibCoef/data_SBOS.dat', status='old')
		DO k=1,58
		 read(98,*) Narea(k), Katc(k), Err_Kligt(k)
		 IF(k.EQ.i) THEN
		  Klight=Katc(k)
		 ENDIF
		END DO
		close(98)
c        	print *,'   i=',i,' Klight = ',Klight
        	ENDIF
	END DO
c
	ENDIF
c
      IF ( NLEVEL.EQ.6 .OR. (NLEVEL.EQ.5 .AND. ADEL.NE.1)) THEN
c      	  print *,'SHIFTER NLEVEL = ',NLEVEL,' IC = ',IC
          Ash_SC=0.005*DESTEP
          Ash_CH=Klight*Ish_CH                                                !Ash_CH=STEP*CHARGE*CHARGE*DNFOTsh
      ENDIF
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c !Result:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  	AMPL=(Aaer_CH+Atef_CH+Ash_CH+Aaer_SC+Atef_SC+Ash_SC)  			       !average number of ph.e.
c       print *,'   IC=',IC,' NLEVEL = ',NLEVEL,' I_area = ',I_area,' ADEL=',ADEL
c      	print *,'   Aaer_CH = ',Aaer_CH,' Ash_CH = ',Ash_CH,' Atef_CH = ',Atef_CH
c      	print *,'   Aaer_SC = ',Aaer_SC,' Ash_SC = ',Ash_SC,' Atef_SC = ',Atef_SC
c       print *,'   AMPL = ',AMPL
c
      END
