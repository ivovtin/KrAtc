c
c    Array dimension constants
c
        integer ATC_NTRK, ATC_NCROS
        parameter (ATC_NTRK=16, ATC_NCROS=16)
        integer ATC_NCNT
        parameter (ATC_NCNT=160)
c
c    Global ATC data
c
c      ATC geometry (atc_geom common-block)
        real gatc_r(ATC_NCNT)        ! rho-position of counter centers
	real gatc_phi(ATC_NCNT)      ! phi-position of counter centers
	real gatc_z(ATC_NCNT)        ! z-position of counter centers
        real gatc_halfr(ATC_NCNT)    ! half rho-size of counters
	real gatc_halfphi(ATC_NCNT)  ! half phi-size of counters
	real gatc_halfz(ATC_NCNT)    ! half z-size of counters
c      ATC calibrations (atc_data common-block)
        integer gatc_bad(ATC_NCNT)   ! bad counter signatures
        real gatc_amp_err(ATC_NCNT)  ! amplitude measurement errors
        real gatc_a0(ATC_NCNT)       ! array of channel pedestals
        real gatc_tau(ATC_NCNT)      ! array of tau parameters
        real gatc_p(ATC_NCNT)        ! array of p parameters
        real gatc_a1pe(ATC_NCNT)     ! single photoelectron amplitudes
        real gatc_da1pe(ATC_NCNT)    ! single photoelectron amplitude errors
	real gatc_minDeltaT          ! DeltaT range lower limit, counts
	real gatc_maxDeltaT          ! DeltaT range upper limit, counts
	real gatc_dcSigmaT           ! transversal coordinate resolution at ATC radius, mm
	real gatc_dcSigmaZ           ! longitudinal coordinate resolution at ATC radius, mm
c
c    Binding DC tracks to ATC counters (atc_track common-block)
c
        integer   atc_ntracks                       ! number of tracks (same as eTracksAll)
        integer*1 atc_trk_ncnt(ATC_NTRK)            ! number of counters crossed by the track      !INTEGER*1            [-128  +127]
        integer atc_trk_cnts(ATC_NCROS,ATC_NTRK)  ! list of counter indexes on the track 1..160       !изменил тип данных - integer*1
        integer*1 atc_illtrack(ATC_NTRK)            ! whether DC track is wrongly reconstructed as determined by ATC
        integer*1 atc_wls_hit(ATC_NCROS,ATC_NTRK)   ! whether WLS is crossed by the track in the counter
        integer*1 atc_near_wls(ATC_NCROS,ATC_NTRK)  ! whether the track had possibility to cross WLS
        integer*1 atc_track_end(ATC_NCROS,ATC_NTRK) ! which end of the track crossed the counter: 0 or 1 (for cosmic run only)
        real atc_trk_aer_len(ATC_NCROS,ATC_NTRK)    ! track length in aerogel for the counter, mm
        real atc_trk_wls_len(ATC_NCROS,ATC_NTRK)    ! track length in WLS for the counter, mm
c      local coordinates of crossing counters
        real atc_rin(ATC_NCROS,ATC_NTRK)            ! local rho coordinate of inner crossing of the counter, mm
        real atc_rout(ATC_NCROS,ATC_NTRK)           ! local rho coordinate of outer crossing of the counter, mm
        real atc_phiin(ATC_NCROS,ATC_NTRK)          ! local phi coordinate of inner crossing of the counter, rad
        real atc_phiout(ATC_NCROS,ATC_NTRK)         ! local phi coordinate of outer crossing of the counter, rad
        real atc_zin(ATC_NCROS,ATC_NTRK)            ! local z coordinate of inner crossing of the counter, mm
        real atc_zout(ATC_NCROS,ATC_NTRK)           ! local z coordinate of outer crossing of the counter, mm
        real atc_ph_in(ATC_NCROS,ATC_NTRK)          ! track phase for inner crossings of the counter, rad
        real atc_ph_out(ATC_NCROS,ATC_NTRK)         ! track phase for outer crossings of the counter, rad
c      global coordinates of crossing counter central plane (WLS)
        real atc_rwls(ATC_NCROS,ATC_NTRK)           ! rho coordinate of crossing WLS, mm
	real atc_phiwls(ATC_NCROS,ATC_NTRK)         ! phi coordinate of crossing WLS, rad
        real atc_zwls(ATC_NCROS,ATC_NTRK)           ! z coordinate of crossing WLS, mm
c      global coordinates of crossing barrel&endcap layers
        real atc_rB(8,ATC_NTRK)                     ! rho coordinate of crossing barrel surfaces, mm
	real atc_phiB(8,ATC_NTRK)                   ! phi coordinate of crossing barrel surfaces, rad
	real atc_zB(8,ATC_NTRK)                     ! z coordinate of crossing barrel surfaces, mm
        real atc_rE(8,ATC_NTRK)                     ! rho coordinate of crossing endcap surfaces, mm
	real atc_phiE(8,ATC_NTRK)                   ! phi coordinate of crossing endcap surfaces, rad
	real atc_zE(8,ATC_NTRK)                     ! z coordinate of crossing endcap surfaces, mm
c      bind counters to tracks
        integer*1 atc_cnt_ntrk(ATC_NCNT)            ! number of tracks in counter
	integer*1 atc_cnt_tracks(ATC_NTRK,ATC_NCNT) ! track indexes (from 1) bound to counter
c
c    Hits in ATC counters (atc_rec common-block)
c    Indexed by physical counter number, not hit number.
c
        real atc_npe(ATC_NCNT), atc_dnpe(ATC_NCNT)  ! amplitudes in number of photoelectrons
        real atc_time(ATC_NCNT)                     ! times of hits (DeltaT corrected)
        real atc_amp(ATC_NCNT)                      ! amplitudes in ADC counts
        real atc_rawtime(ATC_NCNT)                  ! raw times of the hit (not corrected)
        real atc_chi2(ATC_NCNT)                     ! Chi-squares of pulse fit
        integer*1 atc_trg(ATC_NCNT)                 ! whether the counter is present in the event
        integer   atc_htyp(ATC_NCNT)                ! hit type determined by AtcHit
        integer   atc_eventdamage                   ! event damage type, 0 if no damage
        integer   atc_rawDeltaT                     ! DeltaT in counts 0..255
        real atc_deltaT                             ! DeltaT correction in ticks -0.5..0.5
c
c       ATC information to track
c
	integer atcinfo_t;                                                    !number of tracks (same as eTracksAll)
	integer atcinfo_ncnt;                                                 !number of counters crossed by the track
	integer atcinfo_cnt(ATC_NCROS);                                       !list of counter numbers on the track 1..160
	real atcinfo_npe(ATC_NCNT),atcinfo_npen(ATC_NCNT);         	      !the hit amplitude in number of photoelectrons for each counte
	integer atcinfo_wlshit(ATC_NCNT);                                     !WLS crossed
	integer atcinfo_nearwls(ATC_NCNT);				      !track goes close to WLS
	real atcinfo_tlen(ATC_NCNT), atcinfo_pathwls(ATC_NCNT);	              !track length in aerogel volume and that in WLS
c       пересечение треком различных областей счетчика (0 или 1)
	integer aerogel_REGION(ATC_NCNT);
	integer aerogel_REGION0(ATC_NCNT);                                    !ATC_DOUBLE_CROSS
	integer aerogel_REGION5(ATC_NCNT);
	integer aerogel_REGION20(ATC_NCNT);
	integer active_REGION(ATC_NCNT);
	integer active_REGION0(ATC_NCNT);
	integer active_REGION5(ATC_NCNT);
	integer active_REGION20(ATC_NCNT);
	integer testreg(ATC_NCNT);				              !cosmic test region of counter crossed at total thickness
	integer single_aerogel_REGION(ATC_NCNT);
	integer single_aerogel_REGION0(ATC_NCNT);           		      !ATC_SINGLE_CROSS
	integer single_aerogel_REGION5(ATC_NCNT);
	integer single_aerogel_REGION20(ATC_NCNT);
	integer single_active_REGION(ATC_NCNT);
	integer single_active_REGION0(ATC_NCNT);
	integer single_active_REGION5(ATC_NCNT);
	integer single_active_REGION20(ATC_NCNT);
	integer single_testreg(ATC_NCNT);
c
c
        common/atc_geom/gatc_r, gatc_phi, gatc_z,
     +   gatc_halfr, gatc_halfphi, gatc_halfz
c
        common/atc_data/gatc_bad, gatc_amp_err, gatc_a0, gatc_tau,
     +   gatc_p, gatc_a1pe, gatc_da1pe, gatc_minDeltaT, gatc_maxDeltaT,
     +   gatc_dcSigmaT, gatc_dcSigmaZ
c
        common/atc_track/atc_ntracks, atc_trk_ncnt, atc_trk_cnts,
     +   atc_illtrack, atc_wls_hit, atc_near_wls, atc_track_end,
     +   atc_trk_aer_len, atc_trk_wls_len,
     +   atc_rin, atc_rout, atc_phiin, atc_phiout, atc_zin, atc_zout,
     +   atc_ph_in, atc_ph_out, atc_rwls, atc_phiwls, atc_zwls,
     +   atc_rB, atc_phiB, atc_zB, atc_rE, atc_phiE, atc_zE,
     +   atc_cnt_ntrk, atc_cnt_tracks
c
        common/atc_rec/atc_npe, atc_dnpe, atc_time, atc_amp, atc_rawtime,
     +   atc_chi2, atc_trg, atc_htyp, atc_eventdamage,
     +   atc_rawDeltaT, atc_deltaT
c
        common/atctrackinfo/atcinfo_t, atcinfo_ncnt, atcinfo_cnt,
     +   atcinfo_npe, atcinfo_npen, atcinfo_wlshit, atcinfo_nearwls,
     +   atcinfo_tlen, atcinfo_pathwls, aerogel_REGION, aerogel_REGION0,
     +   aerogel_REGION5, aerogel_REGION20, active_REGION, active_REGION0,
     +   active_REGION5, active_REGION20, testreg, single_aerogel_REGION,
     +   single_aerogel_REGION0, single_aerogel_REGION5, single_aerogel_REGION20,
     +   single_active_REGION, single_active_REGION0, single_active_REGION5,
     +   single_active_REGION20, single_testreg
