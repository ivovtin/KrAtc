        program test_atcrec    !подключение реконструкции амплитуд  !текст программы начинается с оператора заголовка
        implicit none
#include "VDDCRec/ktrrec.inc"
#include "KrAtc/atcrec.inc"
c
        integer lenocc
c        integer atc_info
        integer atc_trreg
        integer RecMode/0/          !режим
c        integer RunNumber/5461/     !номер захода для обработки
c        integer RunNumber/19654/     !номер захода для обработки  19691 - exp
        integer RunNumber/19449/     !номер захода для обработки  19691 - exp
        integer Result
        integer IERR
        character fname*128
        character chDir*80
        character chPrefix*20
        character ext*20
        character offcnt*200
        integer t, i, c

c       chDir='/space/users/skononov/runs/psi2s/'
c        chDir='/space/data2/KEDR_RUNS/runs/'                         !директория где лежит файл с заходом
        chDir='/space/runs/'                         !директория где лежит файл с заходом
        chPrefix='daq'
        ext = '.nat'
c
        write(fname,'(a,a,i6.6,a)') chDir(:lenocc(chDir)),
     &          chPrefix(:lenocc(chPrefix)),RunNumber,ext
        print *,'fname:',fname                                      !печатаем директорию и номер файла с форматом
c
        CALL KEDR_OPEN_NAT(fname(:lenocc(fname)),IERR)              !открываем нат файл
        if(IERR.ne.0) then
          print *, 'No such file'
	  goto 100
        endif
c
        call atc_init                                                    !вызов инициализации реконструкции
        offcnt=''
        print *, lenocc(offcnt)
        do c=1,160                                                        !по 160-ти счетчикам цикл            - изменил(80)
          if(gatc_bad(c).ne.0) write(offcnt(lenocc(offcnt)+1:),'(I3)'), c
        enddo
        print *, 'ATC offline: ', offcnt
c
        call atc_run(%val(RunNumber))                                    !инициализация для захода run
        print *, '1pe calibration:'                                      !вывод 1 ф.э. амплитуд
        print '(20I5)', (c,c=1,20)
        print '(20F5.1)', (gatc_a1pe(c),c=1,20)
        print '(20I5)', (c,c=21,40)
        print '(20F5.1)', (gatc_a1pe(c),c=21,40)
        print '(20I5)', (c,c=41,60)
        print '(20F5.1)', (gatc_a1pe(c),c=41,60)
        print '(20I5)', (c,c=61,80)
        print '(20F5.1)', (gatc_a1pe(c),c=61,80)
        print '(20I5)', (c,c=81,100)                            !добавил
        print '(20F5.1)', (gatc_a1pe(c),c=81,100)
        print '(20I5)', (c,c=101,120)
        print '(20F5.1)', (gatc_a1pe(c),c=101,120)
        print '(20I5)', (c,c=121,140)
        print '(20F5.1)', (gatc_a1pe(c),c=121,140)
        print '(20I5)', (c,c=141,160)
        print '(20F5.1)', (gatc_a1pe(c),c=141,160)
c
 10     CALL KEDR_READ_NAT(IERR)
        if(IERR.ne.0) then
          print *, 'Error reading file'
	  goto 100
        endif
c
        call kDCVDrec(RecMode,Result)
        print *,'Ev:',EvNum,EvNumDAQ,NTRAK,NtrackIP,NtrackR
        print *,' eVertex:',Xvert,Yvert,Zvert                    !координаты вершины
        print *,' eSigma:',sXvert,sYvert,sZvert                  !сигма вершины
c
        call atc_event                                           !Реконструкция одного события
c
        do t=1,NTRAK                                             !цикл с 1 по како-то трек                                   !!!!!!!!
          print *,'  track ',t                                     !печатаем номер трека
          print *,'  p:',PTRAK(t),TrackIP(t),CHTRAK(t),SGTRAK(t)
          print *,'  direction:',VxTRAK(t),VyTRAK(t),VzTRAK(t)     !направление трека
c          print *,'  nATC=',atc_trk_ncnt(t)
c          do i=1,atc_trk_ncnt(t)                                   !цикл по кол-ву счетчиков пересеченных треком
c          c=atc_trk_cnts(i,t)                                    ! list of counter indexes on the track 1..160
c
cc integer*1                                                !если не менять типы данных
cc          if(c<0)  then
cc          c=c+1+128+127
cc	  endif
c	  print *, '    Counter ',c,': A=',atc_npe(c),' pe,  Laer=',atc_trk_aer_len(i,t),' mm'         !выводим номер сч-ка, амплитуду и длину трека в аэрогеле
c            if(atc_trreg(t,i,0,0).ne.0) print *, '    track crosses aerogel region'
c            if(atc_trreg(t,i,1,0).ne.0) print *, '    track crosses active region'
c            if(atc_trreg(t,i,2,0).ne.0) print *, '    track crosses CRT region'
c
c          print *,'Local coord. --> phiin=',atc_phiin(t,i)
c          print *,'Local coord. --> rin=',atc_rin(t,i)*0.1
c
c	  enddo
c
          call atc_info(t-1)
c
          do i=1,atcinfo_ncnt                                   !цикл по кол-ву счетчиков пересеченных треком
          c=atcinfo_cnt(i)                                    ! list of counter indexes on the track 1..160
	  print *,' *******>>>> atcinfo_cnt=',atcinfo_cnt(i),' atcinfo_npe=',atcinfo_npe(c+1)
	  enddo
                   print *,'****** End track *******'
        enddo
                   print *,'***************** End event ******************'
c
        goto 10
 100    call KEDR_CLOSE_NAT(Ierr)             !закрываем нат файл
c        call atc_stop                        !Освобождение используемой памяти по окончании работы
c
        end                                   !заканчиваем программу
c
