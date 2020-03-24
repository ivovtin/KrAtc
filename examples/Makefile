#--------------------------------------------------------------------------
# следующие 3 строки лучше не трогать
ifeq (,$(WORKDIR))
WORKDIR=../..
endif

# Всё что делается делается вслух :) - желательно оставить как есть
VERBOSEMAKE=""
# Если ключик -g при сборки не нужен - если уверен в себе
#NDEBUG=""

# Компилятор для *.с файлов
#CC_LOCAL := gcc
# Компилятор для *.cc и *.cpp файлов
#CXX_LOCAL := g++
# Компилятор для *.f, *.for и *.F файлов
#F77_LOCAL := g77
# Линкер для этого пакета по умолчанию
#LD_LOCAL := g++
LD_LOCAL := g77
# Раскомментировать, если необходимо задавать линкер для каждого
# бинарника в отдельности. Обычно не требуется.
#LD_MULTI := ""

# Дополнительные опции для компилятора C/C++
COPTS = -I. -fno-access-control `root-config --cflags` #-ansi -pedantic
NOOPT = ""
# Дополнительные опции для компилятора Fortran
#FOPTS = -I.
FOPTS  =  -g -fvxt -Wall -fno-automatic -finit-local-zero \
-fno-second-underscore -ffixed-line-length-120 -Wno-globals \
-DCERNLIB_LINUX -DCERNLIB_UNIX -DCERNLIB_LNX -DCERNLIB_QMGLIBC -DCERNLIB_BLDLIB \
-DCOMPDATE="'$(COMPDATE)'" -I$(CERN)/pro/include -I$(CERN)/pro/include/geant321

# Дополнительные опции для линкера
ROOTLIBDIR := `root-config --libdir`
#LDOPTS = -Xlinker -rpath -Xlinker $(ROOTLIBDIR) -L$(HOME)/lib/X11
LDOPTS = -L$(shell root-config --libdir) -Wl,-rpath,$(shell root-config --libdir)


#Если определена переменная ONLYBINARY, то библиотека в пакете отсутствует
ONLYBINARY = ""

# Если необходимо добавить CERNLIB, то лучше воспользоваться этим
# ключиком. Церновские библиотеки подшиваются в самом конце. По
# умалчанию подшиваются следущие библиотеки jetset74 mathlib graflib
# geant321 grafX11 packlib
CERNLIBRARY = ""

# Список библиотек, если вам не нравится стандартный набор, который
# расширяем по требованию. В любом случае необходимо раскомментировать
# CERNLIBRARY
#CERNLIBS = jetset74 mathlib graflib geant lapack3 blas
CERNLIBS = geant321 jetset74 pawlib graflib grafX11 mathlib lapack3 blas
#CERNLIBS = mathlib graflib #geant lapack3 blas

# где вываливать исполняемые файлы
BINDIR := ./

# Дополнительные либы (вставляются после либ)
LIB_LOCAL= `root-config --libs` -lstdc++ -lMinuit -ldl -lg2c -lpq -lcrypt -lbz2
#LIB_LOCAL= `root-config --libs` -lMinuit -lpq -lcrypt -lbz2 -lg2c -ldl -lstdc++

# Определим, какие программы мы будем собирать
BINARIES = test_atcrecF
##BINARIES = test_atcrec test_atcrecF test_readnat
##BINARIES = test_readnat test_atcrec

# укажем, из каких модулей этого пакета они состоят
# (эти модули не будут включены в библиотеку)
# и какие библиотеки надо подключить при сборке
#test_atcrec_MODULES := test_atcrec
#test_atcrec_LIBS := KaFramework KrAtc ReadNat AppFramework KDB
#test_atcrec_LD := g++

#test_readnat_MODULES := test_readnat
#test_readnat_LIBS := KaFramework KrAtc ReadNat AppFramework KDB
#test_readnat_LD := g++

test_atcrecF_MODULES := test_atcrecF
test_atcrecF_LIBS := KrAtc VDDCRec KEmcRec LKrTools KrToF KsToF VDRec KDB \
KdDCSim VDDCRec KrMu KrVDDCMu FitTools DchGeom AppFramework ReadNat g2c stdc++
test_atcrecF_LD := g77

# VDDCRec KrMu KrVDDCMu

#test_MODULES := test
#test_LIBS := KrAtc AppFramework KDB
#test_LD := g++

# следующую строку лучше не трогать
include $(WORKDIR)/KcReleaseTools/rules.mk

