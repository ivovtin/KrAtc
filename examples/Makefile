#--------------------------------------------------------------------------
# ��������� 3 ������ ����� �� �������
ifeq (,$(WORKDIR))
WORKDIR=../..
endif

# �ӣ ��� �������� �������� ����� :) - ���������� �������� ��� ����
VERBOSEMAKE=""
# ���� ������ -g ��� ������ �� ����� - ���� ������ � ����
#NDEBUG=""

# ���������� ��� *.� ������
#CC_LOCAL := gcc
# ���������� ��� *.cc � *.cpp ������
#CXX_LOCAL := g++
# ���������� ��� *.f, *.for � *.F ������
#F77_LOCAL := g77
# ������ ��� ����� ������ �� ���������
#LD_LOCAL := g++
LD_LOCAL := g77
# �����������������, ���� ���������� �������� ������ ��� �������
# ��������� � �����������. ������ �� ���������.
#LD_MULTI := ""

# �������������� ����� ��� ����������� C/C++
COPTS = -I. -fno-access-control `root-config --cflags` #-ansi -pedantic
NOOPT = ""
# �������������� ����� ��� ����������� Fortran
#FOPTS = -I.
FOPTS  =  -g -fvxt -Wall -fno-automatic -finit-local-zero \
-fno-second-underscore -ffixed-line-length-120 -Wno-globals \
-DCERNLIB_LINUX -DCERNLIB_UNIX -DCERNLIB_LNX -DCERNLIB_QMGLIBC -DCERNLIB_BLDLIB \
-DCOMPDATE="'$(COMPDATE)'" -I$(CERN)/pro/include -I$(CERN)/pro/include/geant321

# �������������� ����� ��� �������
ROOTLIBDIR := `root-config --libdir`
#LDOPTS = -Xlinker -rpath -Xlinker $(ROOTLIBDIR) -L$(HOME)/lib/X11
LDOPTS = -L$(shell root-config --libdir) -Wl,-rpath,$(shell root-config --libdir)


#���� ���������� ���������� ONLYBINARY, �� ���������� � ������ �����������
ONLYBINARY = ""

# ���� ���������� �������� CERNLIB, �� ����� ��������������� ����
# ��������. ���������� ���������� ����������� � ����� �����. ��
# ��������� ����������� �������� ���������� jetset74 mathlib graflib
# geant321 grafX11 packlib
CERNLIBRARY = ""

# ������ ���������, ���� ��� �� �������� ����������� �����, �������
# ��������� �� ����������. � ����� ������ ���������� �����������������
# CERNLIBRARY
#CERNLIBS = jetset74 mathlib graflib geant lapack3 blas
CERNLIBS = geant321 jetset74 pawlib graflib grafX11 mathlib lapack3 blas
#CERNLIBS = mathlib graflib #geant lapack3 blas

# ��� ���������� ����������� �����
BINDIR := ./

# �������������� ���� (����������� ����� ���)
LIB_LOCAL= `root-config --libs` -lstdc++ -lMinuit -ldl -lg2c -lpq -lcrypt -lbz2
#LIB_LOCAL= `root-config --libs` -lMinuit -lpq -lcrypt -lbz2 -lg2c -ldl -lstdc++

# ���������, ����� ��������� �� ����� ��������
BINARIES = test_atcrecF
##BINARIES = test_atcrec test_atcrecF test_readnat
##BINARIES = test_readnat test_atcrec

# ������, �� ����� ������� ����� ������ ��� �������
# (��� ������ �� ����� �������� � ����������)
# � ����� ���������� ���� ���������� ��� ������
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

# ��������� ������ ����� �� �������
include $(WORKDIR)/KcReleaseTools/rules.mk

