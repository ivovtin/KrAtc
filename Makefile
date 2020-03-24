# -*- mode: Makefile; -*- ��������� Makefile ���� ��� Emacs
#--------------------------------------------------------------------------
#      $Id: Makefile,v 1.17 2019/06/25 06:49:08 ovtin Exp $
#
# Description:
#      Makefile for KcTemplate package
#
# Environment:
#      Software developed for the KEDR Detector at BINP.
#
# Author:
#      Baldin Evgueni Mihilovich (E.M.Baldin@inp.nsk.su)
#      Kononov Sergey (S.A.Kononov@inp.nsk.su)
# Copyright Information:
#      Copyright (C) 2000-2004  Budker's Institute of Nuclear Physics
#
#--------------------------------------------------------------------------
# ��������� 3 ������ ����� �� �������
ifeq (,$(WORKDIR))
WORKDIR = $(shell pushd ../ 1>/dev/null && pwd )
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
F77_LOCAL := g77
# ������ ��� ����� ������ �� ���������
#LD_LOCAL := g++
# �����������������, ���� ���������� �������� ������ ��� �������
# ��������� � �����������. ������ �� ���������.
#LD_MULTI := ""

# �������������� ����� ��� ����������� C/C++
COPTS  = `root-config --cflags` #-ansi -pedantic

#ATC reconstruction customization
ifdef ATC_NOFIT
COPTS += -DATC_NOFIT=1
endif
ifdef ATC_NOROOT
COPTS += -DATC_NOROOT=1
endif

# �������������� ����� ��� ����������� Fortran
FOPTS  = -g -fvxt -Wall -fno-automatic -finit-local-zero \
-fno-second-underscore -ffixed-line-length-120 -Wno-globals \
-DCERNLIB_LINUX -DCERNLIB_UNIX -DCERNLIB_LNX -DCERNLIB_QMGLIBC -DCERNLIB_BLDLIB \
-DCOMPDATE="'$(COMPDATE)'" -I$(CERN)/pro/include -I$(CERN)/pro/include/geant321

# ���������� ����������� - ������-�� ����� ��� �������������
NOOPT = ""

# �������������� ����� ��� �������
LDOPTS = #-shared #���� ������� ������������ ����������� ����������

# ���� ���������� �������� CERNLIB, �� ����� ��������������� ����
# ��������. ���������� ���������� ����������� � ����� �����. ��
# ��������� ����������� �������� ���������� jetset74 mathlib graflib
# geant321 grafX11 packlib
CERNLIBRARY = ""

# ������ ���������, ���� ��� �� �������� ����������� �����, �������
# ��������� �� ����������. � ����� ������ ���������� �����������������
# CERNLIBRARY
CERNLIBS = geant jetset74 pawlib graflib grafX11 mathlib lapack3 blas

# ���� ���������� ���������� ONLYBINARY, �� ���������� � ������
# ����������� �, ��������������, �� ����������
#ONLYBINARY=""
# �������������� ���������� ������.
BINDIR=./
#����������������, ���� ������� ������� ������������ ����������
#LIB_SHARED = ""

# �������������� ����������. ����������� �� ��� ������.
#LIB_LOCAL=
#IB_LOCAL= `root-config --libs` -lMinuit -lpq -lcrypt -lbz2 -ldl -lg2c
LIB_LOCAL= `root-config --libs` -lstdc++ -lg2c -lpq -lcrypt -lbz2

# ���������, ����� ��������� �� ����� ��������.  ����� �������
# ��������� ������ ������������ ��� ������ ������, ����� �� ��������
# bin - ��� ����������, ������� ������� ��� ������� � ������������
# ������� ���������� � ���������� example. ��� ������ ������ ������ ��
# ���������� ���������� ��� ������� ��� ���������� ������� make,
# ��������: make templateTest
BINARIES =

# ������, �� ����� ������� ����� ������ ������� ������������� �
# BINARIES ��������� (��� ������ �� ����� �������� � ����������),
# ����� ���������� ���� ���������� ��� ������ � ����� ����������
# ������� ������ ������ (����������������� LD_MULTI).

# ��������� ������ ����� �� �������
#-----------------------------------
include $(WORKDIR)/KcReleaseTools/rules.mk
#-----------------------------------

#-----------------------------------
# ���� ������� ������� ���-�� �������������, �� �� �����������
#-----------------------------------
# ���� ���� ������� ���������� �����.
# ����� ����������� �������� ������� ��������
# ������ ��������� <Tab>.

.PHONY: example cleanobj reduced

REDUCEDLIBNAME := $(PACKAGE)-reduced

example:
	$(MAKE) -C examples

cleanobj:
	@echo "Removing fickle object files.."
	$(CMDPREFIX)rm -f $(TEMPDIR)/.depend $(TEMPDIR)/{AtcRec,AtcHit,Atc*Fitter}.o $(CMDSUFFIX)

reduced:
	@echo -e "\nMake reduced version of library.."
	@$(MAKE) ATC_NOFIT=1 ATC_NOROOT=1 LIBNAME=$(REDUCEDLIBNAME) cleanobj lib

#���� ������ ������� ���-�� ����� �������
#����������� ��� ���������� make ��� ����������
PRE_MAKE: cleanobj


#���� ������ ������� ���-�� ����� ������ ������
#����������� ��� ���������� make ��� ����������
POST_MAKE: reduced

# ���� ������� ������� ���-�� �ݣ - �������� ������ ������������
#doc:
#	@echo "making docs"
#	@a2ps -X koi8-r Makefile -o Makefile.ps
#	@perldoc -t fixinc.pl > fixinc.txt
# DO NOT DELETE
