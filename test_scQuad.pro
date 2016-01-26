#-------------------------------------------------
#
# Project created by QtCreator 2015-10-21T21:39:19
#
#-------------------------------------------------

TARGET = beamWaistDrift



SOURCES += \
    statvec.cpp \
    pass_sc_equad.cpp \
    test_scQuad.cpp

HEADERS += \
    statvec.h \
    pass_sc_equad.h \


# remove possible other optimization flags
QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

# add the desired -O3 if not present
QMAKE_CXXFLAGS_RELEASE *= -O3
