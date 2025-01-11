QT       += core gui
QT += charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17
INCLUDEPATH += D:\Qt_project\eigen-3.4.0
INCLUDEPATH += D:\Qt_project\fftw-3.3.5-dll64
LIBS += -L"D:\Qt_project\fftw-3.3.5-dll64" -lfftw3-3

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
    mcthread.cpp \
    random_flow.cpp \
    random_flow_window.cpp \
    set_fdm.cpp \
    set_hydrogeological_parameter.cpp \
    set_new_wave.cpp \
    set_new_wave_h.cpp \
    tchartview.cpp

HEADERS += \
    Random_flow.h \
    mcthread.h \
    random_flow_window.h \
    set_fdm.h \
    set_hydrogeological_parameter.h \
    set_new_wave.h \
    set_new_wave_h.h \
    tchartview.h

FORMS += \
    random_flow_window.ui \
    set_fdm.ui \
    set_hydrogeological_parameter.ui \
    set_new_wave.ui \
    set_new_wave_h.ui


TRANSLATIONS += \
    Random_flow_zh_CN.ts
CONFIG += lrelease
CONFIG += embed_translations

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

