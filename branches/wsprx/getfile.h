#ifndef GETFILE_H
#define GETFILE_H
#include <QString>
#include <QFile>
#include <QDebug>
#include <QMessageBox>
#include "commons.h"

void getfile(QString fname, int ntrperiod);
float gran();
int ptt(int nport, int ntx, int* iptt, int* nopen);

#endif // GETFILE_H
