///////////////////////////////////////////////////////////////////////////
// Some code in this file and accompanying files is based on work by
// Moe Wheatley, AE4Y, released under the "Simplified BSD License".
// For more details see the accompanying file LICENSE_WHEATLEY.TXT
///////////////////////////////////////////////////////////////////////////

#ifndef PLOTTER_H
#define PLOTTER_H

#include <QtGui>
#include <QFrame>
#include <QImage>
#include <cstring>
#include "commons.h"

#define VERT_DIVS 7	//specify grid screen divisions
#define HORZ_DIVS 20

class CPlotter : public QFrame
{
  Q_OBJECT
public:
  explicit CPlotter(QWidget *parent = 0);
  ~CPlotter();

  QSize minimumSizeHint() const;
  QSize sizeHint() const;
  QColor  m_ColorTbl[256];
  double  m_dialFreq;
  float   m_fSpan;
  qint32  m_TxFreq;
  qint32  m_w;
  qint32  m_BFO;
  int     m_plotZero;
  int     m_plotGain;
  bool    m_bRFscale;
  bool    m_bCumulative;
  bool    m_transmitted;

  void draw(float sw[]);		               //Update the waterfall
  void SetRunningState(bool running);
  void setPlotZero(int plotZero);
  int  getPlotZero();
  void setPlotGain(int plotGain);
  int  getPlotGain();
  int  plotWidth();
  void UpdateOverlay();
  void setDataFromDisk(bool b);
  void setTxFreq(int n, bool bf);
  void DrawOverlay();
  int  TxFreq();
  void setPalette(QString palette);
  void setFsample(int n);
  void setNsps(int ntrperiod, int nsps);
  void SetPercent2DScreen(int percent){m_Percent2DScreen=percent;}

signals:
  void freezeDecode0(int n);
  void freezeDecode1(int n);

protected:
  //re-implemented widget event handlers
  void paintEvent(QPaintEvent *event);
  void resizeEvent(QResizeEvent* event);

private:

  void MakeFrequencyStrs();
  void UTCstr();
  int XfromFreq(float f);
  float FreqfromX(int x);
  qint64 RoundFreq(qint64 freq, int resolution);

  QPixmap m_WaterfallPixmap;
  QPixmap m_2DPixmap;
  QPixmap m_ScalePixmap;
  QPixmap m_OverlayPixmap;
  QSize   m_Size;
  QString m_Str;
  QString m_HDivText[483];

  double  m_fftBinWidth;
  qint64  m_StartFreq;
  qint32  m_dBStepSize;
  qint32  m_hdivs;
  qint32  m_line;
  qint32  m_fSample;
  qint32  m_xClick;
  qint32  m_freqPerDiv;
  qint32  m_nsps;
  qint32  m_Percent2DScreen;
  qint32  m_h;
  qint32  m_h1;
  qint32  m_h2;
  qint32  m_tol;
  qint32  m_TRperiod;
  qint32  m_RFHz;
  bool    m_Running;
  bool    m_paintEventBusy;
  bool    m_dataFromDisk;
  char    m_sutc[5];


private slots:
  void mousePressEvent(QMouseEvent *event);
  void mouseDoubleClickEvent(QMouseEvent *event);
};

#endif // PLOTTER_H
