#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#ifdef QT5
#include <QtWidgets>
#else
#include <QtGui>
#endif
#include <QThread>
#include <QTimer>
#include <QDateTime>
#include <QList>
#include <QAudioDeviceInfo>
#include <QScopedPointer>

#include "soundin.h"
#include "soundout.h"
#include "commons.h"
#include "psk_reporter.h"
#include "rigclass.h"
#include "signalmeter.h"
#include "logbook/logbook.h"
#include "Detector.hpp"
#include "Modulator.hpp"
#include "decodedtext.h"


#define NUM_JT65_SYMBOLS 126
#define NUM_JT9_SYMBOLS 85
#define NUM_CW_SYMBOLS 250
#define TX_SAMPLE_RATE 48000

extern int itone[NUM_JT65_SYMBOLS]; //Audio tones for all Tx symbols
extern int icw[NUM_CW_SYMBOLS];	    //Dits for CW ID


//--------------------------------------------------------------- MainWindow
namespace Ui {
    class MainWindow;
}

class QSettings;
class WideGraph;
class LogQSO;

class MainWindow : public QMainWindow
{
  Q_OBJECT

// Multiple instances: call MainWindow() with *thekey
public:
  explicit MainWindow(QSettings *, QSharedMemory *shdmem, QString *thekey,
                      qint32 fontSize2, qint32 fontWeight2, unsigned downSampleFactor,
                      QWidget *parent = 0);
  ~MainWindow();

public slots:
  void showSoundInError(const QString& errorMsg);
  void showSoundOutError(const QString& errorMsg);
  void showStatusMessage(const QString& statusMsg);
  void dataSink(qint64 frames);
  void diskDat();
  void diskWriteFinished();
  void freezeDecode(int n);
  void guiUpdate();
  void doubleClickOnCall(bool shift, bool ctrl);
  void doubleClickOnCall2(bool shift, bool ctrl);
  void readFromStdout();
  void readFromStderr();
  void jt9_error(QProcess::ProcessError);
  void setXIT(int n);
  void setFreq4(int rxFreq, int txFreq);

protected:
  virtual void keyPressEvent( QKeyEvent *e );
  void  closeEvent(QCloseEvent*);
  virtual bool eventFilter(QObject *object, QEvent *event);

private slots:
  void on_tx1_editingFinished();
  void on_tx2_editingFinished();
  void on_tx3_editingFinished();
  void on_tx4_editingFinished();
  void on_tx5_editingFinished();
  void on_tx6_editingFinished();
  void on_actionDeviceSetup_triggered();
  void on_monitorButton_clicked();
  void on_actionExit_triggered();
  void on_actionAbout_triggered();
  void OnExit();
  void on_autoButton_clicked();
  void on_stopTxButton_clicked();
  void on_stopButton_clicked();
  void on_actionOnline_Users_Guide_triggered();
  void on_actionWide_Waterfall_triggered();
  void on_actionOpen_triggered();
  void on_actionOpen_next_in_directory_triggered();
  void on_actionDecode_remaining_files_in_directory_triggered();
  void on_actionDelete_all_wav_files_in_SaveDir_triggered();
  void on_actionNone_triggered();
  void on_actionSave_all_triggered();
  void on_actionKeyboard_shortcuts_triggered();
  void on_actionSpecial_mouse_commands_triggered();
  void on_DecodeButton_clicked();
  void decode();
  void decodeBusy(bool b);
  void on_EraseButton_clicked();
  void on_txb1_clicked();
  void on_txFirstCheckBox_stateChanged(int arg1);
  void set_ntx(int n);
  void on_txb2_clicked();
  void on_txb3_clicked();
  void on_txb4_clicked();
  void on_txb5_clicked();
  void on_txb6_clicked();
  void on_lookupButton_clicked();
  void on_addButton_clicked();
  void on_dxCallEntry_textChanged(const QString &arg1);
  void on_dxGridEntry_textChanged(const QString &arg1);
  void on_genStdMsgsPushButton_clicked();
  void on_logQSOButton_clicked();
  void on_actionJT9_1_triggered();
  void on_actionJT65_triggered();
  void on_actionJT9_JT65_triggered();
  void on_TxFreqSpinBox_valueChanged(int arg1);
  void on_actionSave_decoded_triggered();
  void on_actionQuickDecode_triggered();
  void on_actionMediumDecode_triggered();
  void on_actionDeepestDecode_triggered();
  void on_inGain_valueChanged(int n);
  void bumpFqso(int n);
  void on_actionMonitor_OFF_at_startup_triggered();
  void on_actionErase_ALL_TXT_triggered();
  void on_actionErase_wsjtx_log_adi_triggered();
  void showMacros(const QPoint& pos);
  void onPopup1();
  void onPopup2();
  void onPopup3();
  void onPopup4();
  void onPopup5();
  void onPopup6();
  void onPopup7();
  void onPopup8();
  void onPopup9();
  void onPopup10();
  void on_actionConvert_JT9_x_to_RTTY_triggered(bool checked);
  void on_actionLog_dB_reports_to_Comments_triggered(bool checked);
  void startTx2();
  void stopTx();
  void stopTx2();
  void on_actionPrompt_to_log_QSO_triggered(bool checked);
  void on_actionBlank_line_between_decoding_periods_triggered(bool checked);
  void on_actionEnable_DXCC_entity_triggered(bool checked);
  void on_actionClear_DX_Call_and_Grid_after_logging_triggered(bool checked);
  void on_actionDisplay_distance_in_miles_triggered(bool checked);
  void on_pbCallCQ_clicked();
  void on_pbAnswerCaller_clicked();
  void on_pbSendRRR_clicked();
  void on_pbAnswerCQ_clicked();
  void on_pbSendReport_clicked();
  void on_pbSend73_clicked();
  void on_rbGenMsg_toggled(bool checked);
  void on_rbFreeText_toggled(bool checked);
  void on_freeTextMsg_editingFinished();
  void on_actionDouble_click_on_call_sets_Tx_Enable_triggered(bool checked);
  void on_rptSpinBox_valueChanged(int n);
  void on_action_73TxDisable_triggered(bool checked);
  void on_actionRunaway_Tx_watchdog_triggered(bool checked);
  void killFile();
  void on_tuneButton_clicked();
  void on_actionAllow_multiple_instances_triggered(bool checked);
  void on_pbR2T_clicked();
  void on_pbT2R_clicked();
  void acceptQSO2(bool accepted);
  void on_bandComboBox_activated(int index);
  void on_readFreq_clicked();
  void on_pbTxMode_clicked();
  void on_RxFreqSpinBox_valueChanged(int n);
  void on_cbTxLock_clicked(bool checked);
  void on_actionTx2QSO_triggered(bool checked);  
  void on_cbPlus2kHz_toggled(bool checked);
  void on_outAttenuation_valueChanged (int);
  void on_actionAstronomical_data_triggered();
  void on_actionShort_list_of_add_on_prefixes_and_suffixes_triggered();
  void getpfx();
  void on_actionJT9W_1_triggered();

private:
  Q_SIGNAL void startAudioOutputStream (QAudioDeviceInfo, unsigned channels, unsigned msBuffered);
  Q_SIGNAL void stopAudioOutputStream ();

  Q_SIGNAL void startAudioInputStream (QAudioDeviceInfo const&, unsigned channels, int framesPerBuffer, QIODevice * sink, unsigned downSampleFactor);
  Q_SIGNAL void stopAudioInputStream ();

  Q_SIGNAL void startDetector (AudioDevice::Channel);
  Q_SIGNAL void detectorSetMonitoring (bool);
  Q_SIGNAL void detectorClose ();

  Q_SIGNAL void finished ();
  Q_SIGNAL void muteAudioOutput (bool = true);
  Q_SIGNAL void transmitFrequency (unsigned);
  Q_SIGNAL void endTransmitMessage ();
  Q_SIGNAL void tune (bool = true);
  Q_SIGNAL void sendMessage (unsigned symbolsLength, double framesPerSymbol, unsigned frequency, AudioDevice::Channel, bool synchronize = true, double dBSNR = 99.);
  Q_SIGNAL void outAttenuationChanged (qreal);

private:
  QSettings * m_settings;
    Ui::MainWindow *ui;

    // other windows
    QScopedPointer<WideGraph> m_wideGraph;
    QScopedPointer<LogQSO> m_logDlg;

    double  m_dialFreq;
    double  m_toneSpacing;
    double  m_fSpread;

    float   m_DTmin;
    float   m_DTmax;

    qint64  m_msErase;
    qint64  m_secBandChanged;

    qint32  m_idInt;
    qint32  m_waterfallAvg;
    qint32  m_pttMethodIndex;
    qint32  m_ntx;
    qint32  m_pttPort;
    qint32  m_timeout;
    qint32  m_rxFreq;
    qint32  m_txFreq;
    int     m_XIT;
    qint32  m_setftx;
    qint32  m_ndepth;
    qint32  m_sec0;
    qint32  m_RxLog;
    qint32  m_nutc0;
    qint32  m_nrx;
    qint32  m_hsym;

    Detector m_detector;
    QAudioDeviceInfo m_audioInputDevice;
    AudioDevice::Channel m_audioInputChannel;
    SoundInput m_soundInput;

    Modulator m_modulator;
    QAudioDeviceInfo m_audioOutputDevice;
    AudioDevice::Channel m_audioOutputChannel;
    SoundOutput m_soundOutput;
    QThread m_audioThread;

    qint32  m_TRperiod;
    qint32  m_nsps;
    qint32  m_hsymStop;
    qint32  m_len1;
    qint32  m_inGain;
    qint32  m_nsave;
    qint32  m_catPortIndex;
    qint32  m_rig;
    qint32  m_rigIndex;
    qint32  m_serialRate;
    qint32  m_serialRateIndex;
    qint32  m_dataBits;
    qint32  m_dataBitsIndex;
    qint32  m_stopBits;
    qint32  m_stopBitsIndex;
    qint32  m_handshakeIndex;
    qint32  m_ncw;
    qint32  m_secID;
    qint32  m_band;
    qint32  m_repeatMsg;
    qint32  m_watchdogLimit;
    qint32  m_poll;
    qint32  m_fMax;
    qint32  m_bad;
    qint32  m_EMEbandIndex;
    qint32  m_toneMultIndex;
    qint32  m_astroFont;

    bool    m_monitoring;
    bool    m_btxok;		//True if OK to transmit
    bool    m_btxMute;		//True if transmit should be muted
    bool    m_transmitting;
    bool    m_diskData;
    bool    m_loopall;
    bool    m_decoderBusy;
    bool    m_txFirst;
    bool    m_auto;
    bool    m_restart;
    bool    m_startAnother;
    bool    m_saveDecoded;
    bool    m_saveAll;
    bool    m_widebandDecode;
    bool    m_call3Modified;
    bool    m_dataAvailable;
    bool    m_killAll;
    bool    m_bdecoded;
    bool    m_monitorStartOFF;
    bool    m_pskReporter;
    bool    m_pskReporterInit;
    bool    m_noSuffix;
    bool    m_toRTTY;
    bool    m_dBtoComments;
    bool    m_catEnabled;
    bool    m_After73;
    bool    m_promptToLog;
    bool    m_blankLine;
    bool    m_insertBlank;
    bool    m_displayDXCCEntity;
    bool    m_clearCallGrid;
    bool    m_bMiles;
    bool    m_decodedText2;
    bool    m_freeText;
    bool    m_quickCall;
    bool    m_73TxDisable;
    bool    m_sent73;
    bool    m_runaway;
    bool    m_tune;
    bool    m_bRigOpen;
    bool    m_bMultipleOK;
    bool    m_bDTR;
    bool    m_bRTS;
    bool    m_pttData;
    bool    m_dontReadFreq;
    bool    m_lockTxFreq;
    bool    m_tx2QSO;
    bool    m_CATerror;
    bool    m_bSplit;
    bool    m_bXIT;
    bool    m_plus2kHz;
    bool    m_bAstroData;

    char    m_decoded[80];

    float   m_pctZap;

    QLabel* lab1;                            // labels in status bar
    QLabel* lab2;
    QLabel* lab3;
    QLabel* lab4;
    QLabel* lab5;
    QLabel* lab6;

    QMessageBox msgBox0;

    QFuture<void>* future1;
    QFuture<void>* future2;
    QFuture<void>* future3;
    QFutureWatcher<void>* watcher1;
    QFutureWatcher<void>* watcher2;
    QFutureWatcher<void>* watcher3;

    QProcess proc_jt9;

    QTimer m_guiTimer;
    QTimer* ptt1Timer;                 //StartTx delay
    QTimer* ptt0Timer;                 //StopTx delay
    QTimer* logQSOTimer;
    QTimer* killFileTimer;
    QTimer* tuneButtonTimer;

    QString m_path;
    QString m_pbdecoding_style1;
    QString m_pbmonitor_style;
    QString m_pbAutoOn_style;
    QString m_pbTune_style;
    QString m_myCall;
    QString m_myGrid;
    QString m_baseCall;
    QString m_hisCall;
    QString m_hisGrid;
    QString m_appDir;
    QString m_saveDir;
    QString m_dxccPfx;
    QString m_palette;
    QString m_dateTime;
    QString m_mode;
    QString m_modeTx;
    QString m_fname;
    QString m_rpt;
    QString m_rptSent;
    QString m_rptRcvd;
    QString m_qsoStart;
    QString m_qsoStop;
    QString m_catPort;
    QString m_handshake;
    QString m_cmnd;
    QString m_msgSent0;
    QString m_fileToSave;
    QString m_azelDir;

    QStringList m_macro;
    QStringList m_dFreq;           // per band frequency in MHz as a string
    QStringList m_antDescription;  // per band antenna description
    QStringList m_bandDescription; // per band description
    QStringList m_prefix;
    QStringList m_suffix;

    QHash<QString,bool> m_pfx;
    QHash<QString,bool> m_sfx;

    QDateTime m_dateTimeQSO;
    QRect   m_astroGeom;

    QSharedMemory *mem_jt9;
 // Multiple instances:
    QString       *mykey_jt9;
    PSK_Reporter *psk_Reporter;
    SignalMeter *signalMeter;
    LogBook m_logBook;
    DecodedText m_QSOText;
    unsigned m_msAudioOutputBuffered;
    unsigned m_framesAudioInputBuffered;
    unsigned m_downSampleFactor;
    QThread::Priority m_audioThreadPriority;


//---------------------------------------------------- private functions
    void readSettings();
    void writeSettings();
    void createStatusBar();
    void updateStatusBar();
    void msgBox(QString t);
    void genStdMsgs(QString rpt);
    void lookup();
    void ba2msg(QByteArray ba, char* message);
    void msgtype(QString t, QLineEdit* tx);
    void stub();
    void statusChanged();
    void dialFreqChanged2(double f);
    void freeText();
    void rigOpen();
    void pollRigFreq();
    bool gridOK(QString g);
    bool shortList(QString callsign);
    QString baseCall(QString t);
    void transmit (double snr = 99.);
};

extern void getfile(QString fname, int ntrperiod);
extern void savewav(QString fname, int ntrperiod);
extern int killbyname(const char* progName);
extern void getDev(int* numDevices,char hostAPI_DeviceName[][50],
                   int minChan[], int maxChan[],
                   int minSpeed[], int maxSpeed[]);
extern int ptt(int nport, int ntx, int* iptt, int* nopen);

extern "C" {
//----------------------------------------------------- C and Fortran routines
void symspec_(int* k, int* ntrperiod, int* nsps, int* ingain, int* nflatten,
              float* px, float s[], float* df3, int* nhsym, int* npts8);

void genjt9_(char* msg, int* ichk, char* msgsent, int itone[],
             int* itext, int len1, int len2);

void gen65_(char* msg, int* ichk, char* msgsent, int itone[],
             int* itext, int len1, int len2);

bool stdmsg_(const char* msg, int len);

void azdist_(char* MyGrid, char* HisGrid, double* utch, int* nAz, int* nEl,
             int* nDmiles, int* nDkm, int* nHotAz, int* nHotABetter,
             int len1, int len2);

void morse_(char* msg, int* icw, int* ncw, int len);

int ptt_(int nport, int ntx, int* iptt, int* nopen);

}

#endif // MAINWINDOW_H
