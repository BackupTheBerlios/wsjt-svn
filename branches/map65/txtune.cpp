#include "txtune.h"
#include "ui_txtune.h"
#include <QDebug>

extern int txPower;
extern int iqAmp;
extern int iqPhase;
extern bool bTune;

TxTune::TxTune(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::TxTune)
{
    ui->setupUi(this);
}

TxTune::~TxTune()
{
    delete ui;
}

void TxTune::on_pwrSlider_valueChanged(int n)
{
  txPower=0.01*n;
}

void TxTune::on_ampSlider_valueChanged(int n)
{
  m_iqAmp1=n;
  iqAmp=10*m_iqAmp1 + m_iqAmp2;
  QString t;
  t.sprintf("%.4f",1.0 + 0.0001*iqAmp);
  ui->labAmpReal->setText(t);
}

void TxTune::on_fineAmpSlider_valueChanged(int n)
{
  m_iqAmp2=n;
  iqAmp=10*m_iqAmp1 + m_iqAmp2;
  QString t;
  t.sprintf("%.4f",1.0 + 0.0001*iqAmp);
  ui->labAmpReal->setText(t);}

void TxTune::on_phaSlider_valueChanged(int n)
{
  m_iqPha1=n;
  iqPhase=10*m_iqPha1 + m_iqPha2;
  QString t;
  t.sprintf("%.2f",0.01*iqPhase);
  ui->labPhaReal->setText(t);
}

void TxTune::on_finePhaSlider_valueChanged(int n)
{
  m_iqPha2=n;
  iqPhase=10*m_iqPha1 + m_iqPha2;
  QString t;
  t.sprintf("%.2f",0.01*iqPhase);
  ui->labPhaReal->setText(t);
}

void TxTune::set_iqAmp(int n)
{
  m_iqAmp1=n/10;
  m_iqAmp2=n%10;
  ui->ampSlider->setValue(m_iqAmp1);
  ui->fineAmpSlider->setValue(m_iqAmp2);
}

void TxTune::set_iqPhase(int n)
{
  m_iqPha1=n/10;
  m_iqPha2=n%10;
  ui->phaSlider->setValue(m_iqPha1);
  ui->finePhaSlider->setValue(m_iqPha2);
}

void TxTune::on_cbTxImage_toggled(bool b)
{
  ui->ampSlider->setEnabled(b);
  ui->fineAmpSlider->setEnabled(b);
  ui->labAmp->setEnabled(b);
  ui->labFineAmp->setEnabled(b);
  ui->phaSlider->setEnabled(b);
  ui->finePhaSlider->setEnabled(b);
  ui->labPha->setEnabled(b);
  ui->labFinePha->setEnabled(b);
}

void TxTune::on_pbTune_clicked()
{
  bTune = !bTune;
  if(bTune) {
    QString style="QPushButton{background-color: red;}";
    ui->pbTune->setStyleSheet(style);
  } else {
    ui->pbTune->setStyleSheet("");
  }
}
