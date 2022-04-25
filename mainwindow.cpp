#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFileDialog>
#include <QDir>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QFuture> // for async computation
#include <QtConcurrent/QtConcurrent> // for async computation

#include "Analyzer.h"
#include "audiorecorder.h"
#include "result.h"

#include <iostream>
#include <stdlib.h>

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow)
  , m_audioRecorder(new AudioRecorder())
  , m_analyzer(new Analyzer())
{
  ui->setupUi(this);
  QObject::connect(m_analyzer,&Analyzer::msgToConsole,this,&MainWindow::updateConsole);
  QObject::connect(m_analyzer, &Analyzer::updateProgressBar, this, &MainWindow::setProgressBar);
  QObject::connect(m_analyzer, &Analyzer::alertToWindow, this, &MainWindow::CreateAlert);
}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::updateConsole(const QString &text)
{
  ui->textBrowserOutput->append(text);
  std::cerr << text.toStdString() << std::endl;
}

void MainWindow::setProgressBar(int val)
{
  ui->progressBar->setValue(val);
}

void MainWindow::CreateAlert(const QString &text)
{
  QMessageBox::information(this, "Alert", text);
}

void MainWindow::on_pushButtonRecordAudio_clicked()
{
  m_audioRecorder->show();
}

void MainWindow::on_pushButtonSaveFeedback_clicked()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save Feedback", QDir::homePath(), tr("Text files (*.txt);;All Files (*)"));
  if (fileName.isEmpty()) {
      return;
    }
  else {
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly)) {
      QMessageBox::information(this, tr("Unable to open file"),
                               file.errorString());
      return;
    }
    QTextStream out(&file);
    QString text = ui->textEdit->toPlainText();
    out << text;
    file.flush();
    file.close();
  }
}

void MainWindow::on_pushButtonAnalyze_clicked()
{
  if (ui->radioButtonAutoAlignment->isChecked()) {
    param.useManualTransform = false;
  } else {
    bool ok = false;
    param.useManualTransform = true;
    param.transformMatrix[0][0] = ui->lineEdit_00->text().toDouble(&ok);
    param.transformMatrix[0][1] = ui->lineEdit_01->text().toDouble(&ok);
    param.transformMatrix[0][2] = ui->lineEdit_02->text().toDouble(&ok);
    param.transformMatrix[0][3] = ui->lineEdit_03->text().toDouble(&ok);
    param.transformMatrix[1][0] = ui->lineEdit_10->text().toDouble(&ok);
    param.transformMatrix[1][1] = ui->lineEdit_11->text().toDouble(&ok);
    param.transformMatrix[1][2] = ui->lineEdit_12->text().toDouble(&ok);
    param.transformMatrix[1][3] = ui->lineEdit_13->text().toDouble(&ok);
    param.transformMatrix[2][0] = ui->lineEdit_20->text().toDouble(&ok);
    param.transformMatrix[2][1] = ui->lineEdit_21->text().toDouble(&ok);
    param.transformMatrix[2][2] = ui->lineEdit_22->text().toDouble(&ok);
    param.transformMatrix[2][3] = ui->lineEdit_23->text().toDouble(&ok);
    param.transformMatrix[3][0] = ui->lineEdit_30->text().toDouble(&ok);
    param.transformMatrix[3][1] = ui->lineEdit_31->text().toDouble(&ok);
    param.transformMatrix[3][2] = ui->lineEdit_32->text().toDouble(&ok);
    param.transformMatrix[3][3] = ui->lineEdit_33->text().toDouble(&ok);
  }

  if (ui->checkBoxDivision->isChecked()) {
    param.divisionEnabled = true;
  } else {
    param.divisionEnabled = false;
  }

  m_analyzer->param = param;
  QtConcurrent::run(this->m_analyzer, &Analyzer::analyze);
}


void MainWindow::on_pushButtonStudentFolder_clicked()
{
    QString dir = QFileDialog::getExistingDirectory(this, "Open Student Folder", QDir::homePath(), QFileDialog::ShowDirsOnly
                                                    | QFileDialog::DontResolveSymlinks);
    if (dir != "") {
      ui->lineEditStudentFolder->setText(dir);
      string base = dir.toStdString();
      param.studentModel = base + "/model.off";
      param.studentCenterPoint = base + "/center_point.pp";
      param.studentMidpoint =  base + "/mid_point.pp";
      param.studentMarginPoints =  base + "/margin_points.pp";
      param.studentAxialPoints = base + "/axial_points.pp";
      param.studentOcclusalPoints = base + "/occlusal_points.pp";
      param.studentGingivaPoints = base + "/gingiva_points.pp";
      ui->lineEditStudentModel->setText(QString::fromStdString(param.studentModel));
      ui->lineEditStudentCenter->setText(QString::fromStdString(param.studentCenterPoint));
      ui->lineEditStudentMidpoint->setText(QString::fromStdString(param.studentMidpoint));
      ui->lineEditStudentMarginPoints->setText(QString::fromStdString(param.studentMarginPoints));
      ui->lineEditStudentAxialPoints->setText(QString::fromStdString(param.studentAxialPoints));
      ui->lineEditStudentOcclusalPoints->setText(QString::fromStdString(param.studentOcclusalPoints));
      ui->lineEditStudentGingivaPoints->setText(QString::fromStdString(param.studentGingivaPoints));
    }
}

void MainWindow::on_pushButtonStudentModel_clicked()
{
  QString filename = QFileDialog::getOpenFileName(this, "Open Student Model", QDir::homePath());
  if (filename != "") {
    param.studentModel = filename.toStdString();
    ui->lineEditStudentModel->setText(filename);
  }
}

void MainWindow::on_pushButtonStudentCenter_clicked()
{
  QString filename = QFileDialog::getOpenFileName(this, "Open Student Model Center", QDir::homePath());
  if (filename != "") {
    param.studentCenterPoint = filename.toStdString();
    ui->lineEditStudentCenter->setText(filename);
  }
}

void MainWindow::on_pushButtonStudentMidpoint_clicked()
{
  QString filename = QFileDialog::getOpenFileName(this, "Open Student Model Midpoint", QDir::homePath());
  if (filename != "") {
    param.studentMidpoint = filename.toStdString();
    ui->lineEditStudentMidpoint->setText(filename);
  }
}

void MainWindow::on_pushButtonStudentMarginPoints_clicked()
{
  QString filename = QFileDialog::getOpenFileName(this, "Open Student Margin Points", QDir::homePath());
  if (filename != "") {
    param.studentMarginPoints = filename.toStdString();
    ui->lineEditStudentMarginPoints->setText(filename);
  }
}

void MainWindow::on_pushButtonStudentAxialPoints_clicked()
{
  QString filename = QFileDialog::getOpenFileName(this, "Open Student Axial Points", QDir::homePath());
  if (filename != "") {
    param.studentAxialPoints = filename.toStdString();
    ui->lineEditStudentAxialPoints->setText(filename);
  }
}

void MainWindow::on_pushButtonStudentOcclusalPoints_clicked()
{
  QString filename = QFileDialog::getOpenFileName(this, "Open Student Occlusal Points", QDir::homePath());
  if (filename != "") {
    param.studentOcclusalPoints = filename.toStdString();
    ui->lineEditStudentOcclusalPoints->setText(filename);
  }
}

void MainWindow::on_pushButtonStudentGingivaPoints_clicked()
{
  QString filename = QFileDialog::getOpenFileName(this, "Open Student Gingiva Points", QDir::homePath());
  if (filename != "") {
    param.studentGingivaPoints = filename.toStdString();
    ui->lineEditStudentGingivaPoints->setText(filename);
  }
}


void MainWindow::on_pushButtonOriginalModel_clicked()
{
  QString filename = QFileDialog::getOpenFileName(this, "Open Original Model", QDir::homePath());
  if (filename != "") {
    param.originalModel = filename.toStdString();
    ui->lineEditOriginalModel->setText(filename);
  }
}

void MainWindow::on_pushButtonCSVExportPath_clicked()
{
  QString filename = QFileDialog::getOpenFileName(this, "Select Export Path", QDir::homePath());
  if (filename != "") {
    param.originalModel = filename.toStdString();
    ui->lineEditCSVExportPath->setText(filename);
  }
}


void MainWindow::on_radioButtonManualAlignment_toggled(bool checked)
{
    if (checked) {
      ui->labelTransformationMatrix->setEnabled(true);
      ui->lineEdit_00->setEnabled(true);
      ui->lineEdit_01->setEnabled(true);
      ui->lineEdit_02->setEnabled(true);
      ui->lineEdit_03->setEnabled(true);
      ui->lineEdit_10->setEnabled(true);
      ui->lineEdit_11->setEnabled(true);
      ui->lineEdit_12->setEnabled(true);
      ui->lineEdit_13->setEnabled(true);
      ui->lineEdit_20->setEnabled(true);
      ui->lineEdit_21->setEnabled(true);
      ui->lineEdit_22->setEnabled(true);
      ui->lineEdit_23->setEnabled(true);
      ui->lineEdit_30->setEnabled(true);
      ui->lineEdit_31->setEnabled(true);
      ui->lineEdit_32->setEnabled(true);
      ui->lineEdit_33->setEnabled(true);
    } else {
      ui->labelTransformationMatrix->setDisabled(true);
      ui->lineEdit_00->setDisabled(true);
      ui->lineEdit_01->setDisabled(true);
      ui->lineEdit_02->setDisabled(true);
      ui->lineEdit_03->setDisabled(true);
      ui->lineEdit_10->setDisabled(true);
      ui->lineEdit_11->setDisabled(true);
      ui->lineEdit_12->setDisabled(true);
      ui->lineEdit_13->setDisabled(true);
      ui->lineEdit_20->setDisabled(true);
      ui->lineEdit_21->setDisabled(true);
      ui->lineEdit_22->setDisabled(true);
      ui->lineEdit_23->setDisabled(true);
      ui->lineEdit_30->setDisabled(true);
      ui->lineEdit_31->setDisabled(true);
      ui->lineEdit_32->setDisabled(true);
      ui->lineEdit_33->setDisabled(true);
    }
}

void MainWindow::on_checkBoxDivision_toggled(bool checked)
{
    if (checked) {
      ui->labelStudentCenter->setEnabled(true);
      ui->labelStudentMidpoint->setEnabled(true);
      ui->lineEditStudentCenter->setEnabled(true);
      ui->lineEditStudentMidpoint->setEnabled(true);
      ui->pushButtonStudentCenter->setEnabled(true);
      ui->pushButtonStudentMidpoint->setEnabled(true);
    } else {
      ui->labelStudentCenter->setDisabled(true);
      ui->labelStudentMidpoint->setDisabled(true);
      ui->lineEditStudentCenter->setDisabled(true);
      ui->lineEditStudentMidpoint->setDisabled(true);
      ui->pushButtonStudentCenter->setDisabled(true);
      ui->pushButtonStudentMidpoint->setDisabled(true);
    }
}

void MainWindow::on_checkBoxCSVExport_toggled(bool checked)
{
  if (checked) {
    ui->labelCSVExport->setEnabled(true);
    ui->lineEditCSVExportPath->setEnabled(true);
    ui->pushButtonCSVExportPath->setEnabled(true);
  } else {
    ui->labelCSVExport->setDisabled(true);
    ui->lineEditCSVExportPath->setDisabled(true);
    ui->pushButtonCSVExportPath->setDisabled(true);
  }
}
