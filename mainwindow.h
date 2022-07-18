#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "Analyzer.h"
#include "parameter.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void updateConsole(const QString& text);

    void setProgressBar(int val);

    void CreateAlert(const QString& text);

    void on_pushButtonSaveFeedback_clicked();

    void on_pushButtonAnalyze_clicked();

    void on_pushButtonStudentModel_clicked();

    void on_pushButtonStudentCenter_clicked();

    void on_pushButtonStudentMidpoint_clicked();

    void on_pushButtonStudentMarginPoints_clicked();

    void on_pushButtonStudentAxialPoints_clicked();

    void on_pushButtonStudentOcclusalPoints_clicked();

    void on_pushButtonStudentGingivaPoints_clicked();

    void on_pushButtonOriginalModel_clicked();

    void on_pushButtonStudentFolder_clicked();

    void on_checkBoxDivision_toggled(bool checked);

    // export results as csv
    void on_pushButtonExport_clicked();

private:
    Ui::MainWindow *ui;
    Analyzer* m_analyzer;
    Parameter param;
};
#endif // MAINWINDOW_H
