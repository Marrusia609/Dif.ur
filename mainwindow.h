#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include <QtWidgets>

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
    void MatDate();
    void GlobalCircle();
    void GraphOde45();
    void GraphRungeKutta();
    void GraphDormanPrince();
    void DormanPrince(int *i);
    void RungeKutta(int *i);
    void calcVbasRK(double t0global);
    void calcVbasDP(double t0global);

private:
    Ui::MainWindow *ui;

    QVector <double> Kvect;
    // векторы для выгрузки на графики
    QVector<double> tick;
    QVector<double> dptick;
    QVector<double> tickBas;
    QVector<double> tickBasDP;
    QVector<double> CGB;
    QVector<double> Ipg;
    QVector<double> dpCGB;
    QVector<double> dpIpg;
    QVector<double> cgbMatlab;
    QVector<double> insMatlab;
    QVector<double> grafVbas;
    QVector<double> grafVbasDP;

    double Vbas = 0;
    double VbasDP = 0;
    double bas[7];          // для расчёта Vbas
    double basDP[7];          // для расчёта Vbas
    double Vb;
    double VbDP;
    double t = 0;
    double tDP = 0;
    double dZ[14]; // RK
    double DdZ[14]; // DP
};
#endif // MAINWINDOW_H
