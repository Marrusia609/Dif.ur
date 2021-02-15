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
    void DormanPrince(int i);
    void RungeKutta();

private:
    Ui::MainWindow *ui;

    QVector <double> Kvect;
    // векторы для выгрузки на графики
    QVector<double> tick;
    QVector<double> tickBas;
    QVector<double> CGB;
    QVector<double> Ipg;
    QVector<double> cgbMatlab;
    QVector<double> insMatlab;

    double Vbas = 0;
    double t = 0;
    double dZ[14];
};
#endif // MAINWINDOW_H
