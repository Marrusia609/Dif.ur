#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include <QtWidgets>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

extern double t;
// Initial values & variables
extern double Dig;
extern double ti1;
extern double Dbol1;      // доза болюса начальная
extern double Dbol2;      // доза болюса начальная
extern double Dbol3;      // доза болюса начальная
extern double Vbas;
extern double g0;   // plasma glucose, mg/kg START gp?
extern double gt0;  // tissue glucose, mg/kg  START
extern double Ip0;  // portal vein insulin, pmol/kg START
extern double m2;
extern double del;
extern double tm1; // tm
extern double Ti1;
extern double Ti2;
extern double ti2;
extern double Ti3;
extern double ti3;
extern double tms;
extern double Tm;
extern double Vg;

struct DataZ
{
    double t;
    double gp;  // plasma glucose, mg/kg
    double Il;  // liver insulin, pmol/kg
    double Ip;  // plasma insulin, pmol/kg
    double fgut;    // gut carbs, mg
    double fliq;    // vliquid stomach carbs, mg
    double fsol;    // solid stomach carbs, mg
    double gt;  // tissue glucose, mg/kg
    double I1;  // 1st compartment insulin signal, pmol/l
    double Id;  // delayed insulin signal, pmol/l
    double Xt;  // interstitial insulin, pmol/l
    double Ipo; // portal vein insulin, pmol/kg
    double Yt;  // ??? , pmol/kg/min
    double Ii;  // pmol
    double It;  // pmol
};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
