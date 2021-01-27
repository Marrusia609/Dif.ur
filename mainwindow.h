#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

extern double t;
extern double Dig;
extern double ti1;
extern double Dbol;
extern double Vbas;
extern double mt;
extern double kp1;
extern double g0;
extern double gt0;
extern double Ip0;
extern double b;
extern double c;
extern double Fcns;
extern double m2;
extern double Vmxx;
extern double del;
extern double tm1; // tm
extern double Ti1;
extern double Ti2;
extern double ti2;
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
