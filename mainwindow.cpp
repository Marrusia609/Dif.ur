 #include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <vector>
#include <QVector>
#include <math.h>
#include "qcustomplot.h"

using namespace std;
                        // + значит что значение перенесено с новых файлов
double t = 0.0;         // время начинается с нуля

    /* meal data */
double tm1 = 60; // tm  // время приёма пищи  +
double Tm = 20;         // длительность приёма пищи +
double Dig = 1e-05;     // масса углеводов начальная

    /* bolus */
double Ti1 = 10;        // наверное?
double Ti2 = 10;        // что из этого длительность ввода а что время ввода
double ti2 = 60;
double ti1 = 30;
double Ti3 = 10;
double ti3 = 10;

double del = 0.4;       // result of bolus calculation +

double Dbol1 = 1.8;      // доза болюса начальная
double Dbol2 = 1.8;      // доза болюса начальная
double Dbol3 = 1.8;      // доза болюса начальная

    /* bazal ins */
double Vbas = 1.2238;   // доза базального инсулина  +

double Vg = 1.8;        // plasma per BW, dl/kg +

double m2 = 90;         // + М ?

/* Разные названия переменных с начальным массивом Z ? */
double g0 = 0;           // gp?
double gt0 = 169.72;     // +
double Ip0 = 3.08; //5.2040;      // 0.15??  +

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QVector <double> Kvect;

    /* Обработка трекера */
    QTextStream in;                         // поток данных
    QString read;

    int ind = 0;                            // индекс для записи в массивы

    QString constPat ("D:/Qt/pr/untitled/Dif.ur/matlab/PatConst.csv");
    QFile file(constPat);
    if (!file.open(QFile::ReadOnly)){
        qDebug() << "File not exists";
    }
    else{

        in.setDevice(&file);
        in >> read;
        file.close();

        double tr=0;                        // переменные для записи данных из csv файла в разные массивы

        for ( int i=0; i < 37; i++ ){

                ind=read.indexOf(",");
                tr=read.left(ind).toDouble();
                Kvect.append(tr);                // массив времени
                read.remove(0,ind+1);

        }
    }

    double  a = 0.0,    // начальное время
            b = 720.0;    // конечное время
    double  step = 5;   // шаг НЕ ШАГ)) периодичность!


    vector<DataZ> data;
    DataZ start;
    start.t = a;

    int stepCount = (b / step) + 1;
    double currentTime = a + step;

    double K[7];
    double dz;
    double err;
    double s;

    // убрал домножение на шаг в следуещей точке
    double h = 1.5;
    double hmin = 0.01;
    double hmax = 1.5;
    double eps=0.000001; //error allowance in one step calculation.
    double h1 = h;
    double h2 = h;
    double h3 = h;
    double h4 = h;
    double h5 = h;
    double h6 = h;
    double h7 = h;
    double h8 = h;
    double h9 = h;
    double h10 = h;
    double h11 = h;
    double h12 = h;
    double h13 = h;
    double h14 = h;
    double hopt1;
    double hopt2;
    double hopt3;
    double hopt4;
    double hopt5;
    double hopt6;
    double hopt7;
    double hopt8;
    double hopt9;
    double hopt10;
    double hopt11;
    double hopt12;
    double hopt13;
    double hopt14;

    double OB = 6.63;       // result of bolus calculation +

    Dig = 1176*m2+1; // вместо m2 массу углеводов
    Dbol1=del*OB;
    Dbol2=(1-del)*OB;
    Dbol3=0;

    double bas[7];
    bas[0] = 1.2238;
    bas[1] = 140.0;
    bas[2] = 110.0;
    bas[3] = 70.0;
    bas[4] = 1.5;
    //bas[5] = bas[0]*(1+bas[4]);
    bas[6] = 180.0;
    double Vb = bas[0];
    Vbas = Vb;
    QVector <double> grafVbas;
    grafVbas.append(Vbas);

    /************************************************                        Дано                                           *************************************************/

    /* bgDynam */               // присвоение стартовых значений в структуру для расчёта
    double Ginit=122.00;        // начальный уровень гликемии
    g0 = Ginit*1.8;                 // +
    start.gp = g0;                  // +
    start.Il = 2.61;                // +
    start.Ip = 5.2040;                 //55.12/20 +
    start.fgut = 0;                 // mg +
    start.fliq = 0;                 // mg +
    start.fsol = 0;                 // mg +
    start.gt = gt0;                 // +
    start.I1 = 104.08;              // +
    start.Id = 104.08;              // +
    start.Xt = 0;                   // +
    start.Ipo = 3.08;               // +
    start.Yt = 0;                   // +
    start.Ii = 4120;                // +
    start.It = 10830;               // +
    data.push_back(start);

    // Fluctuations
    double vbas;
    double bol1;
    double bol2;
    double bol3;
    double vbol;

    double vm;

    // выходные значения
    double dgp;
    double dgt;
    double dI1;
    double dXt;
    double dIl;
    double dIp;
    double dfgut;
    double dfliq;
    double dfsol;
    double dIpo;
    double dYt;
    double dIt;
    double dIi;
    double dId;
    // присваиваем стартовые значения
    double gp = start.gp;
    double Id = start.Id;
    double Ipo = start.Ipo;
    double fgut = start.fgut;
    double fsol = start.fsol;
    double fliq = start.fliq;
    double gt = start.gt;
    double Xt = start.Xt;
    double I1 = start.I1;
    double Ip = start.Ip;
    double It = start.It;
    double Il = start.Il;
    double Yt = start.Yt;
    double Ii = start.Ii;
    /*************************************************                      Цикл решения ДУ                   ********************************************************/

    for (int j = 1; j <= stepCount; j++)
    {

        //hmin = 16*eps*t;
        /* */
        // Fluctuations
        vbas=Vbas*100; // pmol/min
        //точки из матлаб убрал в последующих двух уравнениях после второй скобки множителя
        bol1=1/Ti1*Dbol1*(1./(1+exp(-3*(t+10-ti1))))*(1./(1+exp(-3*(-10+ti1-t+Ti1))));
        bol2=1/Ti2*Dbol2*(1./(1+exp(-3*(t+10-ti2))))*(1./(1+exp(-3*(-10+ti2-t+Ti2))));
        bol3=1/Ti3*Dbol3*(1./(1+exp(-3*(t+10-ti3))))*(1./(1+exp(-3*(-10+ti3-t+Ti3))));
        vbol=6000*(bol1+bol2+bol3);

        vm=Dig/Tm*(1./(1+exp(-3*(t-tm1))))*(1./(1+exp(-3*(-(t-tm1-Tm)))));

        /* ДУ */

        // Решаем dgp через Рунге-Кутта
        //dgp=EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*fgut-Kvect.at(26)-Kvect.at(31)*(gp-Kvect.at(32))*Heavi1-Kvect.at(24)*gp+Kvect.at(25)*gt; // plasma glucose mg/dl
        // dgp = f(fgut,gt,Id,Ipo,fgut,gp)
        
//        double EGP;
//        EGP=(Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)*Heavi3;

//        double Heavi1;
//        if( (gp-Kvect.at(32))>= 0){
//             Heavi1 = 1;}
//        if( (gp-Kvect.at(32)) < 0){
//             Heavi1 = 0;}

//        double Heavi3;
//        if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)  >= 0){
//             Heavi3 = 1;}
//        if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo) < 0){
//             Heavi3 = 0;}

//        K[0] = (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*fgut-Kvect.at(26)-Kvect.at(31)*(gp-Kvect.at(32))*Heavi1-Kvect.at(24)*gp+Kvect.at(25)*gt;

//        if( ((gp+ h*K[0]/2.0)-Kvect.at(32))>= 0){
//             Heavi1 = 1;}
//        if( ((gp+ h*K[0]/2.0)-Kvect.at(32)) < 0){
//             Heavi1 = 0;}
//        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[0]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0))  >= 0){
//             Heavi3 = 1;}
//        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[0]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0)) < 0){
//             Heavi3 = 0;}
//        K[1] = (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[0]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h/2.0)-Kvect.at(26)-Kvect.at(31)*((gp+ h*K[0]/2.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h*K[0]/2.0)+Kvect.at(25)*(gt+ h/2.0);

//        if( ((gp+ h*K[1]/2.0)-Kvect.at(32))>= 0){
//             Heavi1 = 1;}
//        if( ((gp+ h*K[1]/2.0)-Kvect.at(32)) < 0){
//             Heavi1 = 0;}
//        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[1]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0))  >= 0){
//             Heavi3 = 1;}
//        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[1]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0)) < 0){
//             Heavi3 = 0;}
//        K[2] = (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[1]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h/2.0)-Kvect.at(26)-Kvect.at(31)*((gp+ h*K[1]/2.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h*K[1]/2.0)+Kvect.at(25)*(gt+ h/2.0);

//        if( ((gp+ h*K[2])-Kvect.at(32))>= 0){
//             Heavi1 = 1;}
//        if( ((gp+ h*K[2])-Kvect.at(32)) < 0){
//             Heavi1 = 0;}
//        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[2])-Kvect.at(13)*(Id+ h)-Kvect.at(14)*(Ipo+ h))  >= 0){
//             Heavi3 = 1;}
//        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[2])-Kvect.at(13)*(Id+ h)-Kvect.at(14)*(Ipo+ h)) < 0){
//             Heavi3 = 0;}
//        K[3] = (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[2])-Kvect.at(13)*(Id+ h)-Kvect.at(14)*(Ipo+ h))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h)-Kvect.at(26)-Kvect.at(31)*((gp+ h*K[2])-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h*K[2])+Kvect.at(25)*(gt+ h);

//        dgp   = gp + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dgp через Дормана-Принса

        double Heavi1;
        double Heavi3;
        if( (gp-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( (gp-Kvect.at(32)) < 0){
             Heavi1 = 0;}

        if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo) < 0){
             Heavi3 = 0;}

        K[0] = h1*((Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*fgut-Kvect.at(26)-Kvect.at(31)*(gp-Kvect.at(32))*Heavi1-Kvect.at(24)*gp+Kvect.at(25)*gt);

        if( ((gp+ h1*K[0]/5.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ h1*K[0]/5.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h1*K[0]/5.0)-Kvect.at(13)*(Id+ h1/5.0)-Kvect.at(14)*(Ipo+ h1/5.0))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h1*K[0]/5.0)-Kvect.at(13)*(Id+ h1/5.0)-Kvect.at(14)*(Ipo+ h1/5.0)) < 0){
             Heavi3 = 0;}
        K[1] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ h1*K[0]/5.0)-Kvect.at(13)*(Id+ h1/5.0)-Kvect.at(14)*(Ipo+ h1/5.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1/5.0)-Kvect.at(26)-Kvect.at(31)*((gp+ h1*K[0]/5.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h1*K[0]/5.0)+Kvect.at(25)*(gt+ h1/5.0));

        if( ((gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(13)*(Id+ 3*h1/10.0)-Kvect.at(14)*(Ipo+ 3*h1/10.0))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(13)*(Id+ 3*h1/10.0)-Kvect.at(14)*(Ipo+ 3*h1/10.0)) < 0){
             Heavi3 = 0;}
        K[2] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(13)*(Id+ 3*h1/10.0)-Kvect.at(14)*(Ipo+ 3*h1/10.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ 3*h1/10.0)-Kvect.at(26)-Kvect.at(31)*((gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)+Kvect.at(25)*(gt+ 3*h1/10.0));

        if( ((gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(13)*(Id+ 4*h1/5.0)-Kvect.at(14)*(Ipo+ 4*h1/5.0))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(13)*(Id+ 4*h1/5.0)-Kvect.at(14)*(Ipo+ 4*h1/5.0)) < 0){
             Heavi3 = 0;}
        K[3] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(13)*(Id+ 4*h1/5.0)-Kvect.at(14)*(Ipo+ 4*h1/5.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ 4*h1/5.0)-Kvect.at(26)-Kvect.at(31)*((gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)+Kvect.at(25)*(gt+ 4*h1/5.0));
        
        if( ((gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(13)*(Id+ 8*h1/9.0)-Kvect.at(14)*(Ipo+ 8*h1/9.0))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(13)*(Id+ 8*h1/9.0)-Kvect.at(14)*(Ipo+ 8*h1/9.0)) < 0){
             Heavi3 = 0;}
        K[4] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(13)*(Id+ 8*h1/9.0)-Kvect.at(14)*(Ipo+ 8*h1/9.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ 8*h1/9.0)-Kvect.at(26)-Kvect.at(31)*((gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )+Kvect.at(25)*(gt+ 8*h1/9.0));

        if( ((gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1)) < 0){
             Heavi3 = 0;}
        K[5] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1)-Kvect.at(26)-Kvect.at(31)*((gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))+Kvect.at(25)*(gt+ h1));

        if( ((gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1)) < 0){
             Heavi3 = 0;}
        K[6] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1)-Kvect.at(26)-Kvect.at(31)*((gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)+Kvect.at(25)*(gt+ h1));

        dgp   = gp + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);     // во многих источниках домножается на h, но в одном нету такого.
        dz = gp + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dgp);
        s = pow(eps*h1/(2*err),1/5);
        hopt1 = s*h1;
        if( hopt1 < hmin) hopt1 = hmin;
        else if(hopt1 > hmax) hopt1 = hmax;

        // Решаем dgt через Рунге-Кутта

        // dgt=-((Kvect.at(27)+Kvect.at(28)*Xt)*gt)/(Kvect.at(29)+gt)+Kvect.at(24)*gp-Kvect.at(25)*gt;
        // dgt = f(Xt,gp,gt)

//        K[0] = -((Kvect.at(27)+Kvect.at(28)*Xt)*gt)/(Kvect.at(29)+gt)+Kvect.at(24)*gp-Kvect.at(25)*gt;
//        K[1] = -((Kvect.at(27)+Kvect.at(28)*(Xt + h/2.0))*(gt + h*K[0]/2.0))/(Kvect.at(29)+(gt + h*K[0]/2.0))+Kvect.at(24)*(gp + h/2.0)-Kvect.at(25)*(gt + h*K[0]/2.0);
//        K[2] = -((Kvect.at(27)+Kvect.at(28)*(Xt + h/2.0))*(gt + h*K[1]/2.0))/(Kvect.at(29)+(gt + h*K[1]/2.0))+Kvect.at(24)*(gp + h/2.0)-Kvect.at(25)*(gt + h*K[1]/2.0);
//        K[3] = -((Kvect.at(27)+Kvect.at(28)*(Xt + h))*(gt + h*K[2]))/(Kvect.at(29)+(gt + h*K[2]))+Kvect.at(24)*(gp + h)-Kvect.at(25)*(gt + h*K[2]);
//        dgt   = gt + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dgt через Дормана-Принса

        K[0] = h2*(-((Kvect.at(27)+Kvect.at(28)*Xt)*gt)/(Kvect.at(29)+gt)+Kvect.at(24)*gp-Kvect.at(25)*gt);
        K[1] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + h2/5.0))*(gt + h2*K[0]/5.0))/(Kvect.at(29)+(gt + h2*K[0]/5.0))+Kvect.at(24)*(gp + h2/5.0)-Kvect.at(25)*(gt + h2*K[0]/5.0));
        K[2] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + 3*h2/10.0))*(gt + 3*h2*K[0]/40.0 + 9*h2*K[1]/40.0))/(Kvect.at(29)+(gt + 3*h2*K[0]/40.0 + 9*h2*K[1]/40.0))+Kvect.at(24)*(gp + 3*h2/10.0)-Kvect.at(25)*(gt + 3*h2*K[0]/40.0 + 9*h2*K[1]/40.0));
        K[3] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + 4*h2/5.0))*(gt + 44*h2*K[0]/45.0 + (-56*h2*K[1]/15.0) + 32*h2*K[2]/9.0))/(Kvect.at(29)+(gt + 44*h2*K[0]/45.0 + (-56*h2*K[1]/15.0) + 32*h2*K[2]/9.0))+Kvect.at(24)*(gp + 4*h2/5.0)-Kvect.at(25)*(gt + 44*h2*K[0]/45.0 + (-56*h2*K[1]/15.0) + 32*h2*K[2]/9.0));
        K[4] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + 8*h2/9.0))*(gt + 19372*h2*K[0]/6561.0 + (-25360*h2*K[1]/2187.0) + 64448*h2*K[2]/6561.0 + (-212*h2*K[3]/729.0) ))/(Kvect.at(29)+(gt + 19372*h2*K[0]/6561.0 + (-25360*h2*K[1]/2187.0) + 64448*h2*K[2]/6561.0 + (-212*h2*K[3]/729.0) ))+Kvect.at(24)*(gp + 8*h2/9.0)-Kvect.at(25)*(gt + 19372*h2*K[0]/6561.0 + (-25360*h2*K[1]/2187.0) + 64448*h2*K[2]/6561.0 + (-212*h2*K[3]/729.0) ));
        K[5] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + h2))*(gt + 9017*h2*K[0]/3168.0 + (-355*h2*K[1]/33.0) + 46732*h2*K[2]/5247.0 + (49*h2*K[3]/176.0) + (-5103*h2*K[4]/18656.0) ))/(Kvect.at(29)+(gt + 9017*h2*K[0]/3168.0 + (-355*h2*K[1]/33.0) + 46732*h2*K[2]/5247.0 + (49*h2*K[3]/176.0) + (-5103*h2*K[4]/18656.0) ))+Kvect.at(24)*(gp + h2)-Kvect.at(25)*(gt + 9017*h2*K[0]/3168.0 + (-355*h2*K[1]/33.0) + 46732*h2*K[2]/5247.0 + (49*h2*K[3]/176.0) + (-5103*h2*K[4]/18656.0) ));
        K[6] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + h2))*(gt + 35*h2*K[0]/384.0 + 500*h2*K[2]/1113.0 + (125*h2*K[3]/192.0) + (-2187*h2*K[4]/6784.0) + 11*h2*K[5]/84.0))/(Kvect.at(29)+(gt + 35*h2*K[0]/384.0 + 500*h2*K[2]/1113.0 + (125*h2*K[3]/192.0) + (-2187*h2*K[4]/6784.0) + 11*h2*K[5]/84.0))+Kvect.at(24)*(gp + h2)-Kvect.at(25)*(gt + 35*h2*K[0]/384.0 + 500*h2*K[2]/1113.0 + (125*h2*K[3]/192.0) + (-2187*h2*K[4]/6784.0) + 11*h2*K[5]/84.0));

        dgt = gt + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);     // во многих источниках домножается на h, но в одном нету такого.
        dz = gt + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dgp);
        s = pow(eps*h2/(2*err),1/5);
        hopt2 = s*h2;
        if( hopt2 < hmin) hopt2 = hmin;
        else if(hopt2 > hmax) hopt2 = hmax;

        // Решаем dI1 через Рунге-Кутта

        // dI1=-Kvect.at(15)*(I1-Ip/Kvect.at(4));
        // dI1 = f(Ip,I1)

//        K[0] = -Kvect.at(15)*(I1-Ip/Kvect.at(4));
//        K[1] = -Kvect.at(15)*((I1 + h*K[0]/2.0)-(Ip+ h/2.0)/Kvect.at(4));
//        K[2] = -Kvect.at(15)*((I1 + h*K[1]/2.0)-(Ip+ h/2.0)/Kvect.at(4));
//        K[3] = -Kvect.at(15)*((I1 + h*K[2])-(Ip+ h)/Kvect.at(4));
//        dI1   = I1 + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dI1 через Дормана-Принса

        K[0] = h3*(-Kvect.at(15)*(I1-Ip/Kvect.at(4)));
        K[1] = h3*(-Kvect.at(15)*((I1 + h3*K[0]/5.0)-(Ip+ h3/5.0)/Kvect.at(4)));
        K[2] = h3*(-Kvect.at(15)*((I1 + 3*h3*K[0]/40.0 + 9*h3*K[1]/40.0)-(Ip+ 3*h3/10.0)/Kvect.at(4)));
        K[3] = h3*(-Kvect.at(15)*((I1 + 44*h3*K[0]/45.0 + (-56*h3*K[1]/15.0) + 32*h3*K[2]/9.0)-(Ip+ 4*h3/5.0)/Kvect.at(4)));
        K[4] = h3*(-Kvect.at(15)*((I1 + 19372*h3*K[0]/6561.0 + (-25360*h3*K[1]/2187.0) + 64448*h3*K[2]/6561.0 + (-212*h3*K[3]/729.0) )-(Ip+ 8*h3/9.0)/Kvect.at(4)));
        K[5] = h3*(-Kvect.at(15)*((I1 + 9017*h3*K[0]/3168.0 + (-355*h3*K[1]/33.0) + 46732*h3*K[2]/5247.0 + (49*h3*K[3]/176.0) + (-5103*h3*K[4]/18656.0) )-(Ip+ h3)/Kvect.at(4)));
        K[6] = h3*(-Kvect.at(15)*((I1 + 35*h3*K[0]/384.0 + 500*h3*K[2]/1113.0 + (125*h3*K[3]/192.0) + (-2187*h3*K[4]/6784.0) + 11*h3*K[5]/84.0 )-(Ip+ h3)/Kvect.at(4)));

        dI1 = I1 + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);     // во многих источниках домножается на h, но в одном нету такого.
        dz = I1 + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dI1);
        s = pow(eps*h3/(2*err),1/5);
        hopt3 = s*h3;
        if( hopt3 < hmin) hopt3 = hmin;
        else if(hopt3 > hmax) hopt3 = hmax;

        // Решаем dId через Рунге-Кутта
        // dId=-Kvect.at(15)*(Id-I1);
        // dId = f(I1,Id)?
//        K[0] = -Kvect.at(15)*(Id-I1);
//        K[1] = -Kvect.at(15)*((Id + h*K[0]/2.0)-(I1+h/2.0));
//        K[2] = -Kvect.at(15)*((Id + h*K[1]/2.0)-(I1+h/2.0));
//        K[3] = -Kvect.at(15)*((Id + h*K[2])-(I1+h));
//        dId   = Id + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dId через Дормана-Принса

        K[0] = h4*(-Kvect.at(15)*(Id-I1));
        K[1] =  h4*(-Kvect.at(15)*((Id + h4*K[0]/5.0)-(I1+h4/5.0)));
        K[2] =  h4*(-Kvect.at(15)*((Id + 3*h4*K[0]/40.0 + 9*h4*K[1]/40.0)-(I1+3*h4/10.0)));
        K[3] =  h4*(-Kvect.at(15)*((Id + 44*h4*K[0]/45.0 + (-56*h4*K[1]/15.0) + 32*h4*K[2]/9.0)-(I1+4*h4/5.0)));
        K[4] =  h4*(-Kvect.at(15)*((Id + 19372*h4*K[0]/6561.0 + (-25360*h4*K[1]/2187.0) + 64448*h4*K[2]/6561.0 + (-212*h4*K[3]/729.0))-(I1+8*h4/9.0)));
        K[5] =  h4*(-Kvect.at(15)*((Id + 9017*h4*K[0]/3168.0 + (-355*h4*K[1]/33.0) + 46732*h4*K[2]/5247.0 + (49*h4*K[3]/176.0) + (-5103*h4*K[4]/18656.0))-(I1+h4)));
        K[6] =  h4*(-Kvect.at(15)*((Id + 35*h4*K[0]/384.0 + 500*h4*K[2]/1113.0 + (125*h4*K[3]/192.0) + (-2187*h4*K[4]/6784.0) + 11*h4*K[5]/84.0)-(I1+h4)));

        dId   = Id + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);                                                // во многих источниках домножается на h, но в одном нету такого.
        dz = Id + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dId);
        s = pow(eps*h4/(2*err),1/5);
        hopt4 = s*h4;
        if( hopt4 < hmin) hopt4 = hmin;
        else if(hopt4 > hmax) hopt4 = hmax;


        // Решаем dXt через Рунге-Кутта
        // dXt=-Kvect.at(30)*Xt+Kvect.at(30)*((Ip/Kvect.at(4))-Kvect.at(1))*Heavi2;
        // dXt = f(Ip,Xt)?
//        double Heavi2;
//        if( ((Ip/Kvect.at(4))-Kvect.at(1)) >= 0){
//             Heavi2 = 1;}
//        if( ((Ip/Kvect.at(4))-Kvect.at(1)) < 0){
//             Heavi2 = 0;}
//        K[0] = -Kvect.at(30)*Xt+Kvect.at(30)*((Ip/Kvect.at(4))-Kvect.at(1))*Heavi2;

//        if( (((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
//             Heavi2 = 1;}
//        if( (((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1)) < 0){
//             Heavi2 = 0;}
//        K[1] = -Kvect.at(30)*(Xt+ h*K[0]/2.0)+Kvect.at(30)*(((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1))*Heavi2;

//        if( (((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
//             Heavi2 = 1;}
//        if( (((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1)) < 0){
//             Heavi2 = 0;}
//        K[2] = -Kvect.at(30)*(Xt + h*K[1]/2.0)+Kvect.at(30)*(((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1))*Heavi2;

//        if( (((Ip+ h)/Kvect.at(4))-Kvect.at(1)) >= 0){
//             Heavi2 = 1;}
//        if( (((Ip+ h)/Kvect.at(4))-Kvect.at(1)) < 0){
//             Heavi2 = 0;}
//        K[3] = -Kvect.at(30)*(Xt+ h*K[2])+Kvect.at(30)*(((Ip+ h)/Kvect.at(4))-Kvect.at(1))*Heavi2;

//        dXt   = Xt + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dXt через Дормана-Принса

        double Heavi2;
        if( ((Ip/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( ((Ip/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[0] = h5*(-Kvect.at(30)*Xt+Kvect.at(30)*((Ip/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((Ip+ h5/5.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((Ip+ h5/5.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[1] = h5*(-Kvect.at(30)*(Xt+ h5*K[0]/5.0)+Kvect.at(30)*(((Ip+ h5/5.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((Ip+ 3*h5/10.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((Ip+ 3*h5/10.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[2] = h5*(-Kvect.at(30)*(Xt+ 3*h5*K[0]/40.0 + 9*h5*K[1]/40.0)+Kvect.at(30)*(((Ip+ 3*h5/10.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((Ip+ 4*h5/5.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((Ip+ 4*h5/5.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[3] = h5*(-Kvect.at(30)*(Xt+ 44*h5*K[0]/45.0 + (-56*h5*K[1]/15.0) + 32*h5*K[2]/9.0)+Kvect.at(30)*(((Ip+ 4*h5/5.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((Ip+ 8*h5/9.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((Ip+ 8*h5/9.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[4] = h5*(-Kvect.at(30)*(Xt+ 19372*h5*K[0]/6561.0 + (-25360*h5*K[1]/2187.0) + 64448*h5*K[2]/6561.0 + (-212*h5*K[3]/729.0) )+Kvect.at(30)*(((Ip+ 8*h5/9.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((Ip+ h5)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((Ip+ h5)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[5] = h5*(-Kvect.at(30)*(Xt+ 9017*h5*K[0]/3168.0 + (-355*h5*K[1]/33.0) + 46732*h5*K[2]/5247.0 + (49*h5*K[3]/176.0) + (-5103*h5*K[4]/18656.0) )+Kvect.at(30)*(((Ip+ h5)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((Ip+ h5)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((Ip+ h5)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[6] = h5*(-Kvect.at(30)*(Xt+ 35*h5*K[0]/384.0 + 500*h5*K[2]/1113.0 + (125*h5*K[3]/192.0) + (-2187*h5*K[4]/6784.0) + 11*h5*K[5]/84.0)+Kvect.at(30)*(((Ip+ h5)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        dXt   = Xt + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);                                                // во многих источниках домножается на h, но в одном нету такого.
        dz = Xt + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dXt);
        s = pow(eps*h5/(2*err),1/5);
        hopt5 = s*h5;
        if( hopt5 < hmin) hopt5 = hmin;
        else if(hopt5 > hmax) hopt5 = hmax;


        // Решаем dIl через Рунге-Кутта
        // dIl=-(Kvect.at(5)+Kvect.at(7))*Il+Kvect.at(6)*Ip;
        // dIl = f(Ip,Il)?

//        K[0] = -(Kvect.at(5)+Kvect.at(7))*Il+Kvect.at(6)*Ip;
//        K[1] = -(Kvect.at(5)+Kvect.at(7))*(Il+ h*K[0]/2.0)+Kvect.at(6)*(Ip+ h/2.0);
//        K[2] = -(Kvect.at(5)+Kvect.at(7))*(Il+ h*K[1]/2.0)+Kvect.at(6)*(Ip+ h/2.0);
//        K[3] = -(Kvect.at(5)+Kvect.at(7))*(Il+ h*K[2])+Kvect.at(6)*(Ip+ h);

//        dIl   = Il + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIl через Дормана-Принса

        K[0] = h6*(-(Kvect.at(5)+Kvect.at(7))*Il+Kvect.at(6)*Ip);
        K[1] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ h6*K[0]/5.0)+Kvect.at(6)*(Ip+ h6/5.0));
        K[2] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 3*h6*K[0]/40.0 + 9*h6*K[1]/40.0)+Kvect.at(6)*(Ip+ 3*h6/10.0));
        K[3] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 44*h6*K[0]/45.0 + (-56*h6*K[1]/15.0) + 32*h6*K[2]/9.0)+Kvect.at(6)*(Ip+ 4*h6/5.0));
        K[4] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 19372*h6*K[0]/6561.0 + (-25360*h6*K[1]/2187.0) + 64448*h6*K[2]/6561.0 + (-212*h6*K[3]/729.0))+Kvect.at(6)*(Ip+ 8*h6/9.0));
        K[5] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 9017*h6*K[0]/3168.0 + (-355*h6*K[1]/33.0) + 46732*h6*K[2]/5247.0 + (49*h6*K[3]/176.0) + (-5103*h6*K[4]/18656.0))+Kvect.at(6)*(Ip+ h6));
        K[6] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 35*h6*K[0]/384.0 + 500*h6*K[2]/1113.0 + (125*h6*K[3]/192.0) + (-2187*h6*K[4]/6784.0) + 11*h6*K[5]/84.0)+Kvect.at(6)*(Ip+ h6));

        dIl   = Il + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);                                                // во многих источниках домножается на h, но в одном нету такого.
        dz = Il + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIl);
        s = pow(eps*h6/(2*err),1/5);
        hopt6 = s*h6;
        if( hopt6 < hmin) hopt6 = hmin;
        else if(hopt6 > hmax) hopt6 = hmax;

        // Решаем dIp через Рунге-Кутта
        // dIp=-Kvect.at(6)*Ip+Kvect.at(5)*Il+Kvect.at(10)/Kvect.at(0)*It-Kvect.at(9)*Ip;
        // dIp = f(Il,It,Ip)?
//        K[0] = -Kvect.at(6)*Ip+Kvect.at(5)*Il+Kvect.at(10)/Kvect.at(0)*It-Kvect.at(9)*Ip;
//        K[1] = -Kvect.at(6)*(Ip+ h*K[0]/2.0)+Kvect.at(5)*(Il+ h/2.0)+Kvect.at(10)/Kvect.at(0)*(It+ h/2.0)-Kvect.at(9)*(Ip+ h*K[0]/2.0);
//        K[2] = -Kvect.at(6)*(Ip+ h*K[1]/2.0)+Kvect.at(5)*(Il+ h/2.0)+Kvect.at(10)/Kvect.at(0)*(It+ h/2.0)-Kvect.at(9)*(Ip+ h*K[1]/2.0);
//        K[3] = -Kvect.at(6)*(Ip+ h*K[2])+Kvect.at(5)*(Il+ h)+Kvect.at(10)/Kvect.at(0)*(It+ h)-Kvect.at(9)*(Ip+ h*K[2]);

//        dIp   = Ip + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIp через Дормана-Принса

        K[0] = h7*(-Kvect.at(6)*Ip+Kvect.at(5)*Il+Kvect.at(10)/Kvect.at(0)*It-Kvect.at(9)*Ip);
        K[1] = h7*(-Kvect.at(6)*(Ip+ h7*K[0]/5.0)+Kvect.at(5)*(Il+ h7/5.0)+Kvect.at(10)/Kvect.at(0)*(It+ h7/5.0)-Kvect.at(9)*(Ip+ h7*K[0]/5.0));
        K[2] = h7*(-Kvect.at(6)*(Ip+ 3*h7*K[0]/40.0 + 9*h7*K[1]/40.0)+Kvect.at(5)*(Il+ 3*h7/10.0)+Kvect.at(10)/Kvect.at(0)*(It+ 3*h7/10.0)-Kvect.at(9)*(Ip+ 3*h7*K[0]/40.0 + 9*h7*K[1]/40.0));
        K[3] = h7*(-Kvect.at(6)*(Ip+ 44*h7*K[0]/45.0 + (-56*h7*K[1]/15.0) + 32*h7*K[2]/9.0)+Kvect.at(5)*(Il+ 4*h7/5.0)+Kvect.at(10)/Kvect.at(0)*(It+ 4*h7/5.0)-Kvect.at(9)*(Ip+ 44*h7*K[0]/45.0 + (-56*h7*K[1]/15.0) + 32*h7*K[2]/9.0));
        K[4] = h7*(-Kvect.at(6)*(Ip+ 19372*h7*K[0]/6561.0 + (-25360*h7*K[1]/2187.0) + 64448*h7*K[2]/6561.0 + (-212*h7*K[3]/729.0) )+Kvect.at(5)*(Il+ 8*h7/9.0)+Kvect.at(10)/Kvect.at(0)*(It+ 8*h7/9.0)-Kvect.at(9)*(Ip+ 19372*h7*K[0]/6561.0 + (-25360*h7*K[1]/2187.0) + 64448*h7*K[2]/6561.0 + (-212*h7*K[3]/729.0) ));
        K[5] = h7*(-Kvect.at(6)*(Ip+ 9017*h7*K[0]/3168.0 + (-355*h7*K[1]/33.0) + 46732*h7*K[2]/5247.0 + (49*h7*K[3]/176.0) + (-5103*h7*K[4]/18656.0))+Kvect.at(5)*(Il+ h7)+Kvect.at(10)/Kvect.at(0)*(It+ h7)-Kvect.at(9)*(Ip+ 9017*h7*K[0]/3168.0 + (-355*h7*K[1]/33.0) + 46732*h7*K[2]/5247.0 + (49*h7*K[3]/176.0) + (-5103*h7*K[4]/18656.0)));
        K[6] = h7*(-Kvect.at(6)*(Ip+ 35*h7*K[0]/384.0 + 500*h7*K[2]/1113.0 + (125*h7*K[3]/192.0) + (-2187*h7*K[4]/6784.0) + 11*h7*K[5]/84.0)+Kvect.at(5)*(Il+ h7)+Kvect.at(10)/Kvect.at(0)*(It+ h7)-Kvect.at(9)*(Ip+ 35*h7*K[0]/384.0 + 500*h7*K[2]/1113.0 + (125*h7*K[3]/192.0) + (-2187*h7*K[4]/6784.0) + 11*h7*K[5]/84.0));

        dIp   = Ip + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = Ip + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIp);
        s = pow(eps*h7/(2*err),1/5);
        hopt7 = s*h7;
        if( hopt7 < hmin) hopt7 = hmin;
        else if(hopt7 > hmax) hopt7 = hmax;

        // Решаем dfgut через Рунге-Кутта
        // dfgut=-Kvect.at(17)*fgut+kgut*fliq;
        // dfgut = f(fliq,fsol,fgut)?

//        double kgut;
//        kgut=Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2);

//        K[0] = -Kvect.at(17)*fgut+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2)*fliq;
//        K[1] = -Kvect.at(17)*(fgut+ h*K[0]/2.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h/2.0)+(fliq+ h/2.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h/2.0)+(fliq+ h/2.0)-Kvect.at(22)*Dig))+2)*(fliq+ h/2.0);
//        K[2] = -Kvect.at(17)*(fgut+ h*K[1]/2.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h/2.0)+(fliq+ h/2.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h/2.0)+(fliq+ h/2.0)-Kvect.at(22)*Dig))+2)*(fliq+ h/2.0);
//        K[3] = -Kvect.at(17)*(fgut+ h*K[2])+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h)+(fliq+ h)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h)+(fliq+ h)-Kvect.at(22)*Dig))+2)*(fliq+ h);

//        dfgut   = fgut + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfgut через Дормана-Принса

        K[0] = h8*(-Kvect.at(17)*fgut+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2)*fliq);
        K[1] = h8*(-Kvect.at(17)*(fgut+ h8*K[0]/5.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h8/5.0)+(fliq+ h8/5.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h8/5.0)+(fliq+ h8/5.0)-Kvect.at(22)*Dig))+2)*(fliq+ h8/5.0));
        K[2] = h8*(-Kvect.at(17)*(fgut+ 3*h8*K[0]/40.0 + 9*h8*K[1]/40.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 3*h8/10.0)+(fliq+ 3*h8/10.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 3*h8/10.0)+(fliq+ 3*h8/10.0)-Kvect.at(22)*Dig))+2)*(fliq+ 3*h8/10.0));
        K[3] = h8*(-Kvect.at(17)*(fgut+ 44*h8*K[0]/45.0 + (-56*h8*K[1]/15.0) + 32*h8*K[2]/9.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 4*h8/5.0)+(fliq+ 4*h8/5.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 4*h8/5.0)+(fliq+ 4*h8/5.0)-Kvect.at(22)*Dig))+2)*(fliq+ 4*h8/5.0));
        K[4] = h8*(-Kvect.at(17)*(fgut+ 19372*h8*K[0]/6561.0 + (-25360*h8*K[1]/2187.0) + 64448*h8*K[2]/6561.0 + (-212*h8*K[3]/729.0) )+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 8*h8/9.0)+(fliq+ 8*h8/9.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 8*h8/9.0)+(fliq+ 8*h8/9.0)-Kvect.at(22)*Dig))+2)*(fliq+ 8*h8/9.0));
        K[5] = h8*(-Kvect.at(17)*(fgut+ 9017*h8*K[0]/3168.0 + (-355*h8*K[1]/33.0) + 46732*h8*K[2]/5247.0 + (49*h8*K[3]/176.0) + (-5103*h8*K[4]/18656.0))+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h8)+(fliq+ h8)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h8)+(fliq+ h8)-Kvect.at(22)*Dig))+2)*(fliq+ h8));
        K[6] = h8*(-Kvect.at(17)*(fgut+ 35*h8*K[0]/384.0 + 500*h8*K[2]/1113.0 + (125*h8*K[3]/192.0) + (-2187*h8*K[4]/6784.0) + 11*h8*K[5]/84.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h8)+(fliq+ h8)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h8)+(fliq+ h8)-Kvect.at(22)*Dig))+2)*(fliq+ h8));

        dfgut   = fgut + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = fgut + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dfgut);
        s = pow(eps*h8/(2*err),1/5);
        hopt8 = s*h8;
        if( hopt8 < hmin) hopt8 = hmin;
        else if(hopt8 > hmax) hopt8 = hmax;

        // Решаем dfliq через Рунге-Кутта
        // dfliq=-kgut*fliq+Kvect.at(18)*fsol;
        // dfliq = f(fliq,fsol)?

//        K[0] = -(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2))*fliq+Kvect.at(18)*fsol;
//        K[1] = -(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h/2.0)+(fliq+ h*K[0]/2.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h/2.0)+(fliq+ h*K[0]/2.0)-Kvect.at(22)*Dig))+2))*(fliq+ h*K[0]/2.0)+Kvect.at(18)*(fsol+ h/2.0);
//        K[2] = -(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h/2.0)+(fliq+ h*K[1]/2.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h/2.0)+(fliq+ h*K[1]/2.0)-Kvect.at(22)*Dig))+2))*(fliq+ h*K[1]/2.0)+Kvect.at(18)*(fsol+ h/2.0);
//        K[2] = -(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h)+(fliq+ h*K[2])-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h)+(fliq+ h*K[2])-Kvect.at(22)*Dig))+2))*(fliq+ h*K[2])+Kvect.at(18)*(fsol+ h);

//        dfliq   = fliq + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfliq через Дормана-Принса

        K[0] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2))*fliq+Kvect.at(18)*fsol);
        K[1] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h9/5.0)+(fliq+ h9*K[0]/5.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h9/5.0)+(fliq+ h9*K[0]/5.0)-Kvect.at(22)*Dig))+2))*(fliq+ h9*K[0]/5.0)+Kvect.at(18)*(fsol+ h9/5.0));
        K[2] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 3*h9/10.0)+(fliq+ 3*h9*K[0]/40.0 + 9*h9*K[1]/40.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 3*h9/10.0)+(fliq+3*h9*K[0]/40.0 + 9*h9*K[1]/40.0)-Kvect.at(22)*Dig))+2))*(fliq+ 3*h9*K[0]/40.0 + 9*h9*K[1]/40.0)+Kvect.at(18)*(fsol+ 3*h9/10.0));
        K[3] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 4*h9/5.0)+(fliq+ 44*h9*K[0]/45.0 + (-56*h9*K[1]/15.0) + 32*h9*K[2]/9.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 4*h9/5.0)+(fliq+ 44*h9*K[0]/45.0 + (-56*h9*K[1]/15.0) + 32*h9*K[2]/9.0)-Kvect.at(22)*Dig))+2))*(fliq+ 44*h9*K[0]/45.0 + (-56*h9*K[1]/15.0) + 32*h9*K[2]/9.0)+Kvect.at(18)*(fsol+ 4*h9/5.0));
        K[4] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 8*h9/9.0)+(fliq+ 19372*h9*K[0]/6561.0 + (-25360*h9*K[1]/2187.0) + 64448*h9*K[2]/6561.0 + (-212*h9*K[3]/729.0))-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 8*h9/9.0)+(fliq+ 19372*h9*K[0]/6561.0 + (-25360*h9*K[1]/2187.0) + 64448*h9*K[2]/6561.0 + (-212*h9*K[3]/729.0))-Kvect.at(22)*Dig))+2))*(fliq+ 19372*h9*K[0]/6561.0 + (-25360*h9*K[1]/2187.0) + 64448*h9*K[2]/6561.0 + (-212*h9*K[3]/729.0))+Kvect.at(18)*(fsol+ 8*h9/9.0));
        K[5] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h9)+(fliq+ 9017*h9*K[0]/3168.0 + (-355*h9*K[1]/33.0) + 46732*h9*K[2]/5247.0 + (49*h9*K[3]/176.0) + (-5103*h9*K[4]/18656.0))-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h9)+(fliq+ 9017*h9*K[0]/3168.0 + (-355*h9*K[1]/33.0) + 46732*h9*K[2]/5247.0 + (49*h9*K[3]/176.0) + (-5103*h9*K[4]/18656.0))-Kvect.at(22)*Dig))+2))*(fliq+ 9017*h9*K[0]/3168.0 + (-355*h9*K[1]/33.0) + 46732*h9*K[2]/5247.0 + (49*h9*K[3]/176.0) + (-5103*h9*K[4]/18656.0))+Kvect.at(18)*(fsol+ h9));
        K[6] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h9)+(fliq+ 35*h9*K[0]/384.0 + 500*h9*K[2]/1113.0 + (125*h9*K[3]/192.0) + (-2187*h9*K[4]/6784.0) + 11*h9*K[5]/84.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h9)+(fliq+ 35*h9*K[0]/384.0 + 500*h9*K[2]/1113.0 + (125*h9*K[3]/192.0) + (-2187*h9*K[4]/6784.0) + 11*h9*K[5]/84.0)-Kvect.at(22)*Dig))+2))*(fliq+ 35*h9*K[0]/384.0 + 500*h9*K[2]/1113.0 + (125*h9*K[3]/192.0) + (-2187*h9*K[4]/6784.0) + 11*h9*K[5]/84.0)+Kvect.at(18)*(fsol+ h9));

        dfliq   = fliq + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = fliq + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dfliq);
        s = pow(eps*h9/(2*err),1/5);
        hopt9 = s*h9;
        if( hopt9 < hmin) hopt9 = hmin;
        else if(hopt9 > hmax) hopt9 = hmax;

        // Решаем dfsol через Рунге-Кутта
        // dfsol=-Kvect.at(18)*fsol+vm;
        // dfsol = f(fsol,vm)?

//        K[0] = -Kvect.at(18)*fsol+vm;
//        K[1] = -Kvect.at(18)*(fsol+ h*K[0]/2.0)+(vm+ h/2.0);
//        K[2] = -Kvect.at(18)*(fsol+ h*K[1]/2.0)+(vm+ h/2.0);
//        K[3] = -Kvect.at(18)*(fsol+ h*K[2])+(vm+ h);

//        dfsol   = fsol + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfsol через Дормана-Принса

        K[0] = h10*(-Kvect.at(18)*fsol+vm);
        K[1] = h10*(-Kvect.at(18)*(fsol+ h10*K[0]/5.0)+(vm+ h10/5.0));
        K[2] = h10*(-Kvect.at(18)*(fsol+ 3*h10*K[0]/40.0 + 9*h10*K[1]/40.0)+(vm+ 3*h10/10.0));
        K[3] = h10*(-Kvect.at(18)*(fsol+ 44*h10*K[0]/45.0 + (-56*h10*K[1]/15.0) + 32*h10*K[2]/9.0)+(vm+ 4*h10/5.0));
        K[4] = h10*(-Kvect.at(18)*(fsol+ 19372*h10*K[0]/6561.0 + (-25360*h10*K[1]/2187.0) + 64448*h10*K[2]/6561.0 + (-212*h10*K[3]/729.0))+(vm+ 8*h10/9.0));
        K[5] = h10*(-Kvect.at(18)*(fsol+ 9017*h10*K[0]/3168.0 + (-355*h10*K[1]/33.0) + 46732*h10*K[2]/5247.0 + (49*h10*K[3]/176.0) + (-5103*h10*K[4]/18656.0))+(vm+ h10));
        K[6] = h10*(-Kvect.at(18)*(fsol+ 35*h10*K[0]/384.0 + 500*h10*K[2]/1113.0 + (125*h10*K[3]/192.0) + (-2187*h10*K[4]/6784.0) + 11*h10*K[5]/84.0)+(vm+ h10));

        dfsol   = fsol + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = fsol + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dfsol);
        s = pow(eps*h10/(2*err),1/5);
        hopt10 = s*h10;
        if( hopt10 < hmin) hopt10 = hmin;
        else if(hopt10 > hmax) hopt10 = hmax;

        // Решаем dIpo через Рунге-Кутта
        // dIpo=-Kvect.at(33)*Ipo+(Yt+Kvect.at(3))*Heavi5+(Yt+Kvect.at(3))*(Heavi6);
        // dIpo = f(Yt,dgp,Ipo)?
//        double Heavi5;
//        double Heavi6;

//        if( (dgp) >= 0){
//             Heavi5 = 1;}
//        if( (dgp) < 0){
//             Heavi5 = 0;}

//        if( (-dgp) >= 0){
//             Heavi6 = 1;}
//        if( (-dgp) < 0){
//             Heavi6 = 0;}
//        K[0] = -Kvect.at(33)*Ipo+(Yt+Kvect.at(3))*Heavi5+(Yt+Kvect.at(3))*(Heavi6);

//        if( ((dgp+ h/2.0)) >= 0){
//             Heavi5 = 1;}
//        if( ((dgp+ h/2.0)) < 0){
//             Heavi5 = 0;}

//        if( (-(dgp+ h/2.0)) >= 0){
//             Heavi6 = 1;}
//        if( (-(dgp+ h/2.0)) < 0){
//             Heavi6 = 0;}
//        K[1] = -Kvect.at(33)*(Ipo+ h*K[0]/2.0)+((Yt+ h/2.0)+Kvect.at(3))*Heavi5+((Yt+ h/2.0)+Kvect.at(3))*(Heavi6);

//        if( ((dgp+ h/2.0)) >= 0){
//             Heavi5 = 1;}
//        if( ((dgp+ h/2.0)) < 0){
//             Heavi5 = 0;}

//        if( (-(dgp+ h/2.0)) >= 0){
//             Heavi6 = 1;}
//        if( (-(dgp+ h/2.0)) < 0){
//             Heavi6 = 0;}
//        K[2] = -Kvect.at(33)*(Ipo+ h*K[1]/2.0)+((Yt+ h/2.0)+Kvect.at(3))*Heavi5+((Yt+ h/2.0)+Kvect.at(3))*(Heavi6);

//        if( ((dgp+ h)) >= 0){
//             Heavi5 = 1;}
//        if( ((dgp+ h)) < 0){
//             Heavi5 = 0;}

//        if( (-(dgp+ h)) >= 0){
//             Heavi6 = 1;}
//        if( (-(dgp+ h)) < 0){
//             Heavi6 = 0;}
//        K[3] = -Kvect.at(33)*(Ipo+ h*K[2])+((Yt+ h)+Kvect.at(3))*Heavi5+((Yt+ h)+Kvect.at(3))*(Heavi6);

//        dIpo   = Ipo + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIpo через Дормана-Принса

        double Heavi5;
        double Heavi6;

        if( (dgp) >= 0){
             Heavi5 = 1;}
        if( (dgp) < 0){
             Heavi5 = 0;}

        if( (-dgp) >= 0){
             Heavi6 = 1;}
        if( (-dgp) < 0){
             Heavi6 = 0;}
        K[0] = h11*(-Kvect.at(33)*Ipo+(Yt+Kvect.at(3))*Heavi5+(Yt+Kvect.at(3))*(Heavi6));

        if( ((dgp+ h11/5.0)) >= 0){
             Heavi5 = 1;}
        if( ((dgp+ h11/5.0)) < 0){
             Heavi5 = 0;}

        if( (-(dgp+ h11/5.0)) >= 0){
             Heavi6 = 1;}
        if( (-(dgp+ h11/5.0)) < 0){
             Heavi6 = 0;}
        K[1] = h11*(-Kvect.at(33)*(Ipo+ h11*K[0]/5.0)+((Yt+ h11/5.0)+Kvect.at(3))*Heavi5+((Yt+ h11/5.0)+Kvect.at(3))*(Heavi6));

        if( ((dgp+ 3*h11/10.0)) >= 0){
             Heavi5 = 1;}
        if( ((dgp+ 3*h11/10.0)) < 0){
             Heavi5 = 0;}

        if( (-(dgp+ 3*h11/10.0)) >= 0){
             Heavi6 = 1;}
        if( (-(dgp+ 3*h11/10.0)) < 0){
             Heavi6 = 0;}
        K[2] = h11*(-Kvect.at(33)*(Ipo+ 3*h11*K[0]/40.0 + 9*h11*K[1]/40.0)+((Yt+ 3*h11/10.0)+Kvect.at(3))*Heavi5+((Yt+ 3*h11/10.0)+Kvect.at(3))*(Heavi6));

        if( ((dgp+ 4*h11/5.0)) >= 0){
             Heavi5 = 1;}
        if( ((dgp+ 4*h11/5.0)) < 0){
             Heavi5 = 0;}

        if( (-(dgp+ 4*h11/5.0)) >= 0){
             Heavi6 = 1;}
        if( (-(dgp+ 4*h11/5.0)) < 0){
             Heavi6 = 0;}
        K[3] = h11*(-Kvect.at(33)*(Ipo+ 44*h11*K[0]/45.0 + (-56*h11*K[1]/15.0) + 32*h11*K[2]/9.0)+((Yt+ 4*h11/5.0)+Kvect.at(3))*Heavi5+((Yt+ 4*h11/5.0)+Kvect.at(3))*(Heavi6));

        if( ((dgp+ 4*h11/5.0)) >= 0){
             Heavi5 = 1;}
        if( ((dgp+ 4*h11/5.0)) < 0){
             Heavi5 = 0;}

        if( (-(dgp+ 8*h11/9.0)) >= 0){
             Heavi6 = 1;}
        if( (-(dgp+ 8*h11/9.0)) < 0){
             Heavi6 = 0;}
        K[4] = h11*(-Kvect.at(33)*(Ipo+ 19372*h11*K[0]/6561.0 + (-25360*h11*K[1]/2187.0) + 64448*h11*K[2]/6561.0 + (-212*h11*K[3]/729.0))+((Yt+ 8*h11/9.0)+Kvect.at(3))*Heavi5+((Yt+ 8*h11/9.0)+Kvect.at(3))*(Heavi6));

        if( ((dgp+ h11)) >= 0){
             Heavi5 = 1;}
        if( ((dgp+ h11)) < 0){
             Heavi5 = 0;}

        if( (-(dgp+ h11)) >= 0){
             Heavi6 = 1;}
        if( (-(dgp+ h11)) < 0){
             Heavi6 = 0;}
        K[5] = h11*(-Kvect.at(33)*(Ipo+ 9017*h11*K[0]/3168.0 + (-355*h11*K[1]/33.0) + 46732*h11*K[2]/5247.0 + (49*h11*K[3]/176.0) + (-5103*h11*K[4]/18656.0))+((Yt+ h11)+Kvect.at(3))*Heavi5+((Yt+ h11)+Kvect.at(3))*(Heavi6));

        if( ((dgp+ h11)) >= 0){
             Heavi5 = 1;}
        if( ((dgp+ h11)) < 0){
             Heavi5 = 0;}

        if( (-(dgp+ h11)) >= 0){
             Heavi6 = 1;}
        if( (-(dgp+ h11)) < 0){
             Heavi6 = 0;}
        K[6] = h11*(-Kvect.at(33)*(Ipo+ 35*h11*K[0]/384.0 + 500*h11*K[2]/1113.0 + (125*h11*K[3]/192.0) + (-2187*h11*K[4]/6784.0) + 11*h11*K[5]/84.0)+((Yt+ h11)+Kvect.at(3))*Heavi5+((Yt+ h11)+Kvect.at(3))*(Heavi6));

        dIpo   = Ipo + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = Ipo + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIpo);
        s = pow(eps*h11/(2*err),1/5);
        hopt11 = s*h11;
        if( hopt11 < hmin) hopt11 = hmin;
        else if(hopt11 > hmax) hopt11 = hmax;


        // Решаем dYt через Рунге-Кутта
        // dYt=-Kvect.at(34)*(Yt-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*Yt-Kvect.at(34)*Kvect.at(3))*(Heavi8);
        // dYt = f(gp,Yt)?
//        double Heavi7;
//        double Heavi8;

//        if( (Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//             Heavi7 = 1;}
//        if( (Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//             Heavi7 = 0;}

//        if( (-Kvect.at(3)-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))) >= 0){
//             Heavi8 = 1;}
//        if( (-Kvect.at(3)-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))) < 0){
//             Heavi8 = 0;}
//        K[0] = -Kvect.at(34)*(Yt-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*Yt-Kvect.at(34)*Kvect.at(3))*(Heavi8);

//        if( (Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//             Heavi7 = 1;}
//        if( (Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//             Heavi7 = 0;}

//        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
//             Heavi8 = 1;}
//        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))) < 0){
//             Heavi8 = 0;}
//        K[1] = -Kvect.at(34)*((Yt+ h*K[0]/2.0)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ h*K[0]/2.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8);

//        if( (Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//             Heavi7 = 1;}
//        if( (Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//             Heavi7 = 0;}

//        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
//             Heavi8 = 1;}
//        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))) < 0){
//             Heavi8 = 0;}
//        K[2] = -Kvect.at(34)*((Yt+ h*K[1]/2.0)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ h*K[1]/2.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8);

//        if( (Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//             Heavi7 = 1;}
//        if( (Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//             Heavi7 = 0;}

//        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2))) >= 0){
//             Heavi8 = 1;}
//        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2))) < 0){
//             Heavi8 = 0;}
//        K[3] = -Kvect.at(34)*((Yt+ h*K[2])-Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ h*K[2])-Kvect.at(34)*Kvect.at(3))*(Heavi8);

//        dYt   = Yt + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dYt через Дормана-Принса

        double Heavi7;
        double Heavi8;

        if( (Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[0] = h12*(-Kvect.at(34)*(Yt-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*Yt-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[1] = h12*(-Kvect.at(34)*((Yt+ h12*K[0]/5.0)-Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ h12*K[0]/5.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[2] = h12*(-Kvect.at(34)*((Yt+ 3*h12*K[0]/40.0 + 9*h12*K[1]/40.0)-Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 3*h12*K[0]/40.0 + 9*h12*K[1]/40.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[3] = h12*(-Kvect.at(34)*((Yt+ 44*h12*K[0]/45.0 + (-56*h12*K[1]/15.0) + 32*h12*K[2]/9.0)-Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 44*h12*K[0]/45.0 + (-56*h12*K[1]/15.0) + 32*h12*K[2]/9.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[4] = h12*(-Kvect.at(34)*((Yt+ 19372*h12*K[0]/6561.0 + (-25360*h12*K[1]/2187.0) + 64448*h12*K[2]/6561.0 + (-212*h12*K[3]/729.0))-Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 19372*h12*K[0]/6561.0 + (-25360*h12*K[1]/2187.0) + 64448*h12*K[2]/6561.0 + (-212*h12*K[3]/729.0))-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[5] = h12*(-Kvect.at(34)*((Yt+ 9017*h12*K[0]/3168.0 + (-355*h12*K[1]/33.0) + 46732*h12*K[2]/5247.0 + (49*h12*K[3]/176.0) + (-5103*h12*K[4]/18656.0))-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 9017*h12*K[0]/3168.0 + (-355*h12*K[1]/33.0) + 46732*h12*K[2]/5247.0 + (49*h12*K[3]/176.0) + (-5103*h12*K[4]/18656.0))-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[6] = h12*(-Kvect.at(34)*((Yt+ 35*h12*K[0]/384.0 + 500*h12*K[2]/1113.0 + (125*h12*K[3]/192.0) + (-2187*h12*K[4]/6784.0) + 11*h12*K[5]/84.0)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 35*h12*K[0]/384.0 + 500*h12*K[2]/1113.0 + (125*h12*K[3]/192.0) + (-2187*h12*K[4]/6784.0) + 11*h12*K[5]/84.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        dYt   = Yt + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = Yt + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dYt);
        s = pow(eps*h12/(2*err),1/5);
        hopt12 = s*h12;
        if( hopt12 < hmin) hopt12 = hmin;
        else if(hopt12 > hmax) hopt12 = hmax;

        // Решаем dIt через Рунге-Кутта
        // dIt=Kvect.at(8)*Ii-Kvect.at(10)*It;
        // dIt = f(Ii,It)?
//        K[0] = Kvect.at(8)*Ii-Kvect.at(10)*It;
//        K[1] = Kvect.at(8)*(Ii+ h/2.0)-Kvect.at(10)*(It+ h*K[0]/2.0);
//        K[2] = Kvect.at(8)*(Ii+ h/2.0)-Kvect.at(10)*(It+ h*K[1]/2.0);
//        K[3] = Kvect.at(8)*(Ii+ h)-Kvect.at(10)*(It+ h*K[2]);

//        dIt   = It + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIt через Дормана-Принса

        K[0] = h13*(Kvect.at(8)*Ii-Kvect.at(10)*It);
        K[1] = h13*(Kvect.at(8)*(Ii+ h13/5.0)-Kvect.at(10)*(It+ h13*K[0]/5.0));
        K[2] = h13*(Kvect.at(8)*(Ii+ 3*h13/10.0)-Kvect.at(10)*(It+ 3*h13*K[0]/40.0 + 9*h13*K[1]/40.0));
        K[3] = h13*(Kvect.at(8)*(Ii+ 4*h13/5.0)-Kvect.at(10)*(It+ 44*h13*K[0]/45.0 + (-56*h13*K[1]/15.0) + 32*h13*K[2]/9.0));
        K[4] = h13*(Kvect.at(8)*(Ii+ 8*h13/9.0)-Kvect.at(10)*(It+ 19372*h13*K[0]/6561.0 + (-25360*h13*K[1]/2187.0) + 64448*h13*K[2]/6561.0 + (-212*h13*K[3]/729.0)));
        K[5] = h13*(Kvect.at(8)*(Ii+ h13)-Kvect.at(10)*(It+ 9017*h13*K[0]/3168.0 + (-355*h13*K[1]/33.0) + 46732*h13*K[2]/5247.0 + (49*h13*K[3]/176.0) + (-5103*h13*K[4]/18656.0)));
        K[6] = h13*(Kvect.at(8)*(Ii+ h13)-Kvect.at(10)*(It+ 35*h13*K[0]/384.0 + 500*h13*K[2]/1113.0 + (125*h13*K[3]/192.0) + (-2187*h13*K[4]/6784.0) + 11*h13*K[5]/84.0));

        dIt   = It + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = It + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIt);
        s = pow(eps*h13/(2*err),1/5);
        hopt13 = s*h13;
        if( hopt13 < hmin) hopt13 = hmin;
        else if(hopt13 > hmax) hopt13 = hmax;

        // Решаем dIi через Рунге-Кутта
        // dIi=-Kvect.at(8)*Ii+vbas+vbol;
        // dIi = f(vbol,Ii,vbas)?
//        K[0] = -Kvect.at(8)*Ii+vbas+vbol;
//        K[1] = -Kvect.at(8)*(Ii+ h*K[0]/2.0)+(vbas+ h/2.0)+(vbol+ h/2.0);
//        K[2] = -Kvect.at(8)*(Ii+ h*K[1]/2.0)+(vbas+ h/2.0)+(vbol+ h/2.0);
//        K[3] = -Kvect.at(8)*(Ii+ h*K[2])+(vbas+ h)+(vbol+ h);

//        dIi   = Ii + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIi через Дормана-Принса

        K[0] = h14*(-Kvect.at(8)*Ii+vbas+vbol);
        K[1] = h14*(-Kvect.at(8)*(Ii+ h14*K[0]/5.0)+vbas+(vbol+ h14/5.0));
        K[2] = h14*(-Kvect.at(8)*(Ii+ 3*h14*K[0]/40.0 + 9*h14*K[1]/40.0)+(vbas+ 3*h14/10.0)+(vbol+ 3*h14/10.0));
        K[3] = h14*(-Kvect.at(8)*(Ii+ 44*h14*K[0]/45.0 + (-56*h14*K[1]/15.0) + 32*h14*K[2]/9.0)+(vbas+ 4*h14/5.0)+(vbol+ 4*h14/5.0));
        K[4] = h14*(-Kvect.at(8)*(Ii+ 19372*h14*K[0]/6561.0 + (-25360*h14*K[1]/2187.0) + 64448*h14*K[2]/6561.0 + (-212*h14*K[3]/729.0))+(vbas+ 8*h14/9.0)+(vbol+ 8*h14/9.0));
        K[5] = h14*(-Kvect.at(8)*(Ii+ 9017*h14*K[0]/3168.0 + (-355*h14*K[1]/33.0) + 46732*h14*K[2]/5247.0 + (49*h14*K[3]/176.0) + (-5103*h14*K[4]/18656.0))+(vbas+ h14)+(vbol+ h14));
        K[6] = h14*(-Kvect.at(8)*(Ii+ 35*h14*K[0]/384.0 + 500*h14*K[2]/1113.0 + (125*h14*K[3]/192.0) + (-2187*h14*K[4]/6784.0) + 11*h14*K[5]/84.0)+(vbas+ h14)+(vbol+ h14));

        dIi   = Ii + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = Ii + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIi);
        s = pow(eps*h14/(2*err),1/5);
        hopt14 = s*h14;
        if( hopt14 < hmin) hopt14 = hmin;
        else if(hopt14 > hmax) hopt14 = hmax;

        // Расчёты Vbas

        bas[5] = bas[0]*(1+bas[4]);
        double t1 = 30;
        double t2 = 30 +120;
        if ( (t<t1) or (t>t2) ){
            double v;
            double Heavi10;
            double Heavi11;
            double Heavi12;
            double Heavi13;
            double Heavi14;

            if( (gp-bas[1]) >= 0){
                 Heavi10 = 1;}
            if( (gp-bas[1]) < 0){
                 Heavi10 = 0;}
            if( (bas[6]-gp) >= 0){
                 Heavi11 = 1;}
            if( (bas[6]-gp) < 0){
                 Heavi11 = 0;}
            if( ((bas[2]-gp)) >= 0){
                 Heavi12 = 1;}
            if( ((bas[2]-gp)) < 0){
                 Heavi12 = 0;}
            if( (gp-bas[3]) >= 0){
                 Heavi13 = 1;}
            if( (gp-bas[3]) < 0){
                 Heavi13 = 0;}
            if( (gp-bas[6]) >= 0){
                 Heavi14 = 1;}
            if( (gp-bas[6]) < 0){
                 Heavi14 = 0;}
            v = (bas[5]-bas[0])/(bas[6]-bas[1])*(gp-bas[1])*Heavi10*Heavi11-bas[0]/(bas[2]-bas[3])*(bas[2]-gp)*Heavi12*Heavi13+bas[0]*Heavi13+(bas[5]-bas[0])*Heavi14;
            Vbas = v / Kvect.at(23);
        }
        else{

            Vbas = Vb;
        }

        grafVbas.append(Vbas);

        // присвоение конечных значений в начальные
        DataZ currentValues;
        currentValues.I1 = dI1;
        currentValues.Id = dId;
        currentValues.Ii = dIi;
        currentValues.Il = dIl;
        currentValues.Ip = dIp;
        currentValues.It = dIt;
        currentValues.Xt = dXt;
        currentValues.Yt = dYt;
        currentValues.gp = dgp;
        currentValues.gt = dgt;
        currentValues.Ipo = dIpo;
        currentValues.fgut = dfgut;
        currentValues.fliq = dfliq;
        currentValues.fsol = dfsol;
        currentValues.t = currentTime;
        data.push_back(currentValues);

        I1 = dI1;
        Id = dId;
        Ii = dIi;
        Il = dIl;
        Ip = dIp;
        It = dIt;
        Xt = dXt;
        Yt = dYt;
        gp = dgp;
        gt = dgt;
        Ipo = dIpo;
        fgut = dfgut;
        fliq = dfliq;
        fsol = dfsol;

        h1 = hopt1;
        h2 = hopt2;
        h3 = hopt3;
        h4 = hopt4;
        h5 = hopt5;
        h6 = hopt6;
        h7 = hopt7;
        h8 = hopt8;
        h9 = hopt9;
        h10 = hopt10;
        h11 = hopt11;
        h12 = hopt12;
        h13 = hopt13;
        h14 = hopt14;

        currentTime += step;
        t = currentTime; // увеличиваем время

    }

    QVector<double> tick;
    QVector<double> CGB;
    QVector<double> Ipg;

    for (DataZ dt : data)       // чтение и выгрузка нуэных результов из структуры
    {
        cout << dt.t << "\t" <<  dt.gp/1.8 << "\t" <<  dt.Ip*20 << "\t" << endl;
        tick.append(dt.t);
        CGB.append(dt.gp/Kvect.at(23));
        Ipg.append(dt.Ip*20);
    }

    cout << "Готово! Кол-во шагов: " << stepCount << endl;
    cout << "Размер: " << sizeof(DataZ)*data.size() << " bytes" << endl;

    /* настройка графика главного меню*/

    // to do:                                               выгрузить с матлаб графики и сравнить их

    ui->customPlot->addGraph();
    QPen pen0(Qt::black);
    pen0.setWidth(3);
    ui->customPlot->xAxis->setAutoTickStep(true);
    ui->customPlot->graph(0)->setPen(pen0);
    ui->customPlot->graph(0)->setLineStyle(QCPGraph::lsLine);               // График в виде чего-то там стиль линий

    ui->customPlot->setInteraction(QCP::iRangeDrag, true);                 // Отключаем взаимодействие перетаскивания графика
    ui->customPlot->setInteraction(QCP::iRangeZoom,true);
    ui->customPlot->addGraph();
    QPen pen3(Qt::blue);
    pen3.setWidth(5);                                                       // установить нужную толщину
    ui->customPlot->graph(1)->setPen(pen3);                                 // Устанавливаем цвет графика
    ui->customPlot->graph(1)->setLineStyle(QCPGraph::LineStyle::lsLine);    // График в виде чего-то там стиль

    ui->customPlot->addGraph();
    QPen pen4(Qt::red);
    pen3.setWidth(5);                                                       // установить нужную толщину
    ui->customPlot->graph(2)->setPen(pen4);                                 // Устанавливаем цвет графика
    ui->customPlot->graph(2)->setLineStyle(QCPGraph::LineStyle::lsLine);    // График в виде чего-то там стиль

    ui->customPlot->graph(0)->setData(tick ,CGB);
    ui->customPlot->graph(1)->setData(tick ,Ipg);
    ui->customPlot->graph(2)->setData(tick ,grafVbas);
    ui->customPlot->rescaleAxes();
    ui->customPlot->replot();
}


MainWindow::~MainWindow()
{
    delete ui;
}

