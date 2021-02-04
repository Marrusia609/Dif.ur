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
double Vmxx =0.087;     // mg/kg/min per pmol/l   insulin sens +
double kp1 = 2.76;      // mg/kg/min +



double m2 = 90;         // + М ?
double tms = 0;     // +

/* Разные названия переменных с начальным массивом Z ? */
double g0 = 90;           // gp?
double gt0 = 169.72;      // + можно ли их все сразу заменить на те
double Ip0 = 5.2040;      // 0.15??  +

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
    double  h = 2.5;

    vector<DataZ> data;
    DataZ start;
    start.t = a;

    int stepCount = (b / step) + 1;
    double currentTime = a + step;

    double K[7];
    double dz;
    double err;
    double s;
    double hmin = 0.01;
    double hmax = 0.25;
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


    /************************************************                        Дано                                           *************************************************/

    /* bgDynam */               // присвоение стартовых значений в структуру для расчёта
    double Ginit=122.00;        // начальный уровень гликемии
    g0 = Ginit*1.8;                 // +
    start.gp = g0;                  // +
    start.Il = 2.61;                // +
    start.Ip = Ip0;                 //55.12/20 +
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


    // Initial

    double Gpb;

    double Ipb;
    double Ib;
    double Gb;
    // Initial
    Gpb=g0; // mg/kg
    Ipb=Ip0; // pmol/kg
    Ib=Ipb*20; // pmol/l
    Gb=Gpb/1.8; // mg/dl



    // Insulin

    double m1=0.190; //min^-1 +
    double m6=0.6471; //dimensionless +
    double m3t;
    // Insulin
    m3t=(m6*m1)/(1- m6); // +

    Vg=1.8; // plasma per BW, dl/kg +

    // Glucose compartments
    double k1gg=0.065; // +
    double k2gg=0.079; // +

    // Fluctuations
    double vbas;
    //double Abol;
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
    /*************************************************                      Конец дано)                              ********************************************************/

    for (int j = 1; j <= stepCount; j++)
    {

        /* */
        // Fluctuations
        vbas=Vbas*100; // pmol/min
        //точки из матлаб убрал в последующих двух уравнениях после второй скобки множителя
        bol1=1/Ti1*Dbol1*(1./(1+exp(-3*(t+10-ti1))))*(1./(1+exp(-3*(-10+ti1-t+Ti1))));
        bol2=1/Ti2*Dbol2*(1./(1+exp(-3*(t+10-ti2))))*(1./(1+exp(-3*(-10+ti2-t+Ti2))));
        bol3=1/Ti3*Dbol3*(1./(1+exp(-3*(t+10-ti3))))*(1./(1+exp(-3*(-10+ti3-t+Ti3))));
        vbol=6000*(bol1+bol2+bol3);


        vm=Dig/Tm*(1./(1+exp(-3*(t-tm1))))*(1./(1+exp(-3*(-(t-tm1-Tm)))));



       // Gtb=(Fcns-EGPb+k1gg*Gpb)/k2gg;

        /* ДУ */

        // Решаем dgp через Рунге-Кутта
        //dgp=EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*fgut-Kvect.at(26)-Kvect.at(31)*(gp-Kvect.at(32))*Heavi1-Kvect.at(24)*gp+Kvect.at(25)*gt; // plasma glucose mg/dl
        // dgp = f(fgut,gt,Id,Ipo,fgut,gp)


//        double EGP;
//        EGP=(Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)*Heavi3;


        double Heavi1;
        if( (gp-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( (gp-Kvect.at(32)) < 0){
             Heavi1 = 0;}

        double Heavi3;
        if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo) < 0){
             Heavi3 = 0;}

        K[0] = (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*fgut-Kvect.at(26)-Kvect.at(31)*(gp-Kvect.at(32))*Heavi1-Kvect.at(24)*gp+Kvect.at(25)*gt;

        if( ((gp+ h*K[0]/2.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ h*K[0]/2.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[0]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[0]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0)) < 0){
             Heavi3 = 0;}
        K[1] = (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[0]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h/2.0)-Kvect.at(26)-Kvect.at(31)*((gp+ h*K[0]/2.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h*K[0]/2.0)+Kvect.at(25)*(gt+ h/2.0);

        if( ((gp+ h*K[1]/2.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ h*K[1]/2.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[1]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[1]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0)) < 0){
             Heavi3 = 0;}
        K[2] = (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[1]/2.0)-Kvect.at(13)*(Id+ h/2.0)-Kvect.at(14)*(Ipo+ h/2.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h/2.0)-Kvect.at(26)-Kvect.at(31)*((gp+ h*K[1]/2.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h*K[1]/2.0)+Kvect.at(25)*(gt+ h/2.0);

        if( ((gp+ h*K[2])-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ h*K[2])-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[2])-Kvect.at(13)*(Id+ h)-Kvect.at(14)*(Ipo+ h))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[2])-Kvect.at(13)*(Id+ h)-Kvect.at(14)*(Ipo+ h)) < 0){
             Heavi3 = 0;}
        K[3] = (Kvect.at(11)-Kvect.at(12)*(gp+ h*K[2])-Kvect.at(13)*(Id+ h)-Kvect.at(14)*(Ipo+ h))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h)-Kvect.at(26)-Kvect.at(31)*((gp+ h*K[2])-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h*K[2])+Kvect.at(25)*(gt+ h);

        dgp   = gp + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dgp через Дормана-Принса
/*
        double Heavi1;
        if( (gp-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( (gp-Kvect.at(32)) < 0){
             Heavi1 = 0;}

        double Heavi3;
        if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo) < 0){
             Heavi3 = 0;}

        K[0] = (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*fgut-Kvect.at(26)-Kvect.at(31)*(gp-Kvect.at(32))*Heavi1-Kvect.at(24)*gp+Kvect.at(25)*gt;

        if( ((gp+ h1*K[0]/5.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((gp+ h1*K[0]/5.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h1*K[0]/5.0)-Kvect.at(13)*(Id+ h1/5.0)-Kvect.at(14)*(Ipo+ h1/5.0))  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*(gp+ h1*K[0]/5.0)-Kvect.at(13)*(Id+ h1/5.0)-Kvect.at(14)*(Ipo+ h1/5.0)) < 0){
             Heavi3 = 0;}
        K[1] = (Kvect.at(11)-Kvect.at(12)*(gp+ h1*K[0]/5.0)-Kvect.at(13)*(Id+ h1/5.0)-Kvect.at(14)*(Ipo+ h1/5.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1/5.0)-Kvect.at(26)-Kvect.at(31)*((gp+ h1*K[0]/5.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h1*K[0]/5.0)+Kvect.at(25)*(gt+ h1/5.0);
*/

        // Решаем dgt через Рунге-Кутта

        // dgt=-((Kvect.at(27)+Kvect.at(28)*Xt)*gt)/(Kvect.at(29)+gt)+Kvect.at(24)*gp-Kvect.at(25)*gt;
        // dgt = f(Xt,gp,gt)

        K[0] = -((Kvect.at(27)+Kvect.at(28)*Xt)*gt)/(Kvect.at(29)+gt)+Kvect.at(24)*gp-Kvect.at(25)*gt;
        K[1] = -((Kvect.at(27)+Kvect.at(28)*(Xt + h/2.0))*(gt + h*K[0]/2.0))/(Kvect.at(29)+(gt + h*K[0]/2.0))+Kvect.at(24)*(gp + h/2.0)-Kvect.at(25)*(gt + h*K[0]/2.0);
        K[2] = -((Kvect.at(27)+Kvect.at(28)*(Xt + h/2.0))*(gt + h*K[1]/2.0))/(Kvect.at(29)+(gt + h*K[1]/2.0))+Kvect.at(24)*(gp + h/2.0)-Kvect.at(25)*(gt + h*K[1]/2.0);
        K[3] = -((Kvect.at(27)+Kvect.at(28)*(Xt + h))*(gt + h*K[2]))/(Kvect.at(29)+(gt + h*K[2]))+Kvect.at(24)*(gp + h)-Kvect.at(25)*(gt + h*K[2]);
        dgt   = gt + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dgt через Дормана-Принса
/*

*/

        // Решаем dI1 через Рунге-Кутта

        // dI1=-Kvect.at(15)*(I1-Ip/Kvect.at(4));
        // dI1 = f(Ip,I1)

        K[0] = -Kvect.at(15)*(I1-Ip/Kvect.at(4));
        K[1] = -Kvect.at(15)*((I1 + h*K[0]/2.0)-(Ip+ h/2.0)/Kvect.at(4));
        K[2] = -Kvect.at(15)*((I1 + h*K[1]/2.0)-(Ip+ h/2.0)/Kvect.at(4));
        K[3] = -Kvect.at(15)*((I1 + h*K[2])-(Ip+ h)/Kvect.at(4));
        dI1   = I1 + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dI1 через Дормана-Принса
/*

*/

        // Решаем dId через Рунге-Кутта
        // dId=-Kvect.at(15)*(Id-I1);
        // dId = f(I1,Id)?
        K[0] = -Kvect.at(15)*(Id-I1);
        K[1] = -Kvect.at(15)*((Id + h*K[0]/2.0)-(I1+h/2.0));
        K[2] = -Kvect.at(15)*((Id + h*K[1]/2.0)-(I1+h/2.0));
        K[3] = -Kvect.at(15)*((Id + h*K[2])-(I1+h));
        dId   = Id + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dId через Дормана-Принса
/*
        K[0] = -Kvect.at(15)*(Id-I1);
        K[1] = -Kvect.at(15)*((Id + h4*K[0]/5.0)-(I1+h4/5.0));
        K[2] = -Kvect.at(15)*((Id + 3*h4*K[0]/40.0 + 9*h4*K[1]/40.0)-(I1+3*h4/10.0));
        K[3] = -Kvect.at(15)*((Id + 44*h4*K[0]/45.0 + (-56*h4*K[1]/15.0) + 32*h4*K[2]/9.0)-(I1+4*h4/5.0));
        K[4] = -Kvect.at(15)*((Id + 19372*h4*K[0]/6561.0 + (-25360*h4*K[1]/2187.0) + 64448*h4*K[2]/6561.0 + (-212*h4*K[3]/729.0))-(I1+8*h4/9.0));
        K[5] = -Kvect.at(15)*((Id + 9017*h4*K[0]/3168.0 + (-355*h4*K[1]/33.0) + 46732*h4*K[2]/5247.0 + (49*h4*K[3]/176.0) + (-5103*h4*K[4]/18656.0))-(I1+h4));
        K[6] = -Kvect.at(15)*((Id + 35*h4*K[0]/384.0 + 500*h4*K[2]/1113.0 + (125*h4*K[3]/192.0) + (-2187*h4*K[4]/6784.0) + 11*h4*K[5]/84.0)-(I1+h4));
        dId   = Id + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h4;                                                // во многих источниках домножается на h, но в одном нету такого.
        dz = Id + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h4;
        err = abs(dz-dId);
        s = pow(eps*h4/(2*err),1/5);
        hopt4 = s*h4;
        if( hopt4 < hmin) hopt4 = hmin;
        else if(hopt4 > hmax) hopt4 = hmax;
*/

        // Решаем dXt через Рунге-Кутта
        // dXt=-Kvect.at(30)*Xt+Kvect.at(30)*((Ip/Kvect.at(4))-Kvect.at(1))*Heavi2;
        // dXt = f(Ip,Xt)?
        double Heavi2;
        if( ((Ip/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( ((Ip/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[0] = -Kvect.at(30)*Xt+Kvect.at(30)*((Ip/Kvect.at(4))-Kvect.at(1))*Heavi2;

        if( (((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[1] = -Kvect.at(30)*(Xt+ h*K[0]/2.0)+Kvect.at(30)*(((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1))*Heavi2;

        if( (((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[2] = -Kvect.at(30)*(Xt + h*K[1]/2.0)+Kvect.at(30)*(((Ip+ h/2.0)/Kvect.at(4))-Kvect.at(1))*Heavi2;

        if( (((Ip+ h)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((Ip+ h)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[3] = -Kvect.at(30)*(Xt+ h*K[2])+Kvect.at(30)*(((Ip+ h)/Kvect.at(4))-Kvect.at(1))*Heavi2;

        dXt   = Xt + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dXt через Дормана-Принса
/*

*/

        // Решаем dIl через Рунге-Кутта
        // dIl=-(Kvect.at(5)+Kvect.at(7))*Il+Kvect.at(6)*Ip;
        // dIl = f(Ip,Il)?

        K[0] = -(Kvect.at(5)+Kvect.at(7))*Il+Kvect.at(6)*Ip;
        K[1] = -(Kvect.at(5)+Kvect.at(7))*(Il+ h*K[0]/2.0)+Kvect.at(6)*(Ip+ h/2.0);
        K[2] = -(Kvect.at(5)+Kvect.at(7))*(Il+ h*K[1]/2.0)+Kvect.at(6)*(Ip+ h/2.0);
        K[3] = -(Kvect.at(5)+Kvect.at(7))*(Il+ h*K[2])+Kvect.at(6)*(Ip+ h);

        dIl   = Il + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIl через Дормана-Принса
/*

*/

        // Решаем dIp через Рунге-Кутта
        // dIp=-Kvect.at(6)*Ip+Kvect.at(5)*Il+Kvect.at(10)/Kvect.at(0)*It-Kvect.at(9)*Ip;
        // dIp = f(Il,It,Ip)?
        K[0] = -Kvect.at(6)*Ip+Kvect.at(5)*Il+Kvect.at(10)/Kvect.at(0)*It-Kvect.at(9)*Ip;
        K[1] = -Kvect.at(6)*(Ip+ h*K[0]/2.0)+Kvect.at(5)*(Il+ h/2.0)+Kvect.at(10)/Kvect.at(0)*(It+ h/2.0)-Kvect.at(9)*(Ip+ h*K[0]/2.0);
        K[2] = -Kvect.at(6)*(Ip+ h*K[1]/2.0)+Kvect.at(5)*(Il+ h/2.0)+Kvect.at(10)/Kvect.at(0)*(It+ h/2.0)-Kvect.at(9)*(Ip+ h*K[1]/2.0);
        K[3] = -Kvect.at(6)*(Ip+ h*K[2])+Kvect.at(5)*(Il+ h)+Kvect.at(10)/Kvect.at(0)*(It+ h)-Kvect.at(9)*(Ip+ h*K[2]);

        dIp   = Ip + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIp через Дормана-Принса
/*

*/

        // Решаем dfgut через Рунге-Кутта
        // dfgut=-Kvect.at(17)*fgut+kgut*fliq;
        // dfgut = f(fliq,fsol,fgut)?

//        double kgut;
//        kgut=Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2);

        K[0] = -Kvect.at(17)*fgut+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2)*fliq;
        K[1] = -Kvect.at(17)*(fgut+ h*K[0]/2.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h/2.0)+(fliq+ h/2.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h/2.0)+(fliq+ h/2.0)-Kvect.at(22)*Dig))+2)*(fliq+ h/2.0);
        K[2] = -Kvect.at(17)*(fgut+ h*K[1]/2.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h/2.0)+(fliq+ h/2.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h/2.0)+(fliq+ h/2.0)-Kvect.at(22)*Dig))+2)*(fliq+ h/2.0);
        K[3] = -Kvect.at(17)*(fgut+ h*K[2])+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h)+(fliq+ h)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h)+(fliq+ h)-Kvect.at(22)*Dig))+2)*(fliq+ h);

        dfgut   = fgut + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfgut через Дормана-Принса
/*

*/

        // Решаем dfliq через Рунге-Кутта
        // dfliq=-kgut*fliq+Kvect.at(18)*fsol;
        // dfliq = f(fliq,fsol)?

        K[0] = -(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2))*fliq+Kvect.at(18)*fsol;
        K[1] = -(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h/2.0)+(fliq+ h*K[0]/2.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h/2.0)+(fliq+ h*K[0]/2.0)-Kvect.at(22)*Dig))+2))*(fliq+ h*K[0]/2.0)+Kvect.at(18)*(fsol+ h/2.0);
        K[2] = -(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h/2.0)+(fliq+ h*K[1]/2.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h/2.0)+(fliq+ h*K[1]/2.0)-Kvect.at(22)*Dig))+2))*(fliq+ h*K[1]/2.0)+Kvect.at(18)*(fsol+ h/2.0);
        K[2] = -(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h)+(fliq+ h*K[2])-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h)+(fliq+ h*K[2])-Kvect.at(22)*Dig))+2))*(fliq+ h*K[2])+Kvect.at(18)*(fsol+ h);

        dfliq   = fliq + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfliq через Дормана-Принса
/*

*/
        // Решаем dfsol через Рунге-Кутта
        // dfsol=-Kvect.at(18)*fsol+vm;
        // dfsol = f(fsol,vm)?

        K[0] = -Kvect.at(18)*fsol+vm;
        K[1] = -Kvect.at(18)*(fsol+ h*K[0]/2.0)+(vm+ h/2.0);
        K[2] = -Kvect.at(18)*(fsol+ h*K[1]/2.0)+(vm+ h/2.0);
        K[3] = -Kvect.at(18)*(fsol+ h*K[2])+(vm+ h);

        dfsol   = fsol + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfsol через Дормана-Принса
/*

*/
        // Решаем dIpo через Рунге-Кутта
        // dIpo=-Kvect.at(33)*Ipo+(Yt+Kvect.at(3))*Heavi5+(Yt+Kvect.at(3))*(Heavi6);
        // dIpo = f(Yt,dgp,Ipo)?
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
        K[0] = -Kvect.at(33)*Ipo+(Yt+Kvect.at(3))*Heavi5+(Yt+Kvect.at(3))*(Heavi6);

        if( ((dgp+ h/2.0)) >= 0){
             Heavi5 = 1;}
        if( ((dgp+ h/2.0)) < 0){
             Heavi5 = 0;}

        if( (-(dgp+ h/2.0)) >= 0){
             Heavi6 = 1;}
        if( (-(dgp+ h/2.0)) < 0){
             Heavi6 = 0;}
        K[1] = -Kvect.at(33)*(Ipo+ h*K[0]/2.0)+((Yt+ h/2.0)+Kvect.at(3))*Heavi5+((Yt+ h/2.0)+Kvect.at(3))*(Heavi6);

        if( ((dgp+ h/2.0)) >= 0){
             Heavi5 = 1;}
        if( ((dgp+ h/2.0)) < 0){
             Heavi5 = 0;}

        if( (-(dgp+ h/2.0)) >= 0){
             Heavi6 = 1;}
        if( (-(dgp+ h/2.0)) < 0){
             Heavi6 = 0;}
        K[2] = -Kvect.at(33)*(Ipo+ h*K[1]/2.0)+((Yt+ h/2.0)+Kvect.at(3))*Heavi5+((Yt+ h/2.0)+Kvect.at(3))*(Heavi6);

        if( ((dgp+ h)) >= 0){
             Heavi5 = 1;}
        if( ((dgp+ h)) < 0){
             Heavi5 = 0;}

        if( (-(dgp+ h)) >= 0){
             Heavi6 = 1;}
        if( (-(dgp+ h)) < 0){
             Heavi6 = 0;}
        K[3] = -Kvect.at(33)*(Ipo+ h*K[2])+((Yt+ h)+Kvect.at(3))*Heavi5+((Yt+ h)+Kvect.at(3))*(Heavi6);

        dIpo   = Ipo + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfsol через Дормана-Принса
/*

*/

        // Решаем dYt через Рунге-Кутта
        // dYt=-Kvect.at(34)*(Yt-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*Yt-Kvect.at(34)*Kvect.at(3))*(Heavi8);
        // dYt = f(gp,Yt)?
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
        K[0] = -Kvect.at(34)*(Yt-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*Yt-Kvect.at(34)*Kvect.at(3))*(Heavi8);

        if( (Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[1] = -Kvect.at(34)*((Yt+ h*K[0]/2.0)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ h*K[0]/2.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8);

        if( (Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[2] = -Kvect.at(34)*((Yt+ h*K[1]/2.0)-Kvect.at(35)*((gp+ h/2.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ h*K[1]/2.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8);

        if( (Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[3] = -Kvect.at(34)*((Yt+ h*K[2])-Kvect.at(35)*((gp+ h)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ h*K[2])-Kvect.at(34)*Kvect.at(3))*(Heavi8);

        dYt   = Yt + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfsol через Дормана-Принса
/*

*/

        // Решаем dIt через Рунге-Кутта
        // dIt=Kvect.at(8)*Ii-Kvect.at(10)*It;
        // dIt = f(Ii,It)?
        K[0] = Kvect.at(8)*Ii-Kvect.at(10)*It;
        K[1] = Kvect.at(8)*(Ii+ h/2.0)-Kvect.at(10)*(It+ h*K[0]/2.0);
        K[2] = Kvect.at(8)*(Ii+ h/2.0)-Kvect.at(10)*(It+ h*K[1]/2.0);
        K[3] = Kvect.at(8)*(Ii+ h)-Kvect.at(10)*(It+ h*K[2]);

        dIt   = It + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfsol через Дормана-Принса
/*

*/

        // Решаем dIi через Рунге-Кутта
        // dIi=-Kvect.at(8)*Ii+vbas+vbol;
        // dIi = f(vbol,Ii)?
        K[0] = -Kvect.at(8)*Ii+vbas+vbol;
        K[1] = -Kvect.at(8)*(Ii+ h*K[0]/2.0)+vbas+(vbol+ h/2.0);
        K[2] = -Kvect.at(8)*(Ii+ h*K[1]/2.0)+vbas+(vbol+ h/2.0);
        K[3] = -Kvect.at(8)*(Ii+ h*K[2])+vbas+(vbol+ h);

        dIi   = Ii + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfsol через Дормана-Принса
/*

*/

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
    ui->customPlot->addGraph();
    QPen pen0(Qt::black);
    pen0.setWidth(3);
    ui->customPlot->xAxis->setAutoTickStep(true);
    ui->customPlot->xAxis->setTickLabelType(QCPAxis::ltDateTime);           // Подпись координат по Оси X в качестве Даты и Времени
    ui->customPlot->xAxis->setDateTimeFormat("hh:mm:ss");
    ui->customPlot->graph(0)->setPen(pen0);
    ui->customPlot->graph(0)->setLineStyle(QCPGraph::lsLine);               // График в виде чего-то там стиль линий

    ui->customPlot->setInteraction(QCP::iRangeDrag, false);                 // Отключаем взаимодействие перетаскивания графика
    ui->customPlot->addGraph();
    QPen pen3(Qt::blue);
    pen3.setWidth(5);                                                       // установить нужную толщину
    ui->customPlot->graph(1)->setPen(pen3);                                 // Устанавливаем цвет графика
    ui->customPlot->graph(1)->setLineStyle(QCPGraph::LineStyle::lsLine);    // График в виде чего-то там стиль

    ui->customPlot->graph(0)->setData(tick ,CGB);
    ui->customPlot->graph(1)->setData(tick ,Ipg);
    ui->customPlot->rescaleAxes();
    ui->customPlot->replot();
}


MainWindow::~MainWindow()
{
    delete ui;
}

