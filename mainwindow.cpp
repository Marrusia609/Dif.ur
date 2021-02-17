#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <vector>
#include <QVector>
#include <math.h>
#include "qcustomplot.h"

using namespace std;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    MatDate();
}
MainWindow::~MainWindow()
{
    delete ui;
}
void MainWindow::MatDate()
{
    /* MATLAB данные ode45 выгрузка */
    QTextStream in;
    int ind = 0;                            // индекс для записи в массивы

    QFile fileMat("D:/Qt/pr/untitled/Dif.ur/matlab/progMatExplDate1.csv");
    if (!fileMat.open(QFile::ReadOnly)){
        qDebug() << "File not exists";
    }
    else{

        in.setDevice(&fileMat);
        double tr=0;

        while( !in.atEnd()){
                QString read = in.readLine();

                ind=read.indexOf("\n");
                tr=read.left(ind).toDouble();
                cgbMatlab.append(tr);
                read.remove(0,ind+1);
        }
    fileMat.close();
    }
    QFile fileMat2("D:/Qt/pr/untitled/Dif.ur/matlab/progMatExplDate2.csv");
    if (!fileMat2.open(QFile::ReadOnly)){
        qDebug() << "File not exists";
    }
    else{

        in.setDevice(&fileMat2);
        double tr=0;

        while( !in.atEnd()){
                QString read = in.readLine();

                ind=read.indexOf("\n");
                tr=read.left(ind).toDouble();
                insMatlab.append(tr);
                read.remove(0,ind+1);
        }
    fileMat2.close();
    }

    /* Обработка K[37] */
    QString read;
    QString constPat ("D:/Qt/pr/untitled/Dif.ur/matlab/PatConst.csv");
    QFile file(constPat);
    if (!file.open(QFile::ReadOnly)){
        qDebug() << "File not exists";
    }
    else{

        in.setDevice(&file);
        in >> read;
        file.close();

        double tr=0;

        for ( int i=0; i < 37; i++ ){

                ind=read.indexOf(",");
                tr=read.left(ind).toDouble();
                Kvect.append(tr);
                read.remove(0,ind+1);

        }
        GlobalCircle();
    }
}
void MainWindow::GlobalCircle()
{

/************************************************              Дано             *************************************************/

    // векторы для выгрузки на графики
    QVector<double> cgbMatlab;
    QVector<double> insMatlab;

    double t0global = 0;    // начальное время
    double t1global = 720 -5;  // конечное время

    Vbas = 1.2238;   // доза базального инсулина  +
    VbasDP = 1.2238;   // доза базального инсулина  +
    bas[0] = Vbas;
    bas[1] = 140.0;
    bas[2] = 110.0;
    bas[3] = 70.0;
    bas[4] = 1.5;
    bas[6] = 180.0;
    Vb = bas[0];

    basDP[0] = VbasDP;
    basDP[1] = 140.0;
    basDP[2] = 110.0;
    basDP[3] = 70.0;
    basDP[4] = 1.5;
    basDP[6] = 180.0;
    VbDP = basDP[0];

    grafVbas.append(Vbas);
    grafVbasDP.append(VbasDP);
    tickBas.append(t0global);
    tickBasDP.append(t0global);

    double Ginit=122.00;            // начальный уровень гликемии
    dZ[0] = Ginit*1.8;      // gp
    dZ[1] = 104.08;         // Id
    dZ[2] = 3.08;           // Ipo
    dZ[3] = 0.0;            // fgut
    dZ[4] = 0.0;            // fsol
    dZ[5] = 0.0;            // fliq
    dZ[6] = 169.72;         // gt
    dZ[7] = 0.0;            // Xt
    dZ[8] = 104.08;         // I1
    dZ[9] = 5.2040;         // Ip
    dZ[10] = 10830;         // It
    dZ[11] = 2.61;          // Il
    dZ[12] = 0.0;           // Yt
    dZ[13] = 4120;          // Ii

    DdZ[0] = Ginit*1.8;      // gp
    DdZ[1] = 104.08;         // Id
    DdZ[2] = 3.08;           // Ipo
    DdZ[3] = 0.0;            // fgut
    DdZ[4] = 0.0;            // fsol
    DdZ[5] = 0.0;            // fliq
    DdZ[6] = 169.72;         // gt
    DdZ[7] = 0.0;            // Xt
    DdZ[8] = 104.08;         // I1
    DdZ[9] = 5.2040;         // Ip
    DdZ[10] = 10830;         // It
    DdZ[11] = 2.61;          // Il
    DdZ[12] = 0.0;           // Yt
    DdZ[13] = 4120;          // Ii

        /* Глобальный цикл */
    for (int i = 1; t0global <= t1global; i++){

        t = t0global;
        tDP = t0global;
        int *c = &i;
        RungeKutta(c);
        //DormanPrince(c);

//        for (double j = a; j <= b;)

//        {
//                /* bolus */
//            double Ti1 = 10;
//            double Ti2 = 10;
//            double ti2 = 60;
//            double ti1 = 30;
//            double Ti3 = 10;
//            double ti3 = 10;
//            double tm1 = 60; // tm  // время приёма пищи
//            double Tm = 20;         // длительность приёма пищи

//            // Fluctuations
//            vbas=Vbas*100;      // pmol/min
//            bol1=1/Ti1*Dbol1*(1./(1+exp(-3*(t+10-ti1))))*(1./(1+exp(-3*(-10+ti1-t+Ti1))));
//            bol2=1/Ti2*Dbol2*(1./(1+exp(-3*(t+10-ti2))))*(1./(1+exp(-3*(-10+ti2-t+Ti2))));
//            bol3=1/Ti3*Dbol3*(1./(1+exp(-3*(t+10-ti3))))*(1./(1+exp(-3*(-10+ti3-t+Ti3))));

//            vbol=6000*(bol1+bol2+bol3);

//            vm=Dig/Tm*(1./(1+exp(-3*(t-tm1))))*(1./(1+exp(-3*(-(t-tm1-Tm)))));

//            /* ДУ */

//            // Решение dgp через Дормана-Принса

//            double Heavi1;
//            double Heavi3;
//            if( (gp-Kvect.at(32))>= 0){
//                 Heavi1 = 1;}
//            if( (gp-Kvect.at(32)) < 0){
//                 Heavi1 = 0;}

//            if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)  >= 0){
//                 Heavi3 = 1;}
//            if( (Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo) < 0){
//                 Heavi3 = 0;}
//            double EGP;
//            EGP=(Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)*Heavi3;

//            //K[0] = h1*((Kvect.at(11)-Kvect.at(12)*gp-Kvect.at(13)*Id-Kvect.at(14)*Ipo)*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*fgut-Kvect.at(26)-Kvect.at(31)*(gp-Kvect.at(32))*Heavi1-Kvect.at(24)*gp+Kvect.at(25)*gt);
//            K[0] = h1*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*fgut-Kvect.at(26)-Kvect.at(31)*(gp-Kvect.at(32))*Heavi1-Kvect.at(24)*gp+Kvect.at(25)*gt);

//            if( ((gp+ h1*K[0]/5.0)-Kvect.at(32))>= 0){
//                 Heavi1 = 1;}
//            if( ((gp+ h1*K[0]/5.0)-Kvect.at(32)) < 0){
//                 Heavi1 = 0;}
//            /*if( (Kvect.at(11)-Kvect.at(12)*(gp+ h1*K[0]/5.0)-Kvect.at(13)*(Id+ h1/5.0)-Kvect.at(14)*(Ipo+ h1/5.0))  >= 0){
//                 Heavi3 = 1;}
//            if( (Kvect.at(11)-Kvect.at(12)*(gp+ h1*K[0]/5.0)-Kvect.at(13)*(Id+ h1/5.0)-Kvect.at(14)*(Ipo+ h1/5.0)) < 0){
//                 Heavi3 = 0;}*/
//            //K[1] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ h1*K[0]/5.0)-Kvect.at(13)*(Id+ h1/5.0)-Kvect.at(14)*(Ipo+ h1/5.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1/5.0)-Kvect.at(26)-Kvect.at(31)*((gp+ h1*K[0]/5.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h1*K[0]/5.0)+Kvect.at(25)*(gt+ h1/5.0));
//            K[1] = h1*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1/5.0)-Kvect.at(26)-Kvect.at(31)*((gp+ h1*K[0]/5.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ h1*K[0]/5.0)+Kvect.at(25)*(gt+ h1/5.0));

//            if( ((gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(32))>= 0){
//                 Heavi1 = 1;}
//            if( ((gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(32)) < 0){
//                 Heavi1 = 0;}
//            /*if( (Kvect.at(11)-Kvect.at(12)*(gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(13)*(Id+ 3*h1/10.0)-Kvect.at(14)*(Ipo+ 3*h1/10.0))  >= 0){
//                 Heavi3 = 1;}
//            if( (Kvect.at(11)-Kvect.at(12)*(gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(13)*(Id+ 3*h1/10.0)-Kvect.at(14)*(Ipo+ 3*h1/10.0)) < 0){
//                 Heavi3 = 0;}*/
//            //K[2] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(13)*(Id+ 3*h1/10.0)-Kvect.at(14)*(Ipo+ 3*h1/10.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ 3*h1/10.0)-Kvect.at(26)-Kvect.at(31)*((gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)+Kvect.at(25)*(gt+ 3*h1/10.0));
//            K[2] = h1*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ 3*h1/10.0)-Kvect.at(26)-Kvect.at(31)*((gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 3*h1*K[0]/40.0 + 9*h1*K[1]/40.0)+Kvect.at(25)*(gt+ 3*h1/10.0));

//            if( ((gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(32))>= 0){
//                 Heavi1 = 1;}
//            if( ((gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(32)) < 0){
//                 Heavi1 = 0;}
//            /*if( (Kvect.at(11)-Kvect.at(12)*(gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(13)*(Id+ 4*h1/5.0)-Kvect.at(14)*(Ipo+ 4*h1/5.0))  >= 0){
//                 Heavi3 = 1;}
//            if( (Kvect.at(11)-Kvect.at(12)*(gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(13)*(Id+ 4*h1/5.0)-Kvect.at(14)*(Ipo+ 4*h1/5.0)) < 0){
//                 Heavi3 = 0;}*/
//            //K[3] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(13)*(Id+ 4*h1/5.0)-Kvect.at(14)*(Ipo+ 4*h1/5.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ 4*h1/5.0)-Kvect.at(26)-Kvect.at(31)*((gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)+Kvect.at(25)*(gt+ 4*h1/5.0));
//            K[3] = h1*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ 4*h1/5.0)-Kvect.at(26)-Kvect.at(31)*((gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 44*h1*K[0]/45.0 + (-56*h1*K[1]/15.0) + 32*h1*K[2]/9.0)+Kvect.at(25)*(gt+ 4*h1/5.0));

//            if( ((gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(32))>= 0){
//                 Heavi1 = 1;}
//            if( ((gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(32)) < 0){
//                 Heavi1 = 0;}
//            /*if( (Kvect.at(11)-Kvect.at(12)*(gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(13)*(Id+ 8*h1/9.0)-Kvect.at(14)*(Ipo+ 8*h1/9.0))  >= 0){
//                 Heavi3 = 1;}
//            if( (Kvect.at(11)-Kvect.at(12)*(gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(13)*(Id+ 8*h1/9.0)-Kvect.at(14)*(Ipo+ 8*h1/9.0)) < 0){
//                 Heavi3 = 0;}*/
//            //K[4] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(13)*(Id+ 8*h1/9.0)-Kvect.at(14)*(Ipo+ 8*h1/9.0))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ 8*h1/9.0)-Kvect.at(26)-Kvect.at(31)*((gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )+Kvect.at(25)*(gt+ 8*h1/9.0));
//            K[4] = h1*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ 8*h1/9.0)-Kvect.at(26)-Kvect.at(31)*((gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 19372*h1*K[0]/6561.0 + (-25360*h1*K[1]/2187.0) + 64448*h1*K[2]/6561.0 + (-212*h1*K[3]/729.0) )+Kvect.at(25)*(gt+ 8*h1/9.0));

//            if( ((gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(32))>= 0){
//                 Heavi1 = 1;}
//            if( ((gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(32)) < 0){
//                 Heavi1 = 0;}
//            /*if( (Kvect.at(11)-Kvect.at(12)*(gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1))  >= 0){
//                 Heavi3 = 1;}
//            if( (Kvect.at(11)-Kvect.at(12)*(gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1)) < 0){
//                 Heavi3 = 0;}*/
//            //K[5] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1)-Kvect.at(26)-Kvect.at(31)*((gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))+Kvect.at(25)*(gt+ h1));
//            K[5] = h1*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1)-Kvect.at(26)-Kvect.at(31)*((gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 9017*h1*K[0]/3168.0 + (-355*h1*K[1]/33.0) + 46732*h1*K[2]/5247.0 + (49*h1*K[3]/176.0) + (-5103*h1*K[4]/18656.0))+Kvect.at(25)*(gt+ h1));

//            if( ((gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(32))>= 0){
//                 Heavi1 = 1;}
//            if( ((gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(32)) < 0){
//                 Heavi1 = 0;}
//            /*if( (Kvect.at(11)-Kvect.at(12)*(gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1))  >= 0){
//                 Heavi3 = 1;}
//            if( (Kvect.at(11)-Kvect.at(12)*(gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1)) < 0){
//                 Heavi3 = 0;}*/
//            //K[6] = h1*((Kvect.at(11)-Kvect.at(12)*(gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(13)*(Id+ h1)-Kvect.at(14)*(Ipo+ h1))*Heavi3+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1)-Kvect.at(26)-Kvect.at(31)*((gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)+Kvect.at(25)*(gt+ h1));
//            K[6] = h1*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(fgut+ h1)-Kvect.at(26)-Kvect.at(31)*((gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(gp+ 35*h1*K[0]/384.0 + 500*h1*K[2]/1113.0 + (125*h1*K[3]/192.0) + (-2187*h1*K[4]/6784.0) + 11*h1*K[5]/84.0)+Kvect.at(25)*(gt+ h1));

//            dgp   = gp + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);     // во многих источниках домножается на h, но в одном нету такого.
//            dz = gp + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dgp);
//            s = pow(eps*h1/(2*err),1/5);
//            hopt1 = s*h1;
//            if( hopt1 < hmin)  hopt1 = hmin;
//            else if(hopt1 > hmax) hopt1 = hmax;

//            // Решение dgt через Дормана-Принса

//            K[0] = h2*(-((Kvect.at(27)+Kvect.at(28)*Xt)*gt)/(Kvect.at(29)+gt)+Kvect.at(24)*gp-Kvect.at(25)*gt);
//            K[1] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + h2/5.0))*(gt + h2*K[0]/5.0))/(Kvect.at(29)+(gt + h2*K[0]/5.0))+Kvect.at(24)*(gp + h2/5.0)-Kvect.at(25)*(gt + h2*K[0]/5.0));
//            K[2] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + 3*h2/10.0))*(gt + 3*h2*K[0]/40.0 + 9*h2*K[1]/40.0))/(Kvect.at(29)+(gt + 3*h2*K[0]/40.0 + 9*h2*K[1]/40.0))+Kvect.at(24)*(gp + 3*h2/10.0)-Kvect.at(25)*(gt + 3*h2*K[0]/40.0 + 9*h2*K[1]/40.0));
//            K[3] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + 4*h2/5.0))*(gt + 44*h2*K[0]/45.0 + (-56*h2*K[1]/15.0) + 32*h2*K[2]/9.0))/(Kvect.at(29)+(gt + 44*h2*K[0]/45.0 + (-56*h2*K[1]/15.0) + 32*h2*K[2]/9.0))+Kvect.at(24)*(gp + 4*h2/5.0)-Kvect.at(25)*(gt + 44*h2*K[0]/45.0 + (-56*h2*K[1]/15.0) + 32*h2*K[2]/9.0));
//            K[4] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + 8*h2/9.0))*(gt + 19372*h2*K[0]/6561.0 + (-25360*h2*K[1]/2187.0) + 64448*h2*K[2]/6561.0 + (-212*h2*K[3]/729.0) ))/(Kvect.at(29)+(gt + 19372*h2*K[0]/6561.0 + (-25360*h2*K[1]/2187.0) + 64448*h2*K[2]/6561.0 + (-212*h2*K[3]/729.0) ))+Kvect.at(24)*(gp + 8*h2/9.0)-Kvect.at(25)*(gt + 19372*h2*K[0]/6561.0 + (-25360*h2*K[1]/2187.0) + 64448*h2*K[2]/6561.0 + (-212*h2*K[3]/729.0) ));
//            K[5] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + h2))*(gt + 9017*h2*K[0]/3168.0 + (-355*h2*K[1]/33.0) + 46732*h2*K[2]/5247.0 + (49*h2*K[3]/176.0) + (-5103*h2*K[4]/18656.0) ))/(Kvect.at(29)+(gt + 9017*h2*K[0]/3168.0 + (-355*h2*K[1]/33.0) + 46732*h2*K[2]/5247.0 + (49*h2*K[3]/176.0) + (-5103*h2*K[4]/18656.0) ))+Kvect.at(24)*(gp + h2)-Kvect.at(25)*(gt + 9017*h2*K[0]/3168.0 + (-355*h2*K[1]/33.0) + 46732*h2*K[2]/5247.0 + (49*h2*K[3]/176.0) + (-5103*h2*K[4]/18656.0) ));
//            K[6] = h2*(-((Kvect.at(27)+Kvect.at(28)*(Xt + h2))*(gt + 35*h2*K[0]/384.0 + 500*h2*K[2]/1113.0 + (125*h2*K[3]/192.0) + (-2187*h2*K[4]/6784.0) + 11*h2*K[5]/84.0))/(Kvect.at(29)+(gt + 35*h2*K[0]/384.0 + 500*h2*K[2]/1113.0 + (125*h2*K[3]/192.0) + (-2187*h2*K[4]/6784.0) + 11*h2*K[5]/84.0))+Kvect.at(24)*(gp + h2)-Kvect.at(25)*(gt + 35*h2*K[0]/384.0 + 500*h2*K[2]/1113.0 + (125*h2*K[3]/192.0) + (-2187*h2*K[4]/6784.0) + 11*h2*K[5]/84.0));

//            dgt = gt + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);     // во многих источниках домножается на h, но в одном нету такого.
//            dz = gt + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dgt);
//            s = pow(eps*h2/(2*err),1/5);
//            hopt2 = s*h2;
//            if( hopt2 < hmin) hopt2 = hmin;
//            else if(hopt2 > hmax) hopt2 = hmax;

//            // Решение dI1 через Дормана-Принса

//            K[0] = h3*(-Kvect.at(15)*(I1-Ip/Kvect.at(4)));
//            K[1] = h3*(-Kvect.at(15)*((I1 + h3*K[0]/5.0)-(Ip+ h3/5.0)/Kvect.at(4)));
//            K[2] = h3*(-Kvect.at(15)*((I1 + 3*h3*K[0]/40.0 + 9*h3*K[1]/40.0)-(Ip+ 3*h3/10.0)/Kvect.at(4)));
//            K[3] = h3*(-Kvect.at(15)*((I1 + 44*h3*K[0]/45.0 + (-56*h3*K[1]/15.0) + 32*h3*K[2]/9.0)-(Ip+ 4*h3/5.0)/Kvect.at(4)));
//            K[4] = h3*(-Kvect.at(15)*((I1 + 19372*h3*K[0]/6561.0 + (-25360*h3*K[1]/2187.0) + 64448*h3*K[2]/6561.0 + (-212*h3*K[3]/729.0) )-(Ip+ 8*h3/9.0)/Kvect.at(4)));
//            K[5] = h3*(-Kvect.at(15)*((I1 + 9017*h3*K[0]/3168.0 + (-355*h3*K[1]/33.0) + 46732*h3*K[2]/5247.0 + (49*h3*K[3]/176.0) + (-5103*h3*K[4]/18656.0) )-(Ip+ h3)/Kvect.at(4)));
//            K[6] = h3*(-Kvect.at(15)*((I1 + 35*h3*K[0]/384.0 + 500*h3*K[2]/1113.0 + (125*h3*K[3]/192.0) + (-2187*h3*K[4]/6784.0) + 11*h3*K[5]/84.0 )-(Ip+ h3)/Kvect.at(4)));

//            dI1 = I1 + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);     // во многих источниках домножается на h, но в одном нету такого.
//            dz = I1 + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dI1);
//            s = pow(eps*h3/(2*err),1/5);
//            hopt3 = s*h3;
//            if( hopt3 < hmin) hopt3 = hmin;
//            else if(hopt3 > hmax) hopt3 = hmax;

//            // Решение dId через Дормана-Принса

//            K[0] = h4*(-Kvect.at(15)*(Id-I1));
//            K[1] =  h4*(-Kvect.at(15)*((Id + h4*K[0]/5.0)-(I1+h4/5.0)));
//            K[2] =  h4*(-Kvect.at(15)*((Id + 3*h4*K[0]/40.0 + 9*h4*K[1]/40.0)-(I1+3*h4/10.0)));
//            K[3] =  h4*(-Kvect.at(15)*((Id + 44*h4*K[0]/45.0 + (-56*h4*K[1]/15.0) + 32*h4*K[2]/9.0)-(I1+4*h4/5.0)));
//            K[4] =  h4*(-Kvect.at(15)*((Id + 19372*h4*K[0]/6561.0 + (-25360*h4*K[1]/2187.0) + 64448*h4*K[2]/6561.0 + (-212*h4*K[3]/729.0))-(I1+8*h4/9.0)));
//            K[5] =  h4*(-Kvect.at(15)*((Id + 9017*h4*K[0]/3168.0 + (-355*h4*K[1]/33.0) + 46732*h4*K[2]/5247.0 + (49*h4*K[3]/176.0) + (-5103*h4*K[4]/18656.0))-(I1+h4)));
//            K[6] =  h4*(-Kvect.at(15)*((Id + 35*h4*K[0]/384.0 + 500*h4*K[2]/1113.0 + (125*h4*K[3]/192.0) + (-2187*h4*K[4]/6784.0) + 11*h4*K[5]/84.0)-(I1+h4)));

//            dId   = Id + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);                                                // во многих источниках домножается на h, но в одном нету такого.
//            dz = Id + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dId);
//            s = pow(eps*h4/(2*err),1/5);
//            hopt4 = s*h4;
//            if( hopt4 < hmin) hopt4 = hmin;
//            else if(hopt4 > hmax) hopt4 = hmax;

//            // Решение dXt через Дормана-Принса

//            double Heavi2;
//            if( ((Ip/Kvect.at(4))-Kvect.at(1)) >= 0){
//                 Heavi2 = 1;}
//            if( ((Ip/Kvect.at(4))-Kvect.at(1)) < 0){
//                 Heavi2 = 0;}
//            K[0] = h5*(-Kvect.at(30)*Xt+Kvect.at(30)*((Ip/Kvect.at(4))-Kvect.at(1))*Heavi2);

//            if( (((Ip+ h5/5.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
//                 Heavi2 = 1;}
//            if( (((Ip+ h5/5.0)/Kvect.at(4))-Kvect.at(1)) < 0){
//                 Heavi2 = 0;}
//            K[1] = h5*(-Kvect.at(30)*(Xt+ h5*K[0]/5.0)+Kvect.at(30)*(((Ip+ h5/5.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

//            if( (((Ip+ 3*h5/10.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
//                 Heavi2 = 1;}
//            if( (((Ip+ 3*h5/10.0)/Kvect.at(4))-Kvect.at(1)) < 0){
//                 Heavi2 = 0;}
//            K[2] = h5*(-Kvect.at(30)*(Xt+ 3*h5*K[0]/40.0 + 9*h5*K[1]/40.0)+Kvect.at(30)*(((Ip+ 3*h5/10.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

//            if( (((Ip+ 4*h5/5.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
//                 Heavi2 = 1;}
//            if( (((Ip+ 4*h5/5.0)/Kvect.at(4))-Kvect.at(1)) < 0){
//                 Heavi2 = 0;}
//            K[3] = h5*(-Kvect.at(30)*(Xt+ 44*h5*K[0]/45.0 + (-56*h5*K[1]/15.0) + 32*h5*K[2]/9.0)+Kvect.at(30)*(((Ip+ 4*h5/5.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

//            if( (((Ip+ 8*h5/9.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
//                 Heavi2 = 1;}
//            if( (((Ip+ 8*h5/9.0)/Kvect.at(4))-Kvect.at(1)) < 0){
//                 Heavi2 = 0;}
//            K[4] = h5*(-Kvect.at(30)*(Xt+ 19372*h5*K[0]/6561.0 + (-25360*h5*K[1]/2187.0) + 64448*h5*K[2]/6561.0 + (-212*h5*K[3]/729.0) )+Kvect.at(30)*(((Ip+ 8*h5/9.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

//            if( (((Ip+ h5)/Kvect.at(4))-Kvect.at(1)) >= 0){
//                 Heavi2 = 1;}
//            if( (((Ip+ h5)/Kvect.at(4))-Kvect.at(1)) < 0){
//                 Heavi2 = 0;}
//            K[5] = h5*(-Kvect.at(30)*(Xt+ 9017*h5*K[0]/3168.0 + (-355*h5*K[1]/33.0) + 46732*h5*K[2]/5247.0 + (49*h5*K[3]/176.0) + (-5103*h5*K[4]/18656.0) )+Kvect.at(30)*(((Ip+ h5)/Kvect.at(4))-Kvect.at(1))*Heavi2);

//            if( (((Ip+ h5)/Kvect.at(4))-Kvect.at(1)) >= 0){
//                 Heavi2 = 1;}
//            if( (((Ip+ h5)/Kvect.at(4))-Kvect.at(1)) < 0){
//                 Heavi2 = 0;}
//            K[6] = h5*(-Kvect.at(30)*(Xt+ 35*h5*K[0]/384.0 + 500*h5*K[2]/1113.0 + (125*h5*K[3]/192.0) + (-2187*h5*K[4]/6784.0) + 11*h5*K[5]/84.0)+Kvect.at(30)*(((Ip+ h5)/Kvect.at(4))-Kvect.at(1))*Heavi2);

//            dXt   = Xt + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);                                                // во многих источниках домножается на h, но в одном нету такого.
//            dz = Xt + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dXt);
//            s = pow(eps*h5/(2*err),1/5);
//            hopt5 = s*h5;
//            if( hopt5 < hmin) hopt5 = hmin;
//            else if(hopt5 > hmax) hopt5 = hmax;

//            // Решение dIl через Дормана-Принса

//            K[0] = h6*(-(Kvect.at(5)+Kvect.at(7))*Il+Kvect.at(6)*Ip);
//            K[1] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ h6*K[0]/5.0)+Kvect.at(6)*(Ip+ h6/5.0));
//            K[2] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 3*h6*K[0]/40.0 + 9*h6*K[1]/40.0)+Kvect.at(6)*(Ip+ 3*h6/10.0));
//            K[3] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 44*h6*K[0]/45.0 + (-56*h6*K[1]/15.0) + 32*h6*K[2]/9.0)+Kvect.at(6)*(Ip+ 4*h6/5.0));
//            K[4] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 19372*h6*K[0]/6561.0 + (-25360*h6*K[1]/2187.0) + 64448*h6*K[2]/6561.0 + (-212*h6*K[3]/729.0))+Kvect.at(6)*(Ip+ 8*h6/9.0));
//            K[5] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 9017*h6*K[0]/3168.0 + (-355*h6*K[1]/33.0) + 46732*h6*K[2]/5247.0 + (49*h6*K[3]/176.0) + (-5103*h6*K[4]/18656.0))+Kvect.at(6)*(Ip+ h6));
//            K[6] = h6*(-(Kvect.at(5)+Kvect.at(7))*(Il+ 35*h6*K[0]/384.0 + 500*h6*K[2]/1113.0 + (125*h6*K[3]/192.0) + (-2187*h6*K[4]/6784.0) + 11*h6*K[5]/84.0)+Kvect.at(6)*(Ip+ h6));

//            dIl   = Il + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);                                                // во многих источниках домножается на h, но в одном нету такого.
//            dz = Il + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dIl);
//            s = pow(eps*h6/(2*err),1/5);
//            hopt6 = s*h6;
//            if( hopt6 < hmin) hopt6 = hmin;
//            else if(hopt6 > hmax) hopt6 = hmax;

//            // Решение dIp через Дормана-Принса

//            K[0] = h7*(-Kvect.at(6)*Ip+Kvect.at(5)*Il+Kvect.at(10)/Kvect.at(0)*It-Kvect.at(9)*Ip);
//            K[1] = h7*(-Kvect.at(6)*(Ip+ h7*K[0]/5.0)+Kvect.at(5)*(Il+ h7/5.0)+Kvect.at(10)/Kvect.at(0)*(It+ h7/5.0)-Kvect.at(9)*(Ip+ h7*K[0]/5.0));
//            K[2] = h7*(-Kvect.at(6)*(Ip+ 3*h7*K[0]/40.0 + 9*h7*K[1]/40.0)+Kvect.at(5)*(Il+ 3*h7/10.0)+Kvect.at(10)/Kvect.at(0)*(It+ 3*h7/10.0)-Kvect.at(9)*(Ip+ 3*h7*K[0]/40.0 + 9*h7*K[1]/40.0));
//            K[3] = h7*(-Kvect.at(6)*(Ip+ 44*h7*K[0]/45.0 + (-56*h7*K[1]/15.0) + 32*h7*K[2]/9.0)+Kvect.at(5)*(Il+ 4*h7/5.0)+Kvect.at(10)/Kvect.at(0)*(It+ 4*h7/5.0)-Kvect.at(9)*(Ip+ 44*h7*K[0]/45.0 + (-56*h7*K[1]/15.0) + 32*h7*K[2]/9.0));
//            K[4] = h7*(-Kvect.at(6)*(Ip+ 19372*h7*K[0]/6561.0 + (-25360*h7*K[1]/2187.0) + 64448*h7*K[2]/6561.0 + (-212*h7*K[3]/729.0) )+Kvect.at(5)*(Il+ 8*h7/9.0)+Kvect.at(10)/Kvect.at(0)*(It+ 8*h7/9.0)-Kvect.at(9)*(Ip+ 19372*h7*K[0]/6561.0 + (-25360*h7*K[1]/2187.0) + 64448*h7*K[2]/6561.0 + (-212*h7*K[3]/729.0) ));
//            K[5] = h7*(-Kvect.at(6)*(Ip+ 9017*h7*K[0]/3168.0 + (-355*h7*K[1]/33.0) + 46732*h7*K[2]/5247.0 + (49*h7*K[3]/176.0) + (-5103*h7*K[4]/18656.0))+Kvect.at(5)*(Il+ h7)+Kvect.at(10)/Kvect.at(0)*(It+ h7)-Kvect.at(9)*(Ip+ 9017*h7*K[0]/3168.0 + (-355*h7*K[1]/33.0) + 46732*h7*K[2]/5247.0 + (49*h7*K[3]/176.0) + (-5103*h7*K[4]/18656.0)));
//            K[6] = h7*(-Kvect.at(6)*(Ip+ 35*h7*K[0]/384.0 + 500*h7*K[2]/1113.0 + (125*h7*K[3]/192.0) + (-2187*h7*K[4]/6784.0) + 11*h7*K[5]/84.0)+Kvect.at(5)*(Il+ h7)+Kvect.at(10)/Kvect.at(0)*(It+ h7)-Kvect.at(9)*(Ip+ 35*h7*K[0]/384.0 + 500*h7*K[2]/1113.0 + (125*h7*K[3]/192.0) + (-2187*h7*K[4]/6784.0) + 11*h7*K[5]/84.0));

//            dIp   = Ip + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
//            dz = Ip + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dIp);
//            s = pow(eps*h7/(2*err),1/5);
//            hopt7 = s*h7;
//            if( hopt7 < hmin) hopt7 = hmin;
//            else if(hopt7 > hmax) hopt7 = hmax;

//            double kgut;
//            kgut=Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2);

////            // Решение dfgut через Дормана-Принса

////            K[0] = h8*(-Kvect.at(17)*fgut+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2)*fliq);
////            K[1] = h8*(-Kvect.at(17)*(fgut+ h8*K[0]/5.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h8/5.0)+(fliq+ h8/5.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h8/5.0)+(fliq+ h8/5.0)-Kvect.at(22)*Dig))+2)*(fliq+ h8/5.0));
////            K[2] = h8*(-Kvect.at(17)*(fgut+ 3*h8*K[0]/40.0 + 9*h8*K[1]/40.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 3*h8/10.0)+(fliq+ 3*h8/10.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 3*h8/10.0)+(fliq+ 3*h8/10.0)-Kvect.at(22)*Dig))+2)*(fliq+ 3*h8/10.0));
////            K[3] = h8*(-Kvect.at(17)*(fgut+ 44*h8*K[0]/45.0 + (-56*h8*K[1]/15.0) + 32*h8*K[2]/9.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 4*h8/5.0)+(fliq+ 4*h8/5.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 4*h8/5.0)+(fliq+ 4*h8/5.0)-Kvect.at(22)*Dig))+2)*(fliq+ 4*h8/5.0));
////            K[4] = h8*(-Kvect.at(17)*(fgut+ 19372*h8*K[0]/6561.0 + (-25360*h8*K[1]/2187.0) + 64448*h8*K[2]/6561.0 + (-212*h8*K[3]/729.0) )+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 8*h8/9.0)+(fliq+ 8*h8/9.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 8*h8/9.0)+(fliq+ 8*h8/9.0)-Kvect.at(22)*Dig))+2)*(fliq+ 8*h8/9.0));
////            K[5] = h8*(-Kvect.at(17)*(fgut+ 9017*h8*K[0]/3168.0 + (-355*h8*K[1]/33.0) + 46732*h8*K[2]/5247.0 + (49*h8*K[3]/176.0) + (-5103*h8*K[4]/18656.0))+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h8)+(fliq+ h8)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h8)+(fliq+ h8)-Kvect.at(22)*Dig))+2)*(fliq+ h8));
////            K[6] = h8*(-Kvect.at(17)*(fgut+ 35*h8*K[0]/384.0 + 500*h8*K[2]/1113.0 + (125*h8*K[3]/192.0) + (-2187*h8*K[4]/6784.0) + 11*h8*K[5]/84.0)+Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h8)+(fliq+ h8)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h8)+(fliq+ h8)-Kvect.at(22)*Dig))+2)*(fliq+ h8));

//            K[0] = h8*(-Kvect.at(17)*fgut+kgut*fliq);
//            K[1] = h8*(-Kvect.at(17)*(fgut+ h8*K[0]/5.0)+kgut*(fliq+ h8/5.0));
//            K[2] = h8*(-Kvect.at(17)*(fgut+ 3*h8*K[0]/40.0 + 9*h8*K[1]/40.0)+kgut*(fliq+ 3*h8/10.0));
//            K[3] = h8*(-Kvect.at(17)*(fgut+ 44*h8*K[0]/45.0 + (-56*h8*K[1]/15.0) + 32*h8*K[2]/9.0)+kgut*(fliq+ 4*h8/5.0));
//            K[4] = h8*(-Kvect.at(17)*(fgut+ 19372*h8*K[0]/6561.0 + (-25360*h8*K[1]/2187.0) + 64448*h8*K[2]/6561.0 + (-212*h8*K[3]/729.0) )+kgut*(fliq+ 8*h8/9.0));
//            K[5] = h8*(-Kvect.at(17)*(fgut+ 9017*h8*K[0]/3168.0 + (-355*h8*K[1]/33.0) + 46732*h8*K[2]/5247.0 + (49*h8*K[3]/176.0) + (-5103*h8*K[4]/18656.0))+kgut*(fliq+ h8));
//            K[6] = h8*(-Kvect.at(17)*(fgut+ 35*h8*K[0]/384.0 + 500*h8*K[2]/1113.0 + (125*h8*K[3]/192.0) + (-2187*h8*K[4]/6784.0) + 11*h8*K[5]/84.0)+kgut*(fliq+ h8));

//            dfgut   = fgut + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
//            dz = fgut + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dfgut);
//            s = pow(eps*h8/(2*err),1/5);
//            hopt8 = s*h8;
//            if( hopt8 < hmin) hopt8 = hmin;
//            else if(hopt8 > hmax) hopt8 = hmax;

//            // Решение dfliq через Дормана-Принса kgut

//            K[0] = h9*(-kgut*fliq+Kvect.at(18)*fsol);
//            K[1] = h9*(-kgut*(fliq+ h9*K[0]/5.0)+Kvect.at(18)*(fsol+ h9/5.0));
//            K[2] = h9*(-kgut*(fliq+ 3*h9*K[0]/40.0 + 9*h9*K[1]/40.0)+Kvect.at(18)*(fsol+ 3*h9/10.0));
//            K[3] = h9*(-kgut*(fliq+ 44*h9*K[0]/45.0 + (-56*h9*K[1]/15.0) + 32*h9*K[2]/9.0)+Kvect.at(18)*(fsol+ 4*h9/5.0));
//            K[4] = h9*(-kgut*(fliq+ 19372*h9*K[0]/6561.0 + (-25360*h9*K[1]/2187.0) + 64448*h9*K[2]/6561.0 + (-212*h9*K[3]/729.0))+Kvect.at(18)*(fsol+ 8*h9/9.0));
//            K[5] = h9*(-kgut*(fliq+ 9017*h9*K[0]/3168.0 + (-355*h9*K[1]/33.0) + 46732*h9*K[2]/5247.0 + (49*h9*K[3]/176.0) + (-5103*h9*K[4]/18656.0))+Kvect.at(18)*(fsol+ h9));
//            K[6] = h9*(-kgut*(fliq+ 35*h9*K[0]/384.0 + 500*h9*K[2]/1113.0 + (125*h9*K[3]/192.0) + (-2187*h9*K[4]/6784.0) + 11*h9*K[5]/84.0)+Kvect.at(18)*(fsol+ h9));

////            K[0] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(fsol+fliq-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(fsol+fliq-Kvect.at(22)*Dig))+2))*fliq+Kvect.at(18)*fsol);
////            K[1] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h9/5.0)+(fliq+ h9*K[0]/5.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h9/5.0)+(fliq+ h9*K[0]/5.0)-Kvect.at(22)*Dig))+2))*(fliq+ h9*K[0]/5.0)+Kvect.at(18)*(fsol+ h9/5.0));
////            K[2] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 3*h9/10.0)+(fliq+ 3*h9*K[0]/40.0 + 9*h9*K[1]/40.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 3*h9/10.0)+(fliq+3*h9*K[0]/40.0 + 9*h9*K[1]/40.0)-Kvect.at(22)*Dig))+2))*(fliq+ 3*h9*K[0]/40.0 + 9*h9*K[1]/40.0)+Kvect.at(18)*(fsol+ 3*h9/10.0));
////            K[3] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 4*h9/5.0)+(fliq+ 44*h9*K[0]/45.0 + (-56*h9*K[1]/15.0) + 32*h9*K[2]/9.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 4*h9/5.0)+(fliq+ 44*h9*K[0]/45.0 + (-56*h9*K[1]/15.0) + 32*h9*K[2]/9.0)-Kvect.at(22)*Dig))+2))*(fliq+ 44*h9*K[0]/45.0 + (-56*h9*K[1]/15.0) + 32*h9*K[2]/9.0)+Kvect.at(18)*(fsol+ 4*h9/5.0));
////            K[4] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ 8*h9/9.0)+(fliq+ 19372*h9*K[0]/6561.0 + (-25360*h9*K[1]/2187.0) + 64448*h9*K[2]/6561.0 + (-212*h9*K[3]/729.0))-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ 8*h9/9.0)+(fliq+ 19372*h9*K[0]/6561.0 + (-25360*h9*K[1]/2187.0) + 64448*h9*K[2]/6561.0 + (-212*h9*K[3]/729.0))-Kvect.at(22)*Dig))+2))*(fliq+ 19372*h9*K[0]/6561.0 + (-25360*h9*K[1]/2187.0) + 64448*h9*K[2]/6561.0 + (-212*h9*K[3]/729.0))+Kvect.at(18)*(fsol+ 8*h9/9.0));
////            K[5] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h9)+(fliq+ 9017*h9*K[0]/3168.0 + (-355*h9*K[1]/33.0) + 46732*h9*K[2]/5247.0 + (49*h9*K[3]/176.0) + (-5103*h9*K[4]/18656.0))-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h9)+(fliq+ 9017*h9*K[0]/3168.0 + (-355*h9*K[1]/33.0) + 46732*h9*K[2]/5247.0 + (49*h9*K[3]/176.0) + (-5103*h9*K[4]/18656.0))-Kvect.at(22)*Dig))+2))*(fliq+ 9017*h9*K[0]/3168.0 + (-355*h9*K[1]/33.0) + 46732*h9*K[2]/5247.0 + (49*h9*K[3]/176.0) + (-5103*h9*K[4]/18656.0))+Kvect.at(18)*(fsol+ h9));
////            K[6] = h9*(-(Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*((fsol+ h9)+(fliq+ 35*h9*K[0]/384.0 + 500*h9*K[2]/1113.0 + (125*h9*K[3]/192.0) + (-2187*h9*K[4]/6784.0) + 11*h9*K[5]/84.0)-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*((fsol+ h9)+(fliq+ 35*h9*K[0]/384.0 + 500*h9*K[2]/1113.0 + (125*h9*K[3]/192.0) + (-2187*h9*K[4]/6784.0) + 11*h9*K[5]/84.0)-Kvect.at(22)*Dig))+2))*(fliq+ 35*h9*K[0]/384.0 + 500*h9*K[2]/1113.0 + (125*h9*K[3]/192.0) + (-2187*h9*K[4]/6784.0) + 11*h9*K[5]/84.0)+Kvect.at(18)*(fsol+ h9));

//            dfliq   = fliq + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
//            dz = fliq + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dfliq);
//            s = pow(eps*h9/(2*err),1/5);
//            hopt9 = s*h9;
//            if( hopt9 < hmin) hopt9 = hmin;
//            else if(hopt9 > hmax) hopt9 = hmax;

//            // Решение dfsol через Дормана-Принса

//            K[0] = h10*(-Kvect.at(18)*fsol+vm);
//            K[1] = h10*(-Kvect.at(18)*(fsol+ h10*K[0]/5.0)+(vm));
//            K[2] = h10*(-Kvect.at(18)*(fsol+ 3*h10*K[0]/40.0 + 9*h10*K[1]/40.0)+(vm));
//            K[3] = h10*(-Kvect.at(18)*(fsol+ 44*h10*K[0]/45.0 + (-56*h10*K[1]/15.0) + 32*h10*K[2]/9.0)+(vm));
//            K[4] = h10*(-Kvect.at(18)*(fsol+ 19372*h10*K[0]/6561.0 + (-25360*h10*K[1]/2187.0) + 64448*h10*K[2]/6561.0 + (-212*h10*K[3]/729.0))+(vm));
//            K[5] = h10*(-Kvect.at(18)*(fsol+ 9017*h10*K[0]/3168.0 + (-355*h10*K[1]/33.0) + 46732*h10*K[2]/5247.0 + (49*h10*K[3]/176.0) + (-5103*h10*K[4]/18656.0))+(vm));
//            K[6] = h10*(-Kvect.at(18)*(fsol+ 35*h10*K[0]/384.0 + 500*h10*K[2]/1113.0 + (125*h10*K[3]/192.0) + (-2187*h10*K[4]/6784.0) + 11*h10*K[5]/84.0)+(vm));

//            dfsol   = fsol + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
//            dz = fsol + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dfsol);
//            s = pow(eps*h10/(2*err),1/5);
//            hopt10 = s*h10;
//            if( hopt10 < hmin) hopt10 = hmin;
//            else if(hopt10 > hmax) hopt10 = hmax;

//            // Решение dIpo через Дормана-Принса

//            double Heavi5;
//            double Heavi6;

//            if( (dgp) >= 0){
//                 Heavi5 = 1;}
//            if( (dgp) < 0){
//                 Heavi5 = 0;}

//            if( (-dgp) >= 0){
//                 Heavi6 = 1;}
//            if( (-dgp) < 0){
//                 Heavi6 = 0;}
//            K[0] = h11*(-Kvect.at(33)*Ipo+(Yt+Kvect.at(3))*Heavi5+(Yt+Kvect.at(3))*(Heavi6));

//            K[1] = h11*(-Kvect.at(33)*(Ipo+ h11*K[0]/5.0)+((Yt+ h11/5.0)+Kvect.at(3))*Heavi5+((Yt+ h11/5.0)+Kvect.at(3))*(Heavi6));

//            K[2] = h11*(-Kvect.at(33)*(Ipo+ 3*h11*K[0]/40.0 + 9*h11*K[1]/40.0)+((Yt+ 3*h11/10.0)+Kvect.at(3))*Heavi5+((Yt+ 3*h11/10.0)+Kvect.at(3))*(Heavi6));

//            K[3] = h11*(-Kvect.at(33)*(Ipo+ 44*h11*K[0]/45.0 + (-56*h11*K[1]/15.0) + 32*h11*K[2]/9.0)+((Yt+ 4*h11/5.0)+Kvect.at(3))*Heavi5+((Yt+ 4*h11/5.0)+Kvect.at(3))*(Heavi6));

//            K[4] = h11*(-Kvect.at(33)*(Ipo+ 19372*h11*K[0]/6561.0 + (-25360*h11*K[1]/2187.0) + 64448*h11*K[2]/6561.0 + (-212*h11*K[3]/729.0))+((Yt+ 8*h11/9.0)+Kvect.at(3))*Heavi5+((Yt+ 8*h11/9.0)+Kvect.at(3))*(Heavi6));

//            K[5] = h11*(-Kvect.at(33)*(Ipo+ 9017*h11*K[0]/3168.0 + (-355*h11*K[1]/33.0) + 46732*h11*K[2]/5247.0 + (49*h11*K[3]/176.0) + (-5103*h11*K[4]/18656.0))+((Yt+ h11)+Kvect.at(3))*Heavi5+((Yt+ h11)+Kvect.at(3))*(Heavi6));

//            K[6] = h11*(-Kvect.at(33)*(Ipo+ 35*h11*K[0]/384.0 + 500*h11*K[2]/1113.0 + (125*h11*K[3]/192.0) + (-2187*h11*K[4]/6784.0) + 11*h11*K[5]/84.0)+((Yt+ h11)+Kvect.at(3))*Heavi5+((Yt+ h11)+Kvect.at(3))*(Heavi6));

//            dIpo   = Ipo + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
//            dz = Ipo + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dIpo);
//            s = pow(eps*h11/(2*err),1/5);
//            hopt11 = s*h11;
//            if( hopt11 < hmin) hopt11 = hmin;
//            else if(hopt11 > hmax) hopt11 = hmax;

//            // Решение dYt через Дормана-Принса

//            double Heavi7;
//            double Heavi8;

//            if( (Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//                 Heavi7 = 1;}
//            if( (Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//                 Heavi7 = 0;}

//            if( (-Kvect.at(3)-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))) >= 0){
//                 Heavi8 = 1;}
//            if( (-Kvect.at(3)-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2))) < 0){
//                 Heavi8 = 0;}
//            K[0] = h12*(-Kvect.at(34)*(Yt-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*Yt-Kvect.at(34)*Kvect.at(3))*(Heavi8));

//            if( (Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//                 Heavi7 = 1;}
//            if( (Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//                 Heavi7 = 0;}

//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
//                 Heavi8 = 1;}
//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2))) < 0){
//                 Heavi8 = 0;}
//            K[1] = h12*(-Kvect.at(34)*((Yt+ h12*K[0]/5.0)-Kvect.at(35)*((gp+ h12/5.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ h12*K[0]/5.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

//            if( (Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//                 Heavi7 = 1;}
//            if( (Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//                 Heavi7 = 0;}

//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
//                 Heavi8 = 1;}
//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2))) < 0){
//                 Heavi8 = 0;}
//            K[2] = h12*(-Kvect.at(34)*((Yt+ 3*h12*K[0]/40.0 + 9*h12*K[1]/40.0)-Kvect.at(35)*((gp+ 3*h12/10.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 3*h12*K[0]/40.0 + 9*h12*K[1]/40.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

//            if( (Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//                 Heavi7 = 1;}
//            if( (Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//                 Heavi7 = 0;}

//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
//                 Heavi8 = 1;}
//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2))) < 0){
//                 Heavi8 = 0;}
//            K[3] = h12*(-Kvect.at(34)*((Yt+ 44*h12*K[0]/45.0 + (-56*h12*K[1]/15.0) + 32*h12*K[2]/9.0)-Kvect.at(35)*((gp+ 4*h12/5.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 44*h12*K[0]/45.0 + (-56*h12*K[1]/15.0) + 32*h12*K[2]/9.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

//            if( (Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//                 Heavi7 = 1;}
//            if( (Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//                 Heavi7 = 0;}

//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
//                 Heavi8 = 1;}
//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2))) < 0){
//                 Heavi8 = 0;}
//            K[4] = h12*(-Kvect.at(34)*((Yt+ 19372*h12*K[0]/6561.0 + (-25360*h12*K[1]/2187.0) + 64448*h12*K[2]/6561.0 + (-212*h12*K[3]/729.0))-Kvect.at(35)*((gp+ 8*h12/9.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 19372*h12*K[0]/6561.0 + (-25360*h12*K[1]/2187.0) + 64448*h12*K[2]/6561.0 + (-212*h12*K[3]/729.0))-Kvect.at(34)*Kvect.at(3))*(Heavi8));

//            if( (Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//                 Heavi7 = 1;}
//            if( (Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//                 Heavi7 = 0;}

//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))) >= 0){
//                 Heavi8 = 1;}
//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))) < 0){
//                 Heavi8 = 0;}
//            K[5] = h12*(-Kvect.at(34)*((Yt+ 9017*h12*K[0]/3168.0 + (-355*h12*K[1]/33.0) + 46732*h12*K[2]/5247.0 + (49*h12*K[3]/176.0) + (-5103*h12*K[4]/18656.0))-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 9017*h12*K[0]/3168.0 + (-355*h12*K[1]/33.0) + 46732*h12*K[2]/5247.0 + (49*h12*K[3]/176.0) + (-5103*h12*K[4]/18656.0))-Kvect.at(34)*Kvect.at(3))*(Heavi8));

//            if( (Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
//                 Heavi7 = 1;}
//            if( (Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
//                 Heavi7 = 0;}

//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))) >= 0){
//                 Heavi8 = 1;}
//            if( (-Kvect.at(3)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2))) < 0){
//                 Heavi8 = 0;}
//            K[6] = h12*(-Kvect.at(34)*((Yt+ 35*h12*K[0]/384.0 + 500*h12*K[2]/1113.0 + (125*h12*K[3]/192.0) + (-2187*h12*K[4]/6784.0) + 11*h12*K[5]/84.0)-Kvect.at(35)*((gp+ h12)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(Yt+ 35*h12*K[0]/384.0 + 500*h12*K[2]/1113.0 + (125*h12*K[3]/192.0) + (-2187*h12*K[4]/6784.0) + 11*h12*K[5]/84.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

//            dYt   = Yt + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
//            dz = Yt + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dYt);
//            s = pow(eps*h12/(2*err),1/5);
//            hopt12 = s*h12;
//            if( hopt12 < hmin) hopt12 = hmin;
//            else if(hopt12 > hmax) hopt12 = hmax;

//            // Решение dIt через Дормана-Принса

//            K[0] = h13*(Kvect.at(8)*Ii-Kvect.at(10)*It);
//            K[1] = h13*(Kvect.at(8)*(Ii+ h13/5.0)-Kvect.at(10)*(It+ h13*K[0]/5.0));
//            K[2] = h13*(Kvect.at(8)*(Ii+ 3*h13/10.0)-Kvect.at(10)*(It+ 3*h13*K[0]/40.0 + 9*h13*K[1]/40.0));
//            K[3] = h13*(Kvect.at(8)*(Ii+ 4*h13/5.0)-Kvect.at(10)*(It+ 44*h13*K[0]/45.0 + (-56*h13*K[1]/15.0) + 32*h13*K[2]/9.0));
//            K[4] = h13*(Kvect.at(8)*(Ii+ 8*h13/9.0)-Kvect.at(10)*(It+ 19372*h13*K[0]/6561.0 + (-25360*h13*K[1]/2187.0) + 64448*h13*K[2]/6561.0 + (-212*h13*K[3]/729.0)));
//            K[5] = h13*(Kvect.at(8)*(Ii+ h13)-Kvect.at(10)*(It+ 9017*h13*K[0]/3168.0 + (-355*h13*K[1]/33.0) + 46732*h13*K[2]/5247.0 + (49*h13*K[3]/176.0) + (-5103*h13*K[4]/18656.0)));
//            K[6] = h13*(Kvect.at(8)*(Ii+ h13)-Kvect.at(10)*(It+ 35*h13*K[0]/384.0 + 500*h13*K[2]/1113.0 + (125*h13*K[3]/192.0) + (-2187*h13*K[4]/6784.0) + 11*h13*K[5]/84.0));

//            dIt   = It + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
//            dz = It + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dIt);
//            s = pow(eps*h13/(2*err),1/5);
//            hopt13 = s*h13;
//            if( hopt13 < hmin) hopt13 = hmin;
//            else if(hopt13 > hmax) hopt13 = hmax;

//            // Решение dIi через Дормана-Принса

//            K[0] = h14*(-Kvect.at(8)*Ii+vbas+vbol);
//            K[1] = h14*(-Kvect.at(8)*(Ii+ h14*K[0]/5.0)+vbas+vbol);
//            K[2] = h14*(-Kvect.at(8)*(Ii+ 3*h14*K[0]/40.0 + 9*h14*K[1]/40.0)+(vbas)+vbol);
//            K[3] = h14*(-Kvect.at(8)*(Ii+ 44*h14*K[0]/45.0 + (-56*h14*K[1]/15.0) + 32*h14*K[2]/9.0)+(vbas)+vbol);
//            K[4] = h14*(-Kvect.at(8)*(Ii+ 19372*h14*K[0]/6561.0 + (-25360*h14*K[1]/2187.0) + 64448*h14*K[2]/6561.0 + (-212*h14*K[3]/729.0))+(vbas)+(vbol));
//            K[5] = h14*(-Kvect.at(8)*(Ii+ 9017*h14*K[0]/3168.0 + (-355*h14*K[1]/33.0) + 46732*h14*K[2]/5247.0 + (49*h14*K[3]/176.0) + (-5103*h14*K[4]/18656.0))+(vbas)+(vbol));
//            K[6] = h14*(-Kvect.at(8)*(Ii+ 35*h14*K[0]/384.0 + 500*h14*K[2]/1113.0 + (125*h14*K[3]/192.0) + (-2187*h14*K[4]/6784.0) + 11*h14*K[5]/84.0)+(vbas)+(vbol));

//            dIi   = Ii + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
//            dz = Ii + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
//            err = abs(dz-dIi);
//            s = pow(eps*h14/(2*err),1/5);
//            hopt14 = s*h14;
//            if( hopt14 < hmin) hopt14 = hmin;
//            else if(hopt14 > hmax) hopt14 = hmax;

//            /* __________________________________________________________________________ Обработка выходных значений ___________________________________________________________________________ */


//            tick.append(t);
//            CGB.append(dgp/Kvect.at(23));
//            Ipg.append(dIp*20);

//            I1 = dI1;
//            Id = dId;
//            Ii = dIi;
//            Il = dIl;
//            Ip = dIp;
//            It = dIt;
//            Xt = dXt;
//            Yt = dYt;
//            gp = dgp;
//            gt = dgt;
//            Ipo = dIpo;
//            fgut = dfgut;
//            fliq = dfliq;
//            fsol = dfsol;

//            h1 = hopt1;
//            h2 = hopt2;
//            h3 = hopt3;
//            h4 = hopt4;
//            h5 = hopt5;
//            h6 = hopt6;
//            h7 = hopt7;
//            h8 = hopt8;
//            h9 = hopt9;
//            h10 = hopt10;
//            h11 = hopt11;
//            h12 = hopt12;
//            h13 = hopt13;
//            h14 = hopt14;

//            t = t + h; // увеличиваем время
//            j = j + h;
//        } // конец цикла 5 min


        calcVbasRK(t0global);
        //calcVbasDP(t0global);
        tickBas.append(t0global);
        //tickBasDP.append(t0global);
        t0global+= 5; // Увеличиваем глобальое время от 0 до 720

    } // конец глобльного цикла (720)

    GraphOde45();
    GraphRungeKutta();
    //GraphDormanPrince();
}
/* __________________________________________________________________________ Расчёт Vbas ___________________________________________________________________________ */
void MainWindow::calcVbasRK(double t0global)
{
    bas[5] = bas[0]*(1+bas[4]);
    double t1 = 30;
    double t2 = 30 +120;
    if ( (t0global<t1) or (t0global>t2) ){

        double Heavi10;
        double Heavi11;
        double Heavi12;
        double Heavi13;
        double Heavi14;
        double gpVbas = dZ[0]/1.8; //gp/1.8;

        if( (gpVbas-bas[1]) > 0){
             Heavi10 = 1;}
        if( (gpVbas-bas[1]) < 0){
             Heavi10 = 0;}
        if( (gpVbas-bas[1]) == 0){
             Heavi10 = 0.5;}

        if( (bas[6]-gpVbas) > 0){
             Heavi11 = 1;}
        if( (bas[6]-gpVbas) < 0){
             Heavi11 = 0;}
        if( (bas[6]-gpVbas) == 0){
             Heavi11 = 0.5;}

        if( ((bas[2]-gpVbas)) > 0){
             Heavi12 = 1;}
        if( ((bas[2]-gpVbas)) < 0){
             Heavi12 = 0;}
        if( ((bas[2]-gpVbas)) == 0){
             Heavi12 = 0.5;}

        if( (gpVbas-bas[3]) > 0){
             Heavi13 = 1;}
        if( (gpVbas-bas[3]) < 0){
             Heavi13 = 0;}
        if( (gpVbas-bas[3]) == 0){
             Heavi13 = 0.5;}

        if( (gpVbas-bas[6]) > 0){
             Heavi14 = 1;}
        if( (gpVbas-bas[6]) < 0){
             Heavi14 = 0;}
        if( (gpVbas-bas[6]) == 0){
             Heavi14 = 0.5;}

        Vbas = (bas[5]-bas[0])/(bas[6]-bas[1])*(gpVbas-bas[1])*Heavi10*Heavi11-bas[0]/(bas[2]-bas[3])*(bas[2]-gpVbas)*Heavi12*Heavi13+bas[0]*Heavi13+(bas[5]-bas[0])*Heavi14;
        bas[0] = Vbas;
    }
    else{

        Vbas = Vb;
        bas[0] = Vbas;
    }

    grafVbas.append(Vbas);
}
void MainWindow::calcVbasDP(double t0global)
{

    basDP[5] = basDP[0]*(1+basDP[4]);
    double t1 = 30;
    double t2 = 30 +120;
    if ( (t0global<t1) or (t0global>t2) ){

        double Heavi10;
        double Heavi11;
        double Heavi12;
        double Heavi13;
        double Heavi14;
        double gpVbas = DdZ[0]/1.8; //gp/1.8;

        if( (gpVbas-basDP[1]) >= 0){
             Heavi10 = 1;}
        if( (gpVbas-basDP[1]) < 0){
             Heavi10 = 0;}
        if( (basDP[6]-gpVbas) >= 0){
             Heavi11 = 1;}
        if( (basDP[6]-gpVbas) < 0){
             Heavi11 = 0;}
        if( ((basDP[2]-gpVbas)) >= 0){
             Heavi12 = 1;}
        if( ((basDP[2]-gpVbas)) < 0){
             Heavi12 = 0;}
        if( (gpVbas-basDP[3]) >= 0){
             Heavi13 = 1;}
        if( (gpVbas-basDP[3]) < 0){
             Heavi13 = 0;}
        if( (gpVbas-basDP[6]) >= 0){
             Heavi14 = 1;}
        if( (gpVbas-basDP[6]) < 0){
             Heavi14 = 0;}
        VbasDP = (basDP[5]-basDP[0])/(basDP[6]-basDP[1])*(gpVbas-basDP[1])*Heavi10*Heavi11-basDP[0]/(basDP[2]-basDP[3])*(basDP[2]-gpVbas)*Heavi12*Heavi13+basDP[0]*Heavi13+(basDP[5]-basDP[0])*Heavi14;
        basDP[0] = VbasDP;
    }
    else{

        VbasDP = VbDP;
        basDP[0] = Vbas;
    }

    grafVbasDP.append(VbasDP);
}
/* ________________________________________________________ Графики  ______________________________________________________________ */

void MainWindow::GraphRungeKutta()
{

    ui->customPlot->addGraph();
    QPen pen0(Qt::green);
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
    ui->customPlot->graph(2)->setPen(pen4);                                 // Устанавливаем цвет графика
    ui->customPlot->graph(2)->setLineStyle(QCPGraph::LineStyle::lsLine);    // График в виде чего-то там стиль

    ui->customPlot->addGraph();
    ui->customPlot->graph(3)->setPen(QPen(Qt::red, 1, Qt::SolidLine));
    ui->customPlot->graph(3)->addData(0 ,180);
    ui->customPlot->graph(3)->addData(720 ,180);

    ui->customPlot->graph(0)->setData(tick ,CGB);
    ui->customPlot->graph(1)->setData(tick ,Ipg);
    ui->customPlot->graph(2)->setData(tickBas,grafVbas);
    ui->customPlot->rescaleAxes();
    ui->customPlot->replot();
}

void MainWindow::GraphOde45()
{

    ui->matCustomPlot->addGraph();
    ui->matCustomPlot->xAxis->setAutoTickStep(true);
    QPen pen0(Qt::green);
    pen0.setWidth(3);
    ui->matCustomPlot->graph(0)->setPen(pen0);
    ui->matCustomPlot->graph(0)->setLineStyle(QCPGraph::lsLine);               // График в виде чего-то там стиль линий

    ui->matCustomPlot->setInteraction(QCP::iRangeDrag, true);                 // Отключаем взаимодействие перетаскивания графика
    ui->matCustomPlot->setInteraction(QCP::iRangeZoom,true);
    ui->matCustomPlot->addGraph();
    QPen pen3(Qt::blue);
    pen3.setWidth(5);
    ui->matCustomPlot->graph(1)->setPen(pen3);                                 // Устанавливаем цвет графика
    ui->matCustomPlot->graph(1)->setLineStyle(QCPGraph::LineStyle::lsLine);    // График в виде чего-то там стиль

    ui->matCustomPlot->addGraph();
    ui->matCustomPlot->graph(2)->setPen(QPen(Qt::red, 1, Qt::SolidLine));
    ui->matCustomPlot->graph(2)->addData(0 ,180);
    ui->matCustomPlot->graph(2)->addData(720 ,180);

    QVector <double> T720;
    for (int var = 0; var <= 720;) {
        T720.append(var);
        var= var +1;
    }
    ui->matCustomPlot->graph(0)->setData(T720 ,cgbMatlab);
    ui->matCustomPlot->graph(1)->setData(T720 ,insMatlab);

    ui->matCustomPlot->rescaleAxes();
    ui->matCustomPlot->replot();

}
void MainWindow::GraphDormanPrince()
{
    ui->dpCustomPlot->addGraph();
    QPen pen0(Qt::green);
    pen0.setWidth(3);
    ui->dpCustomPlot->xAxis->setAutoTickStep(true);
    ui->dpCustomPlot->graph(0)->setPen(pen0);
    ui->dpCustomPlot->graph(0)->setLineStyle(QCPGraph::lsLine);               // График в виде чего-то там стиль линий

    ui->dpCustomPlot->setInteraction(QCP::iRangeDrag, true);                 // Отключаем взаимодействие перетаскивания графика
    ui->dpCustomPlot->setInteraction(QCP::iRangeZoom,true);
    ui->dpCustomPlot->addGraph();
    QPen pen3(Qt::blue);
    pen3.setWidth(5);                                                       // установить нужную толщину
    ui->dpCustomPlot->graph(1)->setPen(pen3);                                 // Устанавливаем цвет графика
    ui->dpCustomPlot->graph(1)->setLineStyle(QCPGraph::LineStyle::lsLine);    // График в виде чего-то там стиль

    ui->dpCustomPlot->addGraph();
    QPen pen4(Qt::red);
    ui->dpCustomPlot->graph(2)->setPen(pen4);                                 // Устанавливаем цвет графика
    ui->dpCustomPlot->graph(2)->setLineStyle(QCPGraph::LineStyle::lsLine);    // График в виде чего-то там стиль

    ui->dpCustomPlot->addGraph();
    ui->dpCustomPlot->graph(3)->setPen(QPen(Qt::red, 1, Qt::SolidLine));
    ui->dpCustomPlot->graph(3)->addData(0 ,180);
    ui->dpCustomPlot->graph(3)->addData(720 ,180);

    ui->dpCustomPlot->graph(0)->setData(dptick ,dpCGB);
    ui->dpCustomPlot->graph(1)->setData(dptick ,dpIpg);
    ui->dpCustomPlot->graph(2)->setData(tickBasDP,grafVbasDP);
    ui->dpCustomPlot->rescaleAxes();
    ui->dpCustomPlot->replot();
}
void MainWindow::DormanPrince(int *i)
{
    double  a = *i*5 - 5;    // промедуток 5 минут
    double  b = *i*5;        // конечное время

    double eps=0.000001; //error allowance in one step calculation.
    double hmin = 0.01;
    double hmax = 1;

    double h = 0.25;

    double OB = 6.63;//(m2/19) *1.5;      // result of bolus calculation Потом переместить в Глобальный цикл

    double K[7];
    double dz;
    double err;
    double s;

    double del = 0.4;       // result of bolus calculation +
    double m2 = 90;         // + М ?

    double Dig = 1176*m2+1; // вместо m2 массу углеводов
    double Dbol1=del*OB;
    double Dbol2=(1-del)*OB;
    double Dbol3=0.0;

    /* Индивидуальные шаги для Дормана-Принса */
    double hi[14];
    hi[0] = h;
    hi[1] = h;
    hi[2] = h;
    hi[3] = h;
    hi[4] = h;
    hi[5] = h;
    hi[6] = h;
    hi[7] = h;
    hi[8] = h;
    hi[9] = h;
    hi[10] = h;
    hi[11] = h;
    hi[12] = h;
    hi[13] = h;


    /* Индивидуальные оптимизированные шаги для Дормана-Принса */
    double hopt[14];

    for (double j = a; j <= b;)

    {
            /* bolus */
        double Ti1 = 10;
        double Ti2 = 10;
        double ti2 = 60;
        double ti1 = 30;
        double Ti3 = 10;
        double ti3 = 10;
        double tm1 = 60; // tm  // время приёма пищи
        double Tm = 20;         // длительность приёма пищи

        // выходные значения dZout
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

        // Fluctuations
        double vbas=VbasDP*100;      // pmol/min
        double bol1=1/Ti1*Dbol1*(1./(1+exp(-3*(tDP+10-ti1))))*(1./(1+exp(-3*(-10+ti1-tDP+Ti1))));
        double bol2=1/Ti2*Dbol2*(1./(1+exp(-3*(tDP+10-ti2))))*(1./(1+exp(-3*(-10+ti2-tDP+Ti2))));
        double bol3=1/Ti3*Dbol3*(1./(1+exp(-3*(tDP+10-ti3))))*(1./(1+exp(-3*(-10+ti3-tDP+Ti3))));

        double vbol=6000*(bol1+bol2+bol3);

        double vm=Dig/Tm*(1./(1+exp(-3*(tDP-tm1))))*(1./(1+exp(-3*(-(tDP-tm1-Tm)))));

        /* ДУ */

        // Решение dgp через Дормана-Принса

        double Heavi1;
        double Heavi3;
        if( (DdZ[0]-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( (DdZ[0]-Kvect.at(32)) < 0){
             Heavi1 = 0;}

        if( (Kvect.at(11)-Kvect.at(12)*DdZ[0]-Kvect.at(13)*DdZ[1]-Kvect.at(14)*DdZ[2])  >= 0){
             Heavi3 = 1;}
        if( (Kvect.at(11)-Kvect.at(12)*DdZ[0]-Kvect.at(13)*DdZ[1]-Kvect.at(14)*DdZ[2]) < 0){
             Heavi3 = 0;}
        double EGP;
        EGP=(Kvect.at(11)-Kvect.at(12)*DdZ[0]-Kvect.at(13)*DdZ[1]-Kvect.at(14)*DdZ[2])*Heavi3;

        K[0] = hi[0]*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*DdZ[3]-Kvect.at(26)-Kvect.at(31)*(DdZ[0]-Kvect.at(32))*Heavi1-Kvect.at(24)*DdZ[0]+Kvect.at(25)*DdZ[6]);

        if( ((DdZ[0]+ hi[0]*K[0]/5.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((DdZ[0]+ hi[0]*K[0]/5.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        K[1] = hi[0]*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(DdZ[3]+ hi[0]/5.0)-Kvect.at(26)-Kvect.at(31)*((DdZ[0]+ hi[0]*K[0]/5.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(DdZ[0]+ hi[0]*K[0]/5.0)+Kvect.at(25)*(DdZ[6]+ hi[0]/5.0));

        if( ((DdZ[0]+ 3*hi[0]*K[0]/40.0 + 9*hi[0]*K[1]/40.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((DdZ[0]+ 3*hi[0]*K[0]/40.0 + 9*hi[0]*K[1]/40.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        K[2] = hi[0]*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(DdZ[3]+ 3*hi[0]/10.0)-Kvect.at(26)-Kvect.at(31)*((DdZ[0]+ 3*hi[0]*K[0]/40.0 + 9*hi[0]*K[1]/40.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(DdZ[0]+ 3*hi[0]*K[0]/40.0 + 9*hi[0]*K[1]/40.0)+Kvect.at(25)*(DdZ[6]+ 3*hi[0]/10.0));

        if( ((DdZ[0]+ 44*hi[0]*K[0]/45.0 + (-56*hi[0]*K[1]/15.0) + 32*hi[0]*K[2]/9.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((DdZ[0]+ 44*hi[0]*K[0]/45.0 + (-56*hi[0]*K[1]/15.0) + 32*hi[0]*K[2]/9.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        K[3] = hi[0]*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(DdZ[3]+ 4*hi[0]/5.0)-Kvect.at(26)-Kvect.at(31)*((DdZ[0]+ 44*hi[0]*K[0]/45.0 + (-56*hi[0]*K[1]/15.0) + 32*hi[0]*K[2]/9.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(DdZ[0]+ 44*hi[0]*K[0]/45.0 + (-56*hi[0]*K[1]/15.0) + 32*hi[0]*K[2]/9.0)+Kvect.at(25)*(DdZ[6]+ 4*hi[0]/5.0));

        if( ((DdZ[0]+ 19372*hi[0]*K[0]/6561.0 + (-25360*hi[0]*K[1]/2187.0) + 64448*hi[0]*K[2]/6561.0 + (-212*hi[0]*K[3]/729.0) )-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((DdZ[0]+ 19372*hi[0]*K[0]/6561.0 + (-25360*hi[0]*K[1]/2187.0) + 64448*hi[0]*K[2]/6561.0 + (-212*hi[0]*K[3]/729.0) )-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        K[4] = hi[0]*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(DdZ[3]+ 8*hi[0]/9.0)-Kvect.at(26)-Kvect.at(31)*((DdZ[0]+ 19372*hi[0]*K[0]/6561.0 + (-25360*hi[0]*K[1]/2187.0) + 64448*hi[0]*K[2]/6561.0 + (-212*hi[0]*K[3]/729.0) )-Kvect.at(32))*Heavi1-Kvect.at(24)*(DdZ[0]+ 19372*hi[0]*K[0]/6561.0 + (-25360*hi[0]*K[1]/2187.0) + 64448*hi[0]*K[2]/6561.0 + (-212*hi[0]*K[3]/729.0) )+Kvect.at(25)*(DdZ[6]+ 8*hi[0]/9.0));

        if( ((DdZ[0]+ 9017*hi[0]*K[0]/3168.0 + (-355*hi[0]*K[1]/33.0) + 46732*hi[0]*K[2]/5247.0 + (49*hi[0]*K[3]/176.0) + (-5103*hi[0]*K[4]/18656.0))-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((DdZ[0]+ 9017*hi[0]*K[0]/3168.0 + (-355*hi[0]*K[1]/33.0) + 46732*hi[0]*K[2]/5247.0 + (49*hi[0]*K[3]/176.0) + (-5103*hi[0]*K[4]/18656.0))-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        K[5] = hi[0]*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(DdZ[3]+ hi[0])-Kvect.at(26)-Kvect.at(31)*((DdZ[0]+ 9017*hi[0]*K[0]/3168.0 + (-355*hi[0]*K[1]/33.0) + 46732*hi[0]*K[2]/5247.0 + (49*hi[0]*K[3]/176.0) + (-5103*hi[0]*K[4]/18656.0))-Kvect.at(32))*Heavi1-Kvect.at(24)*(DdZ[0]+ 9017*hi[0]*K[0]/3168.0 + (-355*hi[0]*K[1]/33.0) + 46732*hi[0]*K[2]/5247.0 + (49*hi[0]*K[3]/176.0) + (-5103*hi[0]*K[4]/18656.0))+Kvect.at(25)*(DdZ[6]+ hi[0]));

        if( ((DdZ[0]+ 35*hi[0]*K[0]/384.0 + 500*hi[0]*K[2]/1113.0 + (125*hi[0]*K[3]/192.0) + (-2187*hi[0]*K[4]/6784.0) + 11*hi[0]*K[5]/84.0)-Kvect.at(32))>= 0){
             Heavi1 = 1;}
        if( ((DdZ[0]+ 35*hi[0]*K[0]/384.0 + 500*hi[0]*K[2]/1113.0 + (125*hi[0]*K[3]/192.0) + (-2187*hi[0]*K[4]/6784.0) + 11*hi[0]*K[5]/84.0)-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        K[6] = hi[0]*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(DdZ[3]+ hi[0])-Kvect.at(26)-Kvect.at(31)*((DdZ[0]+ 35*hi[0]*K[0]/384.0 + 500*hi[0]*K[2]/1113.0 + (125*hi[0]*K[3]/192.0) + (-2187*hi[0]*K[4]/6784.0) + 11*hi[0]*K[5]/84.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(DdZ[0]+ 35*hi[0]*K[0]/384.0 + 500*hi[0]*K[2]/1113.0 + (125*hi[0]*K[3]/192.0) + (-2187*hi[0]*K[4]/6784.0) + 11*hi[0]*K[5]/84.0)+Kvect.at(25)*(DdZ[6]+ hi[0]));

        dgp   = DdZ[0] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);     // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[0] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dgp);
        s = pow(eps*hi[0]/(2*err),1/5);
        hopt[0] = s*hi[0];
        if( hopt[0] < hmin)  hopt[0] = hmin;
        else if(hopt[0] > hmax) hopt[0] = hmax;

        // Решение dgt через Дормана-Принса

        K[0] = hi[1]*(-((Kvect.at(27)+Kvect.at(28)*DdZ[7])*DdZ[6])/(Kvect.at(29)+DdZ[6])+Kvect.at(24)*DdZ[0]-Kvect.at(25)*DdZ[6]);
        K[1] = hi[1]*(-((Kvect.at(27)+Kvect.at(28)*(DdZ[7] + hi[1]/5.0))*(DdZ[6] + hi[1]*K[0]/5.0))/(Kvect.at(29)+(DdZ[6] + hi[1]*K[0]/5.0))+Kvect.at(24)*(DdZ[0] + hi[1]/5.0)-Kvect.at(25)*(DdZ[6] + hi[1]*K[0]/5.0));
        K[2] = hi[1]*(-((Kvect.at(27)+Kvect.at(28)*(DdZ[7] + 3*hi[1]/10.0))*(DdZ[6] + 3*hi[1]*K[0]/40.0 + 9*hi[1]*K[1]/40.0))/(Kvect.at(29)+(DdZ[6] + 3*hi[1]*K[0]/40.0 + 9*hi[1]*K[1]/40.0))+Kvect.at(24)*(DdZ[0] + 3*hi[1]/10.0)-Kvect.at(25)*(DdZ[6] + 3*hi[1]*K[0]/40.0 + 9*hi[1]*K[1]/40.0));
        K[3] = hi[1]*(-((Kvect.at(27)+Kvect.at(28)*(DdZ[7] + 4*hi[1]/5.0))*(DdZ[6] + 44*hi[1]*K[0]/45.0 + (-56*hi[1]*K[1]/15.0) + 32*hi[1]*K[2]/9.0))/(Kvect.at(29)+(DdZ[6] + 44*hi[1]*K[0]/45.0 + (-56*hi[1]*K[1]/15.0) + 32*hi[1]*K[2]/9.0))+Kvect.at(24)*(DdZ[0] + 4*hi[1]/5.0)-Kvect.at(25)*(DdZ[6] + 44*hi[1]*K[0]/45.0 + (-56*hi[1]*K[1]/15.0) + 32*hi[1]*K[2]/9.0));
        K[4] = hi[1]*(-((Kvect.at(27)+Kvect.at(28)*(DdZ[7] + 8*hi[1]/9.0))*(DdZ[6] + 19372*hi[1]*K[0]/6561.0 + (-25360*hi[1]*K[1]/2187.0) + 64448*hi[1]*K[2]/6561.0 + (-212*hi[1]*K[3]/729.0) ))/(Kvect.at(29)+(DdZ[6] + 19372*hi[1]*K[0]/6561.0 + (-25360*hi[1]*K[1]/2187.0) + 64448*hi[1]*K[2]/6561.0 + (-212*hi[1]*K[3]/729.0) ))+Kvect.at(24)*(DdZ[0] + 8*hi[1]/9.0)-Kvect.at(25)*(DdZ[6] + 19372*hi[1]*K[0]/6561.0 + (-25360*hi[1]*K[1]/2187.0) + 64448*hi[1]*K[2]/6561.0 + (-212*hi[1]*K[3]/729.0) ));
        K[5] = hi[1]*(-((Kvect.at(27)+Kvect.at(28)*(DdZ[7] + hi[1]))*(DdZ[6] + 9017*hi[1]*K[0]/3168.0 + (-355*hi[1]*K[1]/33.0) + 46732*hi[1]*K[2]/5247.0 + (49*hi[1]*K[3]/176.0) + (-5103*hi[1]*K[4]/18656.0) ))/(Kvect.at(29)+(DdZ[6] + 9017*hi[1]*K[0]/3168.0 + (-355*hi[1]*K[1]/33.0) + 46732*hi[1]*K[2]/5247.0 + (49*hi[1]*K[3]/176.0) + (-5103*hi[1]*K[4]/18656.0) ))+Kvect.at(24)*(DdZ[0] + hi[1])-Kvect.at(25)*(DdZ[6] + 9017*hi[1]*K[0]/3168.0 + (-355*hi[1]*K[1]/33.0) + 46732*hi[1]*K[2]/5247.0 + (49*hi[1]*K[3]/176.0) + (-5103*hi[1]*K[4]/18656.0) ));
        K[6] = hi[1]*(-((Kvect.at(27)+Kvect.at(28)*(DdZ[7] + hi[1]))*(DdZ[6] + 35*hi[1]*K[0]/384.0 + 500*hi[1]*K[2]/1113.0 + (125*hi[1]*K[3]/192.0) + (-2187*hi[1]*K[4]/6784.0) + 11*hi[1]*K[5]/84.0))/(Kvect.at(29)+(DdZ[6] + 35*hi[1]*K[0]/384.0 + 500*hi[1]*K[2]/1113.0 + (125*hi[1]*K[3]/192.0) + (-2187*hi[1]*K[4]/6784.0) + 11*hi[1]*K[5]/84.0))+Kvect.at(24)*(DdZ[0] + hi[1])-Kvect.at(25)*(DdZ[6] + 35*hi[1]*K[0]/384.0 + 500*hi[1]*K[2]/1113.0 + (125*hi[1]*K[3]/192.0) + (-2187*hi[1]*K[4]/6784.0) + 11*hi[1]*K[5]/84.0));

        dgt = DdZ[6] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);     // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[6] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dgt);
        s = pow(eps*hi[1]/(2*err),1/5);
        hopt[1] = s*hi[1];
        if( hopt[1] < hmin) hopt[1] = hmin;
        else if(hopt[1] > hmax) hopt[1] = hmax;

        // Решение dI1 через Дормана-Принса

        K[0] = hi[2]*(-Kvect.at(15)*(DdZ[8]-DdZ[9]/Kvect.at(4)));
        K[1] = hi[2]*(-Kvect.at(15)*((DdZ[8] + hi[2]*K[0]/5.0)-(DdZ[9]+ hi[2]/5.0)/Kvect.at(4)));
        K[2] = hi[2]*(-Kvect.at(15)*((DdZ[8] + 3*hi[2]*K[0]/40.0 + 9*hi[2]*K[1]/40.0)-(DdZ[9]+ 3*hi[2]/10.0)/Kvect.at(4)));
        K[3] = hi[2]*(-Kvect.at(15)*((DdZ[8] + 44*hi[2]*K[0]/45.0 + (-56*hi[2]*K[1]/15.0) + 32*hi[2]*K[2]/9.0)-(DdZ[9]+ 4*hi[2]/5.0)/Kvect.at(4)));
        K[4] = hi[2]*(-Kvect.at(15)*((DdZ[8] + 19372*hi[2]*K[0]/6561.0 + (-25360*hi[2]*K[1]/2187.0) + 64448*hi[2]*K[2]/6561.0 + (-212*hi[2]*K[3]/729.0) )-(DdZ[9]+ 8*hi[2]/9.0)/Kvect.at(4)));
        K[5] = hi[2]*(-Kvect.at(15)*((DdZ[8] + 9017*hi[2]*K[0]/3168.0 + (-355*hi[2]*K[1]/33.0) + 46732*hi[2]*K[2]/5247.0 + (49*hi[2]*K[3]/176.0) + (-5103*hi[2]*K[4]/18656.0) )-(DdZ[9]+ hi[2])/Kvect.at(4)));
        K[6] = hi[2]*(-Kvect.at(15)*((DdZ[8] + 35*hi[2]*K[0]/384.0 + 500*hi[2]*K[2]/1113.0 + (125*hi[2]*K[3]/192.0) + (-2187*hi[2]*K[4]/6784.0) + 11*hi[2]*K[5]/84.0 )-(DdZ[9]+ hi[2])/Kvect.at(4)));

        dI1 = DdZ[8] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);     // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[8] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dI1);
        s = pow(eps*hi[2]/(2*err),1/5);
        hopt[2] = s*hi[2];
        if( hopt[2] < hmin) hopt[2] = hmin;
        else if(hopt[2] > hmax) hopt[2] = hmax;

        // Решение dId через Дормана-Принса

        K[0] = hi[3]*(-Kvect.at(15)*(DdZ[1]-DdZ[8]));
        K[1] =  hi[3]*(-Kvect.at(15)*((DdZ[1] + hi[3]*K[0]/5.0)-(DdZ[8]+hi[3]/5.0)));
        K[2] =  hi[3]*(-Kvect.at(15)*((DdZ[1] + 3*hi[3]*K[0]/40.0 + 9*hi[3]*K[1]/40.0)-(DdZ[8]+3*hi[3]/10.0)));
        K[3] =  hi[3]*(-Kvect.at(15)*((DdZ[1] + 44*hi[3]*K[0]/45.0 + (-56*hi[3]*K[1]/15.0) + 32*hi[3]*K[2]/9.0)-(DdZ[8]+4*hi[3]/5.0)));
        K[4] =  hi[3]*(-Kvect.at(15)*((DdZ[1] + 19372*hi[3]*K[0]/6561.0 + (-25360*hi[3]*K[1]/2187.0) + 64448*hi[3]*K[2]/6561.0 + (-212*hi[3]*K[3]/729.0))-(DdZ[8]+8*hi[3]/9.0)));
        K[5] =  hi[3]*(-Kvect.at(15)*((DdZ[1] + 9017*hi[3]*K[0]/3168.0 + (-355*hi[3]*K[1]/33.0) + 46732*hi[3]*K[2]/5247.0 + (49*hi[3]*K[3]/176.0) + (-5103*hi[3]*K[4]/18656.0))-(DdZ[8]+hi[3])));
        K[6] =  hi[3]*(-Kvect.at(15)*((DdZ[1] + 35*hi[3]*K[0]/384.0 + 500*hi[3]*K[2]/1113.0 + (125*hi[3]*K[3]/192.0) + (-2187*hi[3]*K[4]/6784.0) + 11*hi[3]*K[5]/84.0)-(DdZ[8]+hi[3])));

        dId   = DdZ[1] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);                                                // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[1] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dId);
        s = pow(eps*hi[3]/(2*err),1/5);
        hopt[3] = s*hi[3];
        if( hopt[3] < hmin) hopt[3] = hmin;
        else if(hopt[3] > hmax) hopt[3] = hmax;

        // Решение dXt через Дормана-Принса

        double Heavi2;
        if( ((DdZ[9]/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( ((DdZ[9]/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[0] = hi[4]*(-Kvect.at(30)*DdZ[7]+Kvect.at(30)*((DdZ[9]/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((DdZ[9]+ hi[4]/5.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((DdZ[9]+ hi[4]/5.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[1] = hi[4]*(-Kvect.at(30)*(DdZ[7]+ hi[4]*K[0]/5.0)+Kvect.at(30)*(((DdZ[9]+ hi[4]/5.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((DdZ[9]+ 3*hi[4]/10.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((DdZ[9]+ 3*hi[4]/10.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[2] = hi[4]*(-Kvect.at(30)*(DdZ[7]+ 3*hi[4]*K[0]/40.0 + 9*hi[4]*K[1]/40.0)+Kvect.at(30)*(((DdZ[9]+ 3*hi[4]/10.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((DdZ[9]+ 4*hi[4]/5.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((DdZ[9]+ 4*hi[4]/5.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[3] = hi[4]*(-Kvect.at(30)*(DdZ[7]+ 44*hi[4]*K[0]/45.0 + (-56*hi[4]*K[1]/15.0) + 32*hi[4]*K[2]/9.0)+Kvect.at(30)*(((DdZ[9]+ 4*hi[4]/5.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((DdZ[9]+ 8*hi[4]/9.0)/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((DdZ[9]+ 8*hi[4]/9.0)/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[4] = hi[4]*(-Kvect.at(30)*(DdZ[7]+ 19372*hi[4]*K[0]/6561.0 + (-25360*hi[4]*K[1]/2187.0) + 64448*hi[4]*K[2]/6561.0 + (-212*hi[4]*K[3]/729.0) )+Kvect.at(30)*(((DdZ[9]+ 8*hi[4]/9.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((DdZ[9]+ hi[4])/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((DdZ[9]+ hi[4])/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[5] = hi[4]*(-Kvect.at(30)*(DdZ[7]+ 9017*hi[4]*K[0]/3168.0 + (-355*hi[4]*K[1]/33.0) + 46732*hi[4]*K[2]/5247.0 + (49*hi[4]*K[3]/176.0) + (-5103*hi[4]*K[4]/18656.0) )+Kvect.at(30)*(((DdZ[9]+ hi[4])/Kvect.at(4))-Kvect.at(1))*Heavi2);

        if( (((DdZ[9]+ hi[4])/Kvect.at(4))-Kvect.at(1)) >= 0){
             Heavi2 = 1;}
        if( (((DdZ[9]+ hi[4])/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        K[6] = hi[4]*(-Kvect.at(30)*(DdZ[7]+ 35*hi[4]*K[0]/384.0 + 500*hi[4]*K[2]/1113.0 + (125*hi[4]*K[3]/192.0) + (-2187*hi[4]*K[4]/6784.0) + 11*hi[4]*K[5]/84.0)+Kvect.at(30)*(((DdZ[9]+ hi[4])/Kvect.at(4))-Kvect.at(1))*Heavi2);

        dXt   = DdZ[7] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);                                                // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[7] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dXt);
        s = pow(eps*hi[4]/(2*err),1/5);
        hopt[4] = s*hi[4];
        if( hopt[4] < hmin) hopt[4] = hmin;
        else if(hopt[4] > hmax) hopt[4] = hmax;

        // Решение dIl через Дормана-Принса

        K[0] = hi[5]*(-(Kvect.at(5)+Kvect.at(7))*DdZ[11]+Kvect.at(6)*DdZ[9]);
        K[1] = hi[5]*(-(Kvect.at(5)+Kvect.at(7))*(DdZ[11]+ hi[5]*K[0]/5.0)+Kvect.at(6)*(DdZ[9]+ hi[5]/5.0));
        K[2] = hi[5]*(-(Kvect.at(5)+Kvect.at(7))*(DdZ[11]+ 3*hi[5]*K[0]/40.0 + 9*hi[5]*K[1]/40.0)+Kvect.at(6)*(DdZ[9]+ 3*hi[5]/10.0));
        K[3] = hi[5]*(-(Kvect.at(5)+Kvect.at(7))*(DdZ[11]+ 44*hi[5]*K[0]/45.0 + (-56*hi[5]*K[1]/15.0) + 32*hi[5]*K[2]/9.0)+Kvect.at(6)*(DdZ[9]+ 4*hi[5]/5.0));
        K[4] = hi[5]*(-(Kvect.at(5)+Kvect.at(7))*(DdZ[11]+ 19372*hi[5]*K[0]/6561.0 + (-25360*hi[5]*K[1]/2187.0) + 64448*hi[5]*K[2]/6561.0 + (-212*hi[5]*K[3]/729.0))+Kvect.at(6)*(DdZ[9]+ 8*hi[5]/9.0));
        K[5] = hi[5]*(-(Kvect.at(5)+Kvect.at(7))*(DdZ[11]+ 9017*hi[5]*K[0]/3168.0 + (-355*hi[5]*K[1]/33.0) + 46732*hi[5]*K[2]/5247.0 + (49*hi[5]*K[3]/176.0) + (-5103*hi[5]*K[4]/18656.0))+Kvect.at(6)*(DdZ[9]+ hi[5]));
        K[6] = hi[5]*(-(Kvect.at(5)+Kvect.at(7))*(DdZ[11]+ 35*hi[5]*K[0]/384.0 + 500*hi[5]*K[2]/1113.0 + (125*hi[5]*K[3]/192.0) + (-2187*hi[5]*K[4]/6784.0) + 11*hi[5]*K[5]/84.0)+Kvect.at(6)*(DdZ[9]+ hi[5]));

        dIl   = DdZ[11] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);                                                // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[11] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIl);
        s = pow(eps*hi[5]/(2*err),1/5);
        hopt[5] = s*hi[5];
        if( hopt[5] < hmin) hopt[5] = hmin;
        else if(hopt[5] > hmax) hopt[5] = hmax;

        // Решение dIp через Дормана-Принса

        K[0] = hi[6]*(-Kvect.at(6)*DdZ[9]+Kvect.at(5)*DdZ[8]+Kvect.at(10)/Kvect.at(0)*DdZ[10]-Kvect.at(9)*DdZ[9]);
        K[1] = hi[6]*(-Kvect.at(6)*(DdZ[9]+ hi[6]*K[0]/5.0)+Kvect.at(5)*(DdZ[8]+ hi[6]/5.0)+Kvect.at(10)/Kvect.at(0)*(DdZ[10]+ hi[6]/5.0)-Kvect.at(9)*(DdZ[9]+ hi[6]*K[0]/5.0));
        K[2] = hi[6]*(-Kvect.at(6)*(DdZ[9]+ 3*hi[6]*K[0]/40.0 + 9*hi[6]*K[1]/40.0)+Kvect.at(5)*(DdZ[8]+ 3*hi[6]/10.0)+Kvect.at(10)/Kvect.at(0)*(DdZ[10]+ 3*hi[6]/10.0)-Kvect.at(9)*(DdZ[9]+ 3*hi[6]*K[0]/40.0 + 9*hi[6]*K[1]/40.0));
        K[3] = hi[6]*(-Kvect.at(6)*(DdZ[9]+ 44*hi[6]*K[0]/45.0 + (-56*hi[6]*K[1]/15.0) + 32*hi[6]*K[2]/9.0)+Kvect.at(5)*(DdZ[8]+ 4*hi[6]/5.0)+Kvect.at(10)/Kvect.at(0)*(DdZ[10]+ 4*hi[6]/5.0)-Kvect.at(9)*(DdZ[9]+ 44*hi[6]*K[0]/45.0 + (-56*hi[6]*K[1]/15.0) + 32*hi[6]*K[2]/9.0));
        K[4] = hi[6]*(-Kvect.at(6)*(DdZ[9]+ 19372*hi[6]*K[0]/6561.0 + (-25360*hi[6]*K[1]/2187.0) + 64448*hi[6]*K[2]/6561.0 + (-212*hi[6]*K[3]/729.0) )+Kvect.at(5)*(DdZ[8]+ 8*hi[6]/9.0)+Kvect.at(10)/Kvect.at(0)*(DdZ[10]+ 8*hi[6]/9.0)-Kvect.at(9)*(DdZ[9]+ 19372*hi[6]*K[0]/6561.0 + (-25360*hi[6]*K[1]/2187.0) + 64448*hi[6]*K[2]/6561.0 + (-212*hi[6]*K[3]/729.0) ));
        K[5] = hi[6]*(-Kvect.at(6)*(DdZ[9]+ 9017*hi[6]*K[0]/3168.0 + (-355*hi[6]*K[1]/33.0) + 46732*hi[6]*K[2]/5247.0 + (49*hi[6]*K[3]/176.0) + (-5103*hi[6]*K[4]/18656.0))+Kvect.at(5)*(DdZ[8]+ hi[6])+Kvect.at(10)/Kvect.at(0)*(DdZ[10]+ hi[6])-Kvect.at(9)*(DdZ[9]+ 9017*hi[6]*K[0]/3168.0 + (-355*hi[6]*K[1]/33.0) + 46732*hi[6]*K[2]/5247.0 + (49*hi[6]*K[3]/176.0) + (-5103*hi[6]*K[4]/18656.0)));
        K[6] = hi[6]*(-Kvect.at(6)*(DdZ[9]+ 35*hi[6]*K[0]/384.0 + 500*hi[6]*K[2]/1113.0 + (125*hi[6]*K[3]/192.0) + (-2187*hi[6]*K[4]/6784.0) + 11*hi[6]*K[5]/84.0)+Kvect.at(5)*(DdZ[8]+ hi[6])+Kvect.at(10)/Kvect.at(0)*(DdZ[10]+ hi[6])-Kvect.at(9)*(DdZ[9]+ 35*hi[6]*K[0]/384.0 + 500*hi[6]*K[2]/1113.0 + (125*hi[6]*K[3]/192.0) + (-2187*hi[6]*K[4]/6784.0) + 11*hi[6]*K[5]/84.0));

        dIp   = DdZ[9] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[9] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIp);
        s = pow(eps*hi[6]/(2*err),1/5);
        hopt[6] = s*hi[6];
        if( hopt[6] < hmin) hopt[6] = hmin;
        else if(hopt[6] > hmax) hopt[6] = hmax;


        // Решение dfgut через Дормана-Принса
        double kgut;
        kgut=Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(DdZ[4]+DdZ[5]-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(DdZ[4]+DdZ[5]-Kvect.at(22)*Dig))+2);

        K[0] = hi[7]*(-Kvect.at(17)*DdZ[3]+kgut*DdZ[5]);
        K[1] = hi[7]*(-Kvect.at(17)*(DdZ[3]+ hi[7]*K[0]/5.0)+kgut*(DdZ[5]+ hi[7]/5.0));
        K[2] = hi[7]*(-Kvect.at(17)*(DdZ[3]+ 3*hi[7]*K[0]/40.0 + 9*hi[7]*K[1]/40.0)+kgut*(DdZ[5]+ 3*hi[7]/10.0));
        K[3] = hi[7]*(-Kvect.at(17)*(DdZ[3]+ 44*hi[7]*K[0]/45.0 + (-56*hi[7]*K[1]/15.0) + 32*hi[7]*K[2]/9.0)+kgut*(DdZ[5]+ 4*hi[7]/5.0));
        K[4] = hi[7]*(-Kvect.at(17)*(DdZ[3]+ 19372*hi[7]*K[0]/6561.0 + (-25360*hi[7]*K[1]/2187.0) + 64448*hi[7]*K[2]/6561.0 + (-212*hi[7]*K[3]/729.0) )+kgut*(DdZ[5]+ 8*hi[7]/9.0));
        K[5] = hi[7]*(-Kvect.at(17)*(DdZ[3]+ 9017*hi[7]*K[0]/3168.0 + (-355*hi[7]*K[1]/33.0) + 46732*hi[7]*K[2]/5247.0 + (49*hi[7]*K[3]/176.0) + (-5103*hi[7]*K[4]/18656.0))+kgut*(DdZ[5]+ hi[7]));
        K[6] = hi[7]*(-Kvect.at(17)*(DdZ[3]+ 35*hi[7]*K[0]/384.0 + 500*hi[7]*K[2]/1113.0 + (125*hi[7]*K[3]/192.0) + (-2187*hi[7]*K[4]/6784.0) + 11*hi[7]*K[5]/84.0)+kgut*(DdZ[5]+ hi[7]));

        dfgut   = DdZ[3] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[3] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dfgut);
        s = pow(eps*hi[7]/(2*err),1/5);
        hopt[7] = s*hi[7];
        if( hopt[7] < hmin) hopt[7] = hmin;
        else if(hopt[7] > hmax) hopt[7] = hmax;

        // Решение dfliq через Дормана-Принса kgut

        K[0] = hi[8]*(-kgut*DdZ[5]+Kvect.at(18)*DdZ[4]);
        K[1] = hi[8]*(-kgut*(DdZ[5]+ hi[8]*K[0]/5.0)+Kvect.at(18)*(DdZ[4]+ hi[8]/5.0));
        K[2] = hi[8]*(-kgut*(DdZ[5]+ 3*hi[8]*K[0]/40.0 + 9*hi[8]*K[1]/40.0)+Kvect.at(18)*(DdZ[4]+ 3*hi[8]/10.0));
        K[3] = hi[8]*(-kgut*(DdZ[5]+ 44*hi[8]*K[0]/45.0 + (-56*hi[8]*K[1]/15.0) + 32*hi[8]*K[2]/9.0)+Kvect.at(18)*(DdZ[4]+ 4*hi[8]/5.0));
        K[4] = hi[8]*(-kgut*(DdZ[5]+ 19372*hi[8]*K[0]/6561.0 + (-25360*hi[8]*K[1]/2187.0) + 64448*hi[8]*K[2]/6561.0 + (-212*hi[8]*K[3]/729.0))+Kvect.at(18)*(DdZ[4]+ 8*hi[8]/9.0));
        K[5] = hi[8]*(-kgut*(DdZ[5]+ 9017*hi[8]*K[0]/3168.0 + (-355*hi[8]*K[1]/33.0) + 46732*hi[8]*K[2]/5247.0 + (49*hi[8]*K[3]/176.0) + (-5103*hi[8]*K[4]/18656.0))+Kvect.at(18)*(DdZ[4]+ hi[8]));
        K[6] = hi[8]*(-kgut*(DdZ[5]+ 35*hi[8]*K[0]/384.0 + 500*hi[8]*K[2]/1113.0 + (125*hi[8]*K[3]/192.0) + (-2187*hi[8]*K[4]/6784.0) + 11*hi[8]*K[5]/84.0)+Kvect.at(18)*(DdZ[4]+ hi[8]));

        dfliq   = DdZ[5] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[5] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dfliq);
        s = pow(eps*hi[8]/(2*err),1/5);
        hopt[8] = s*hi[8];
        if( hopt[8] < hmin) hopt[8] = hmin;
        else if(hopt[8] > hmax) hopt[8] = hmax;

        // Решение dfsol через Дормана-Принса

        K[0] = hi[9]*(-Kvect.at(18)*DdZ[4]+vm);
        K[1] = hi[9]*(-Kvect.at(18)*(DdZ[4]+ hi[9]*K[0]/5.0)+(vm));
        K[2] = hi[9]*(-Kvect.at(18)*(DdZ[4]+ 3*hi[9]*K[0]/40.0 + 9*hi[9]*K[1]/40.0)+(vm));
        K[3] = hi[9]*(-Kvect.at(18)*(DdZ[4]+ 44*hi[9]*K[0]/45.0 + (-56*hi[9]*K[1]/15.0) + 32*hi[9]*K[2]/9.0)+(vm));
        K[4] = hi[9]*(-Kvect.at(18)*(DdZ[4]+ 19372*hi[9]*K[0]/6561.0 + (-25360*hi[9]*K[1]/2187.0) + 64448*hi[9]*K[2]/6561.0 + (-212*hi[9]*K[3]/729.0))+(vm));
        K[5] = hi[9]*(-Kvect.at(18)*(DdZ[4]+ 9017*hi[9]*K[0]/3168.0 + (-355*hi[9]*K[1]/33.0) + 46732*hi[9]*K[2]/5247.0 + (49*hi[9]*K[3]/176.0) + (-5103*hi[9]*K[4]/18656.0))+(vm));
        K[6] = hi[9]*(-Kvect.at(18)*(DdZ[4]+ 35*hi[9]*K[0]/384.0 + 500*hi[9]*K[2]/1113.0 + (125*hi[9]*K[3]/192.0) + (-2187*hi[9]*K[4]/6784.0) + 11*hi[9]*K[5]/84.0)+(vm));

        dfsol   = DdZ[4] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[4] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dfsol);
        s = pow(eps*hi[9]/(2*err),1/5);
        hopt[9] = s*hi[9];
        if( hopt[9] < hmin) hopt[9] = hmin;
        else if(hopt[9] > hmax) hopt[9] = hmax;

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
        K[0] = hi[10]*(-Kvect.at(33)*DdZ[2]+(DdZ[12]+Kvect.at(3))*Heavi5+(DdZ[12]+Kvect.at(3))*(Heavi6));
        K[1] = hi[10]*(-Kvect.at(33)*(DdZ[2]+ hi[10]*K[0]/5.0)+((DdZ[12]+ hi[10]/5.0)+Kvect.at(3))*Heavi5+((DdZ[12]+ hi[10]/5.0)+Kvect.at(3))*(Heavi6));
        K[2] = hi[10]*(-Kvect.at(33)*(DdZ[2]+ 3*hi[10]*K[0]/40.0 + 9*hi[10]*K[1]/40.0)+((DdZ[12]+ 3*hi[10]/10.0)+Kvect.at(3))*Heavi5+((DdZ[12]+ 3*hi[10]/10.0)+Kvect.at(3))*(Heavi6));
        K[3] = hi[10]*(-Kvect.at(33)*(DdZ[2]+ 44*hi[10]*K[0]/45.0 + (-56*hi[10]*K[1]/15.0) + 32*hi[10]*K[2]/9.0)+((DdZ[12]+ 4*hi[10]/5.0)+Kvect.at(3))*Heavi5+((DdZ[12]+ 4*hi[10]/5.0)+Kvect.at(3))*(Heavi6));
        K[4] = hi[10]*(-Kvect.at(33)*(DdZ[2]+ 19372*hi[10]*K[0]/6561.0 + (-25360*hi[10]*K[1]/2187.0) + 64448*hi[10]*K[2]/6561.0 + (-212*hi[10]*K[3]/729.0))+((DdZ[12]+ 8*hi[10]/9.0)+Kvect.at(3))*Heavi5+((DdZ[12]+ 8*hi[10]/9.0)+Kvect.at(3))*(Heavi6));
        K[5] = hi[10]*(-Kvect.at(33)*(DdZ[2]+ 9017*hi[10]*K[0]/3168.0 + (-355*hi[10]*K[1]/33.0) + 46732*hi[10]*K[2]/5247.0 + (49*hi[10]*K[3]/176.0) + (-5103*hi[10]*K[4]/18656.0))+((DdZ[12]+ hi[10])+Kvect.at(3))*Heavi5+((DdZ[12]+ hi[10])+Kvect.at(3))*(Heavi6));
        K[6] = hi[10]*(-Kvect.at(33)*(DdZ[2]+ 35*hi[10]*K[0]/384.0 + 500*hi[10]*K[2]/1113.0 + (125*hi[10]*K[3]/192.0) + (-2187*hi[10]*K[4]/6784.0) + 11*hi[10]*K[5]/84.0)+((DdZ[12]+ hi[10])+Kvect.at(3))*Heavi5+((DdZ[12]+ hi[10])+Kvect.at(3))*(Heavi6));

        dIpo   = DdZ[2] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[2] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIpo);
        s = pow(eps*hi[10]/(2*err),1/5);
        hopt[10] = s*hi[10];
        if( hopt[10] < hmin) hopt[10] = hmin;
        else if(hopt[10] > hmax) hopt[10] = hmax;

        // Решение dYt через Дормана-Принса

        double Heavi7;
        double Heavi8;

        if( (Kvect.at(35)*(DdZ[0]/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*(DdZ[0]/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*(DdZ[0]/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*(DdZ[0]/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[0] = hi[11]*(-Kvect.at(34)*(DdZ[12]-Kvect.at(35)*(DdZ[0]/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*DdZ[12]-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((DdZ[0]+ hi[11]/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((DdZ[0]+ hi[11]/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ hi[11]/5.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ hi[11]/5.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[1] = hi[11]*(-Kvect.at(34)*((DdZ[12]+ hi[11]*K[0]/5.0)-Kvect.at(35)*((DdZ[0]+ hi[11]/5.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(DdZ[12]+ hi[11]*K[0]/5.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((DdZ[0]+ 3*hi[11]/10.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((DdZ[0]+ 3*hi[11]/10.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ 3*hi[11]/10.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ 3*hi[11]/10.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[2] = hi[11]*(-Kvect.at(34)*((DdZ[12]+ 3*hi[11]*K[0]/40.0 + 9*hi[11]*K[1]/40.0)-Kvect.at(35)*((DdZ[0]+ 3*hi[11]/10.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(DdZ[12]+ 3*hi[11]*K[0]/40.0 + 9*hi[11]*K[1]/40.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((DdZ[0]+ 4*hi[11]/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((DdZ[0]+ 4*hi[11]/5.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ 4*hi[11]/5.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ 4*hi[11]/5.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[3] = hi[11]*(-Kvect.at(34)*((DdZ[12]+ 44*hi[11]*K[0]/45.0 + (-56*hi[11]*K[1]/15.0) + 32*hi[11]*K[2]/9.0)-Kvect.at(35)*((DdZ[0]+ 4*hi[11]/5.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(DdZ[12]+ 44*hi[11]*K[0]/45.0 + (-56*hi[11]*K[1]/15.0) + 32*hi[11]*K[2]/9.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((DdZ[0]+ 8*hi[11]/9.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((DdZ[0]+ 8*hi[11]/9.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ 8*hi[11]/9.0)/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ 8*hi[11]/9.0)/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[4] = hi[11]*(-Kvect.at(34)*((DdZ[12]+ 19372*hi[11]*K[0]/6561.0 + (-25360*hi[11]*K[1]/2187.0) + 64448*hi[11]*K[2]/6561.0 + (-212*hi[11]*K[3]/729.0))-Kvect.at(35)*((DdZ[0]+ 8*hi[11]/9.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(DdZ[12]+ 19372*hi[11]*K[0]/6561.0 + (-25360*hi[11]*K[1]/2187.0) + 64448*hi[11]*K[2]/6561.0 + (-212*hi[11]*K[3]/729.0))-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[5] = hi[11]*(-Kvect.at(34)*((DdZ[12]+ 9017*hi[11]*K[0]/3168.0 + (-355*hi[11]*K[1]/33.0) + 46732*hi[11]*K[2]/5247.0 + (49*hi[11]*K[3]/176.0) + (-5103*hi[11]*K[4]/18656.0))-Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(DdZ[12]+ 9017*hi[11]*K[0]/3168.0 + (-355*hi[11]*K[1]/33.0) + 46732*hi[11]*K[2]/5247.0 + (49*hi[11]*K[3]/176.0) + (-5103*hi[11]*K[4]/18656.0))-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        if( (Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) >= 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}

        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2))) >= 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        K[6] = hi[11]*(-Kvect.at(34)*((DdZ[12]+ 35*hi[11]*K[0]/384.0 + 500*hi[11]*K[2]/1113.0 + (125*hi[11]*K[3]/192.0) + (-2187*hi[11]*K[4]/6784.0) + 11*hi[11]*K[5]/84.0)-Kvect.at(35)*((DdZ[0]+ hi[11])/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(DdZ[12]+ 35*hi[11]*K[0]/384.0 + 500*hi[11]*K[2]/1113.0 + (125*hi[11]*K[3]/192.0) + (-2187*hi[11]*K[4]/6784.0) + 11*hi[11]*K[5]/84.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        dYt   = DdZ[12] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[12] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dYt);
        s = pow(eps*hi[11]/(2*err),1/5);
        hopt[11] = s*hi[11];
        if( hopt[11] < hmin) hopt[11] = hmin;
        else if(hopt[11] > hmax) hopt[11] = hmax;

        // Решение dIt через Дормана-Принса

        K[0] = hi[12]*(Kvect.at(8)*DdZ[13]-Kvect.at(10)*DdZ[10]);
        K[1] = hi[12]*(Kvect.at(8)*(DdZ[13]+ hi[12]/5.0)-Kvect.at(10)*(DdZ[10]+ hi[12]*K[0]/5.0));
        K[2] = hi[12]*(Kvect.at(8)*(DdZ[13]+ 3*hi[12]/10.0)-Kvect.at(10)*(DdZ[10]+ 3*hi[12]*K[0]/40.0 + 9*hi[12]*K[1]/40.0));
        K[3] = hi[12]*(Kvect.at(8)*(DdZ[13]+ 4*hi[12]/5.0)-Kvect.at(10)*(DdZ[10]+ 44*hi[12]*K[0]/45.0 + (-56*hi[12]*K[1]/15.0) + 32*hi[12]*K[2]/9.0));
        K[4] = hi[12]*(Kvect.at(8)*(DdZ[13]+ 8*hi[12]/9.0)-Kvect.at(10)*(DdZ[10]+ 19372*hi[12]*K[0]/6561.0 + (-25360*hi[12]*K[1]/2187.0) + 64448*hi[12]*K[2]/6561.0 + (-212*hi[12]*K[3]/729.0)));
        K[5] = hi[12]*(Kvect.at(8)*(DdZ[13]+ hi[12])-Kvect.at(10)*(DdZ[10]+ 9017*hi[12]*K[0]/3168.0 + (-355*hi[12]*K[1]/33.0) + 46732*hi[12]*K[2]/5247.0 + (49*hi[12]*K[3]/176.0) + (-5103*hi[12]*K[4]/18656.0)));
        K[6] = hi[12]*(Kvect.at(8)*(DdZ[13]+ hi[12])-Kvect.at(10)*(DdZ[10]+ 35*hi[12]*K[0]/384.0 + 500*hi[12]*K[2]/1113.0 + (125*hi[12]*K[3]/192.0) + (-2187*hi[12]*K[4]/6784.0) + 11*hi[12]*K[5]/84.0));

        dIt   = DdZ[10] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[10] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIt);
        s = pow(eps*hi[12]/(2*err),1/5);
        hopt[12] = s*hi[12];
        if( hopt[12] < hmin) hopt[12] = hmin;
        else if(hopt[12] > hmax) hopt[12] = hmax;

        // Решение dIi через Дормана-Принса

        K[0] = hi[13]*(-Kvect.at(8)*DdZ[13]+vbas+vbol);
        K[1] = hi[13]*(-Kvect.at(8)*(DdZ[13]+ hi[13]*K[0]/5.0)+vbas+vbol);
        K[2] = hi[13]*(-Kvect.at(8)*(DdZ[13]+ 3*hi[13]*K[0]/40.0 + 9*hi[13]*K[1]/40.0)+(vbas)+vbol);
        K[3] = hi[13]*(-Kvect.at(8)*(DdZ[13]+ 44*hi[13]*K[0]/45.0 + (-56*hi[13]*K[1]/15.0) + 32*hi[13]*K[2]/9.0)+(vbas)+vbol);
        K[4] = hi[13]*(-Kvect.at(8)*(DdZ[13]+ 19372*hi[13]*K[0]/6561.0 + (-25360*hi[13]*K[1]/2187.0) + 64448*hi[13]*K[2]/6561.0 + (-212*hi[13]*K[3]/729.0))+(vbas)+(vbol));
        K[5] = hi[13]*(-Kvect.at(8)*(DdZ[13]+ 9017*hi[13]*K[0]/3168.0 + (-355*hi[13]*K[1]/33.0) + 46732*hi[13]*K[2]/5247.0 + (49*hi[13]*K[3]/176.0) + (-5103*hi[13]*K[4]/18656.0))+(vbas)+(vbol));
        K[6] = hi[13]*(-Kvect.at(8)*(DdZ[13]+ 35*hi[13]*K[0]/384.0 + 500*hi[13]*K[2]/1113.0 + (125*hi[13]*K[3]/192.0) + (-2187*hi[13]*K[4]/6784.0) + 11*hi[13]*K[5]/84.0)+(vbas)+(vbol));

        dIi   = DdZ[13] + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0);            // во многих источниках домножается на h, но в одном нету такого.
        dz = DdZ[13] + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0);
        err = abs(dz-dIi);
        s = pow(eps*hi[13]/(2*err),1/5);
        hopt[13] = s*hi[13];
        if( hopt[13] < hmin) hopt[13] = hmin;
        else if(hopt[13] > hmax) hopt[13] = hmax;

        DdZ[0] = dgp;
        DdZ[1] = dId;
        DdZ[2] = dIpo;
        DdZ[3] = dfgut;
        DdZ[4] = dfsol;
        DdZ[5] = dfliq;
        DdZ[6] = dgt;
        DdZ[7] = dXt;
        DdZ[8] = dI1;
        DdZ[9] = dIp;
        DdZ[10] = dIt;
        DdZ[11] = dIl;
        DdZ[12] = dYt;
        DdZ[13] = dIi;

        dptick.append(tDP);
        dpCGB.append(DdZ[0]/Kvect.at(23));
        dpIpg.append(DdZ[9]*20);

        hi[0] = hopt[0];
        hi[1] = hopt[1];
        hi[2] = hopt[2];
        hi[3] = hopt[3];
        hi[4] = hopt[4];
        hi[5] = hopt[5];
        hi[6] = hopt[6];
        hi[7] = hopt[7];
        hi[8] = hopt[8];
        hi[9] = hopt[9];
        hi[10] = hopt[10];
        hi[11] = hopt[11];
        hi[12] = hopt[12];
        hi[13] = hopt[13];

        tDP = tDP + h; // увеличиваем время
        j = j + h;
    }

}
void MainWindow::RungeKutta(int *i)
{
    double  a = *i*5 - 5;    // промедуток 5 минут
    double  b = *i*5;        // конечное время

    double h = 0.04;//0.25;

    double OB = 6.63;//(m2/19) *1.5;      // result of bolus calculation Потом переместить в Глобальный цикл

    double KRK[4];

    double del = 0.4;       // result of bolus calculation +
    double m2 = 90;         // + М ?

    double Dig = 1176*m2+1; // вместо m2 массу углеводов
    double Dbol1=del*OB;
    double Dbol2=(1-del)*OB;
    double Dbol3=0;

    double Heavi3;
    if( (Kvect.at(11)-Kvect.at(12)*dZ[0]-Kvect.at(13)*dZ[1]-Kvect.at(14)*dZ[2])  > 0){
         Heavi3 = 1;}
    if( (Kvect.at(11)-Kvect.at(12)*dZ[0]-Kvect.at(13)*dZ[1]-Kvect.at(14)*dZ[2]) < 0){
         Heavi3 = 0;}
    if( (Kvect.at(11)-Kvect.at(12)*dZ[0]-Kvect.at(13)*dZ[1]-Kvect.at(14)*dZ[2]) == 0){
         Heavi3 = 0.5;}

    double EGP;
    EGP=(Kvect.at(11)-Kvect.at(12)*dZ[0]-Kvect.at(13)*dZ[1]-Kvect.at(14)*dZ[2])*Heavi3;
    double kgut;
    kgut=Kvect.at(19)+(Kvect.at(20)-Kvect.at(19))/2*(tanh((5/(2*Dig*(1-Kvect.at(21))))*(dZ[4]+dZ[5]-Kvect.at(21)*Dig))-tanh((5/(2*Dig*Kvect.at(22)))*(dZ[4]+dZ[5]-Kvect.at(22)*Dig))+2);

    /* bolus */
    double ti1 = 30;
    double ti2 = 60;
    double ti3 = 10;

    double Ti1 = 10;
    double Ti2 = 10;
    double Ti3 = 10;


    double tm1 = 60; // tm  // время приёма пищи
    double Tm = 20;         // длительность приёма пищи
    // Fluctuations
    double vbas=Vbas*100;      // pmol/min
    double bol1=1/Ti1*Dbol1*(1./(1+exp(-3*(t+10-ti1))))*(1./(1+exp(-3*(-10+ti1-t+Ti1))));
    double bol2=1/Ti2*Dbol2*(1./(1+exp(-3*(t+10-ti2))))*(1./(1+exp(-3*(-10+ti2-t+Ti2))));
    double bol3=1/Ti3*Dbol3*(1./(1+exp(-3*(t+10-ti3))))*(1./(1+exp(-3*(-10+ti3-t+Ti3))));

    double vbol=6000*(bol1+bol2+bol3);

    double vm=Dig/Tm*(1./(1+exp(-3*(t-tm1))))*(1./(1+exp(-3*(-(t-tm1-Tm)))));

    for (double j = a; j <= b;)

    {

        // выходные значения dZout
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

        /* ДУ */

    // Решаем dgp через Рунге-Кутта
    //dgp=EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*fgut-Kvect.at(26)-Kvect.at(31)*(gp-Kvect.at(32))*Heavi1-Kvect.at(24)*gp+Kvect.at(25)*gt; // plasma glucose mg/dl
    // dgp = f(fgut,gt,Id,Ipo,fgut,gp)

        double Heavi1;
        if( (dZ[0]-Kvect.at(32)) > 0){
             Heavi1 = 1;}
        if( (dZ[0]-Kvect.at(32)) < 0){
             Heavi1 = 0;}
        if( (dZ[0]-Kvect.at(32)) == 0){
             Heavi1 = 0.5;}

        KRK[0] = h*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*dZ[3]-Kvect.at(26)-Kvect.at(31)*(dZ[0]-Kvect.at(32))*Heavi1-Kvect.at(24)*dZ[0]+Kvect.at(25)*dZ[6]);

            if( ((dZ[0]+ h*KRK[0]/2.0)-Kvect.at(32)) > 0){
                 Heavi1 = 1;}
            if( ((dZ[0]+ h*KRK[0]/2.0)-Kvect.at(32)) < 0){
                 Heavi1 = 0;}
            if( ((dZ[0]+ h*KRK[0]/2.0)-Kvect.at(32)) == 0){
                 Heavi1 = 0.5;}

        KRK[1] = h*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(dZ[3]+ h/2.0)-Kvect.at(26)-Kvect.at(31)*((dZ[0]+ h*KRK[0]/2.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(dZ[0]+ h*KRK[0]/2.0)+Kvect.at(25)*(dZ[6]+ h/2.0));

            if( ((dZ[0]+ h*KRK[1]/2.0)-Kvect.at(32)) > 0){
                 Heavi1 = 1;}
            if( ((dZ[0]+ h*KRK[1]/2.0)-Kvect.at(32)) < 0){
                 Heavi1 = 0;}
            if( ((dZ[0]+ h*KRK[1]/2.0)-Kvect.at(32)) == 0){
                 Heavi1 = 0.5;}

        KRK[2] = h*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(dZ[3]+ h/2.0)-Kvect.at(26)-Kvect.at(31)*((dZ[0]+ h*KRK[1]/2.0)-Kvect.at(32))*Heavi1-Kvect.at(24)*(dZ[0]+ h*KRK[1]/2.0)+Kvect.at(25)*(dZ[6]+ h/2.0));

            if( ((dZ[0]+ h*KRK[2])-Kvect.at(32)) > 0){
                 Heavi1 = 1;}
            if( ((dZ[0]+ h*KRK[2])-Kvect.at(32)) < 0){
                 Heavi1 = 0;}
            if( ((dZ[0]+ h*KRK[2])-Kvect.at(32)) == 0){
                 Heavi1 = 0.5;}

        KRK[3] = h*(EGP+Kvect.at(16)/Kvect.at(0)*Kvect.at(17)*(dZ[3]+ h)-Kvect.at(26)-Kvect.at(31)*((dZ[0]+ h*KRK[2])-Kvect.at(32))*Heavi1-Kvect.at(24)*(dZ[0]+ h*KRK[2])+Kvect.at(25)*(dZ[6]+ h));

        dgp   = dZ[0] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dIl через Рунге-Кутта
        // dIl=-(Kvect.at(5)+Kvect.at(7))*Il+Kvect.at(6)*Ip;
        // dIl = f(Ip,Il)?

        KRK[0] = h*(-(Kvect.at(5)+Kvect.at(7))*dZ[11]+Kvect.at(6)*dZ[9]);
        KRK[1] = h*(-(Kvect.at(5)+Kvect.at(7))*(dZ[11]+ h*KRK[0]/2.0)+Kvect.at(6)*(dZ[9]+ h/2.0));
        KRK[2] = h*(-(Kvect.at(5)+Kvect.at(7))*(dZ[11]+ h*KRK[1]/2.0)+Kvect.at(6)*(dZ[9]+ h/2.0));
        KRK[3] = h*(-(Kvect.at(5)+Kvect.at(7))*(dZ[11]+ h*KRK[2])+Kvect.at(6)*(dZ[9]+ h));

        dIl   = dZ[11] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dIp через Рунге-Кутта
        // dIp=-Kvect.at(6)*Ip+Kvect.at(5)*Il+Kvect.at(10)/Kvect.at(0)*It-Kvect.at(9)*Ip;
        // dIp = f(Il,It,Ip)?
        KRK[0] = h*(-Kvect.at(6)*dZ[9]+Kvect.at(5)*dZ[11]+Kvect.at(10)/Kvect.at(0)*dZ[10]-Kvect.at(9)*dZ[9]);
        KRK[1] = h*(-Kvect.at(6)*(dZ[9]+ h*KRK[0]/2.0)+Kvect.at(5)*(dZ[11]+ h/2.0)+Kvect.at(10)/Kvect.at(0)*(dZ[10]+ h/2.0)-Kvect.at(9)*(dZ[9]+ h*KRK[0]/2.0));
        KRK[2] = h*(-Kvect.at(6)*(dZ[9]+ h*KRK[1]/2.0)+Kvect.at(5)*(dZ[11]+ h/2.0)+Kvect.at(10)/Kvect.at(0)*(dZ[10]+ h/2.0)-Kvect.at(9)*(dZ[9]+ h*KRK[1]/2.0));
        KRK[3] = h*(-Kvect.at(6)*(dZ[9]+ h*KRK[2])+Kvect.at(5)*(dZ[11]+ h)+Kvect.at(10)/Kvect.at(0)*(dZ[10]+ h)-Kvect.at(9)*(dZ[9]+ h*KRK[2]));

        dIp   = dZ[9] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dfgut через Рунге-Кутта
        // dfgut=-Kvect.at(17)*fgut+kgut*fliq;
        // dfgut = f(fliq,fsol,fgut)?

        KRK[0] = h*(-Kvect.at(17)*dZ[3]+kgut*dZ[5]);
        KRK[1] = h*(-Kvect.at(17)*(dZ[3]+ h*KRK[0]/2.0)+kgut*(dZ[5]+ h/2.0));
        KRK[2] = h*(-Kvect.at(17)*(dZ[3]+ h*KRK[1]/2.0)+kgut*(dZ[5]+ h/2.0));
        KRK[3] = h*(-Kvect.at(17)*(dZ[3]+ h*KRK[2])+kgut*(dZ[5]+ h));

        dfgut   = dZ[3] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dfliq через Рунге-Кутта
        // dfliq=-kgut*fliq+Kvect.at(18)*fsol;
        // dfliq = f(fliq,fsol)?

        KRK[0] = h*(-kgut*dZ[5]+Kvect.at(18)*dZ[4]);
        KRK[1] = h*(-kgut*(dZ[5]+ h*KRK[0]/2.0)+Kvect.at(18)*(dZ[4]+ h/2.0));
        KRK[2] = h*(-kgut*(dZ[5]+ h*KRK[1]/2.0)+Kvect.at(18)*(dZ[4]+ h/2.0));
        KRK[3] = h*(-kgut*(dZ[5]+ h*KRK[2])+Kvect.at(18)*(dZ[4]+ h));

        dfliq   = dZ[5] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dfsol через Рунге-Кутта
        // dfsol=-Kvect.at(18)*fsol+vm;
        // dfsol = f(fsol,vm)?

        KRK[0] = h*(-Kvect.at(18)*dZ[4]+vm);
        KRK[1] = h*(-Kvect.at(18)*(dZ[4]+ h*KRK[0]/2.0)+(vm));
        KRK[2] = h*(-Kvect.at(18)*(dZ[4]+ h*KRK[1]/2.0)+(vm));
        KRK[3] = h*(-Kvect.at(18)*(dZ[4]+ h*KRK[2])+(vm));

        dfsol   = dZ[4] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dgt через Рунге-Кутта
        // dgt=-((Kvect.at(27)+Kvect.at(28)*Xt)*gt)/(Kvect.at(29)+gt)+Kvect.at(24)*gp-Kvect.at(25)*gt;
        // dgt = f(Xt,gp,gt)

        KRK[0] = h*(-((Kvect.at(27)+Kvect.at(28)*dZ[7])*dZ[6])/(Kvect.at(29)+dZ[6])+Kvect.at(24)*dZ[0]-Kvect.at(25)*dZ[6]);
        KRK[1] = h*(-((Kvect.at(27)+Kvect.at(28)*(dZ[7] + h/2.0))*(dZ[6] + h*KRK[0]/2.0))/(Kvect.at(29)+(dZ[6] + h*KRK[0]/2.0))+Kvect.at(24)*(dZ[0] + h/2.0)-Kvect.at(25)*(dZ[6] + h*KRK[0]/2.0));
        KRK[2] = h*(-((Kvect.at(27)+Kvect.at(28)*(dZ[7] + h/2.0))*(dZ[6] + h*KRK[1]/2.0))/(Kvect.at(29)+(dZ[6] + h*KRK[1]/2.0))+Kvect.at(24)*(dZ[0] + h/2.0)-Kvect.at(25)*(dZ[6] + h*KRK[1]/2.0));
        KRK[3] = h*(-((Kvect.at(27)+Kvect.at(28)*(dZ[7] + h))*(dZ[6] + h*KRK[2]))/(Kvect.at(29)+(dZ[6] + h*KRK[2]))+Kvect.at(24)*(dZ[0] + h)-Kvect.at(25)*(dZ[6] + h*KRK[2]));
        dgt   = dZ[6] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dI1 через Рунге-Кутта
        // dI1=-Kvect.at(15)*(I1-Ip/Kvect.at(4));
        // dI1 = f(Ip,I1)

        KRK[0] = h*(-Kvect.at(15)*(dZ[8]-dZ[9]/Kvect.at(4)));
        KRK[1] = h*(-Kvect.at(15)*((dZ[8] + h*KRK[0]/2.0)-(dZ[9]+ h/2.0)/Kvect.at(4)));
        KRK[2] = h*(-Kvect.at(15)*((dZ[8] + h*KRK[1]/2.0)-(dZ[9]+ h/2.0)/Kvect.at(4)));
        KRK[3] = h*(-Kvect.at(15)*((dZ[8] + h*KRK[2])-(dZ[9]+ h)/Kvect.at(4)));
        dI1   = dZ[8] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dId через Рунге-Кутта
        // dId=-Kvect.at(15)*(Id-I1);
        // dId = f(I1,Id)?

        KRK[0] = h*(-Kvect.at(15)*(dZ[1]-dZ[8]));
        KRK[1] = h*(-Kvect.at(15)*((dZ[1] + h*KRK[0]/2.0)-(dZ[8]+h/2.0)));
        KRK[2] = h*(-Kvect.at(15)*((dZ[1] + h*KRK[1]/2.0)-(dZ[8]+h/2.0)));
        KRK[3] = h*(-Kvect.at(15)*((dZ[1] + h*KRK[2])-(dZ[8]+h)));
        dId   = dZ[1] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dXt через Рунге-Кутта
        // dXt=-Kvect.at(30)*Xt+Kvect.at(30)*((Ip/Kvect.at(4))-Kvect.at(1))*Heavi2;
        // dXt = f(Ip,Xt)?

        double Heavi2;
        if( ((dZ[9]/Kvect.at(4))-Kvect.at(1)) > 0){
             Heavi2 = 1;}
        if( ((dZ[9]/Kvect.at(4))-Kvect.at(1)) < 0){
             Heavi2 = 0;}
        if( ((dZ[9]/Kvect.at(4))-Kvect.at(1)) == 0){
             Heavi2 = 0.5;}

        KRK[0] = h*(-Kvect.at(30)*dZ[7]+Kvect.at(30)*((dZ[9]/Kvect.at(4))-Kvect.at(1))*Heavi2);

            if( (((dZ[9]+ h/2.0)/Kvect.at(4))-Kvect.at(1)) > 0){
                 Heavi2 = 1;}
            if( (((dZ[9]+ h/2.0)/Kvect.at(4))-Kvect.at(1)) < 0){
                 Heavi2 = 0;}
            if( (((dZ[9]+ h/2.0)/Kvect.at(4))-Kvect.at(1)) == 0){
                 Heavi2 = 0.5;}

        KRK[1] = h*(-Kvect.at(30)*(dZ[7]+ h*KRK[0]/2.0)+Kvect.at(30)*(((dZ[9]+ h/2.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

            if( (((dZ[9]+ h/2.0)/Kvect.at(4))-Kvect.at(1)) > 0){
                 Heavi2 = 1;}
            if( (((dZ[9]+ h/2.0)/Kvect.at(4))-Kvect.at(1)) < 0){
                 Heavi2 = 0;}
            if( (((dZ[9]+ h/2.0)/Kvect.at(4))-Kvect.at(1)) == 0){
                 Heavi2 = 0.5;}

        KRK[2] = h*(-Kvect.at(30)*(dZ[7] + h*KRK[1]/2.0)+Kvect.at(30)*(((dZ[9]+ h/2.0)/Kvect.at(4))-Kvect.at(1))*Heavi2);

            if( (((dZ[9]+ h)/Kvect.at(4))-Kvect.at(1)) > 0){
                 Heavi2 = 1;}
            if( (((dZ[9]+ h)/Kvect.at(4))-Kvect.at(1)) < 0){
                 Heavi2 = 0;}
            if( (((dZ[9]+ h)/Kvect.at(4))-Kvect.at(1)) == 0){
                 Heavi2 = 0.5;}

        KRK[3] = h*(-Kvect.at(30)*(dZ[7]+ h*KRK[2])+Kvect.at(30)*(((dZ[9]+ h)/Kvect.at(4))-Kvect.at(1))*Heavi2);

        dXt   = dZ[7] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dIpo через Рунге-Кутта
        // dIpo=-Kvect.at(33)*Ipo+(Yt+Kvect.at(3))*Heavi5+(Yt+Kvect.at(3))*(Heavi6);
        // dIpo = f(Yt,dgp,Ipo)?

        double Heavi5;
        double Heavi6;

        if( (dgp) > 0){
             Heavi5 = 1;}
        if( (dgp) < 0){
             Heavi5 = 0;}
        if( (dgp) == 0){
             Heavi5 = 0.5;}

        if( (-dgp) > 0){
             Heavi6 = 1;}
        if( (-dgp) < 0){
             Heavi6 = 0;}
        if( (-dgp) == 0){
             Heavi6 = 0.5;}

        KRK[0] = h*(-Kvect.at(33)*dZ[2]+(dZ[12]+Kvect.at(3))*Heavi5+(dZ[12]+Kvect.at(3))*(Heavi6));

        KRK[1] = h*(-Kvect.at(33)*(dZ[2]+ h*KRK[0]/2.0)+((dZ[12]+ h/2.0)+Kvect.at(3))*Heavi5+((dZ[12]+ h/2.0)+Kvect.at(3))*(Heavi6));

        KRK[2] = h*(-Kvect.at(33)*(dZ[2]+ h*KRK[1]/2.0)+((dZ[12]+ h/2.0)+Kvect.at(3))*Heavi5+((dZ[12]+ h/2.0)+Kvect.at(3))*(Heavi6));

        KRK[3] = h*(-Kvect.at(33)*(dZ[2]+ h*KRK[2])+((dZ[12]+ h)+Kvect.at(3))*Heavi5+((dZ[12]+ h)+Kvect.at(3))*(Heavi6));

        dIpo   = dZ[2] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dYt через Рунге-Кутта
        // dYt=-Kvect.at(34)*(Yt-Kvect.at(35)*(gp/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*Yt-Kvect.at(34)*Kvect.at(3))*(Heavi8);
        // dYt = f(gp,Yt)?

        double Heavi7;
        double Heavi8;

        if( (Kvect.at(35)*(dZ[0]/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) > 0){
             Heavi7 = 1;}
        if( (Kvect.at(35)*(dZ[0]/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
             Heavi7 = 0;}
        if( (Kvect.at(35)*(dZ[0]/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) == 0){
             Heavi7 = 0.5;}

        if( (-Kvect.at(3)-Kvect.at(35)*(dZ[0]/Kvect.at(23)-Kvect.at(2))) > 0){
             Heavi8 = 1;}
        if( (-Kvect.at(3)-Kvect.at(35)*(dZ[0]/Kvect.at(23)-Kvect.at(2))) < 0){
             Heavi8 = 0;}
        if( (-Kvect.at(3)-Kvect.at(35)*(dZ[0]/Kvect.at(23)-Kvect.at(2))) == 0){
             Heavi8 = 0.5;}

        KRK[0] = h*(-Kvect.at(34)*(dZ[12]-Kvect.at(35)*(dZ[0]/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*dZ[12]-Kvect.at(34)*Kvect.at(3))*(Heavi8));

            if( (Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) > 0){
                 Heavi7 = 1;}
            if( (Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
                 Heavi7 = 0;}
            if( (Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) == 0){
                 Heavi7 = 0.5;}

            if( (-Kvect.at(3)-Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))) > 0){
                 Heavi8 = 1;}
            if( (-Kvect.at(3)-Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))) < 0){
                 Heavi8 = 0;}
            if( (-Kvect.at(3)-Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))) == 0){
                 Heavi8 = 0.5;}

        KRK[1] = h*(-Kvect.at(34)*((dZ[12]+ h*KRK[0]/2.0)-Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(dZ[12]+ h*KRK[0]/2.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

            if( (Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) > 0){
                 Heavi7 = 1;}
            if( (Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
                 Heavi7 = 0;}
            if( (Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) == 0){
                 Heavi7 = 0.5;}

            if( (-Kvect.at(3)-Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))) > 0){
                 Heavi8 = 1;}
            if( (-Kvect.at(3)-Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))) < 0){
                 Heavi8 = 0;}
            if( (-Kvect.at(3)-Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2))) == 0){
                 Heavi8 = 0.5;}

        KRK[2] = h*(-Kvect.at(34)*((dZ[12]+ h*KRK[1]/2.0)-Kvect.at(35)*((dZ[0]+ h/2.0)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(dZ[12]+ h*KRK[1]/2.0)-Kvect.at(34)*Kvect.at(3))*(Heavi8));

            if( (Kvect.at(35)*((dZ[0]+ h)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) > 0){
                 Heavi7 = 1;}
            if( (Kvect.at(35)*((dZ[0]+ h)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) < 0){
                 Heavi7 = 0;}
            if( (Kvect.at(35)*((dZ[0]+ h)/Kvect.at(23)-Kvect.at(2))+Kvect.at(3)) == 0){
                 Heavi7 = 0.5;}

            if( (-Kvect.at(3)-Kvect.at(35)*((dZ[0]+ h)/Kvect.at(23)-Kvect.at(2))) > 0){
                 Heavi8 = 1;}
            if( (-Kvect.at(3)-Kvect.at(35)*((dZ[0]+ h)/Kvect.at(23)-Kvect.at(2))) < 0){
                 Heavi8 = 0;}
            if( (-Kvect.at(3)-Kvect.at(35)*((dZ[0]+ h)/Kvect.at(23)-Kvect.at(2))) == 0){
                 Heavi8 = 0.5;}

        KRK[3] = h*(-Kvect.at(34)*((dZ[12]+ h*KRK[2])-Kvect.at(35)*((dZ[0]+ h)/Kvect.at(23)-Kvect.at(2)))*Heavi7+(-Kvect.at(34)*(dZ[12]+ h*KRK[2])-Kvect.at(34)*Kvect.at(3))*(Heavi8));

        dYt   = dZ[12] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dIi через Рунге-Кутта
        // dIi=-Kvect.at(8)*Ii+vbas+vbol;
        // dIi = f(vbol,Ii,vbas)?
        KRK[0] = h*(-Kvect.at(8)*dZ[13]+vbas+vbol);
        KRK[1] = h*(-Kvect.at(8)*(dZ[13]+ h*KRK[0]/2.0)+(vbas)+(vbol));
        KRK[2] = h*(-Kvect.at(8)*(dZ[13]+ h*KRK[1]/2.0)+(vbas)+(vbol));
        KRK[3] = h*(-Kvect.at(8)*(dZ[13]+ h*KRK[2])+(vbas)+(vbol));

        dIi   = dZ[13] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        // Решаем dIt через Рунге-Кутта
        // dIt=Kvect.at(8)*Ii-Kvect.at(10)*It;
        // dIt = f(Ii,It)?
        KRK[0] = h*(Kvect.at(8)*dZ[13]-Kvect.at(10)*dZ[10]);
        KRK[1] = h*(Kvect.at(8)*(dZ[13]+ h/2.0)-Kvect.at(10)*(dZ[10]+ h*KRK[0]/2.0));
        KRK[2] = h*(Kvect.at(8)*(dZ[13]+ h/2.0)-Kvect.at(10)*(dZ[10]+ h*KRK[1]/2.0));
        KRK[3] = h*(Kvect.at(8)*(dZ[13]+ h)-Kvect.at(10)*(dZ[10]+ h*KRK[2]));

        dIt   = dZ[10] + (KRK[0] + 2.0*KRK[1] + 2.0*KRK[2] + KRK[3])/6.0;

        dZ[0] = dgp;
        dZ[1] = dId;
        dZ[2] = dIpo;
        dZ[3] = dfgut;
        dZ[4] = dfsol;
        dZ[5] = dfliq;
        dZ[6] = dgt;
        dZ[7] = dXt;
        dZ[8] = dI1;
        dZ[9] = dIp;
        dZ[10] = dIt;
        dZ[11] = dIl;
        dZ[12] = dYt;
        dZ[13] = dIi;

        tick.append(t);
        CGB.append(dZ[0]/Kvect.at(23));
        Ipg.append(dZ[9]*20);

        t = t + h; // увеличиваем время
        j = j + h;
    }
}
