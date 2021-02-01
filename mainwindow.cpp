#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <vector>
#include <QVector>
#include <math.h>
#include "qcustomplot.h"


using namespace std;
                            // начальные значения ???
double t = 0.0;         // время сейчас?

    /* meal data */
double tm1 = 30; // tm  // время приёма пищи
double Tm = 15;         // длительность приёма пищи
double Dig = 1e-05;     // масса углеводов начальная
    /* bolus */
double Ti1 = 10;        // наверное?
double Ti2 = 10;        // что из этого длительность ввода а что время ввода
double ti2 = 20;
double ti1 = 20;
double Dbol = 1.8;      // доза болюса начальная
    /* bazal ins */
double Vbas = 0.15;     // доза базального инсулина

double b = 0.72;        // percentage, kgut=(kmin+kmax)/2 right
double c = 0.115;       //percentage, kgut=(kmin+kmax)/2 left обеспечивает 2 волны

double Vg = 1.8;        // plasma per BW, dl/kg
double Vmxx =0.047;     // mg/kg/min per pmol/l   insulin sens
double kp1 = 3.27;      // mg/kg/min

    /* Glucose Utilization */
double Fcns = 0.8;      // mg/kg/min

double m2 = 0.001;      // неивестная переменная c неизвестным значением
double del = 0.001;     // неивестная переменная c неизвестным значением
double tms = 0.001;     // неивестная переменная c неизвестным значением
double mt = 0.001;      // неивестная переменная c неизвестным значением

/* Разные названия переменных с начальным массивом Z ? */
double g0 = 90;      // gp?
double gt0 = 90;     // можно ли их все сразу заменить на те
double Ip0 = 0.20;     // что в структуре Z ?

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // структура объявляется в заголовочном


    double  a = 0.0,    // начальное время
            b = 1440.0;    // конечное время
    double  h = 5;   // шаг

    vector<DataZ> data;
    DataZ start;
    start.t = a;

    int stepCount = (b / h) + 1;
    double currentTime = a + h;

    // double K[4];
    double K[7];

    Dig = (1e-05)*1000; // вместо 1е-05 массу углеводов

     tm1 = 30;
     Tm = 30;
     tms = 0.0001;

    /************************************************                        Дано                                           *************************************************/

    /* bgDynam */               // присвоение стартовых значений в структуру для расчёта
    double Ginit=130.22;        // начальный уровень гликемии
    start.gp=Ginit*1.8;
    start.Il=2.478;
    start.Ip=2.756; //55.12/20;
    start.fgut=0; // mg
    start.fliq=0; // mg
    start.fsol=0; // mg
    start.gt=116.8;
    start.I1=55.12;
    start.Id=55.12;
    start.Xt=29.56;
    start.Ipo=0;
    start.Yt=-0.6926;
    start.Ii=584.6;
    start.It=1796;
    data.push_back(start);

   // DVI=vbasal*115.75; //u/h->pmol/min

    // Initial
    double Sb=1.55; // pmol/kg/min
    double EGPb=1.92; // mg/kg/min
    double Gpb;
    double Gtb;
    double Ipb;
    double Ib;
    double Gb;
    // Initial
    Gpb=g0; // mg/kg
    Gtb=gt0; // mg/kg
    Ipb=Ip0; // pmol/kg
    Ib=Ipb*20; // pmol/l
    Gb=Gpb/1.8; // mg/dl



    // Insulin
    double V1=0.05; //l/kg
    double m1=0.190; //min^-1
    double m6=0.6471; //dimensionless
    double m3t;
    // Insulin
    m3t=(m6*m1)/(1- m6);

    // Insulin Kraegen
    double k21=2.97e-2; // КС всасывания инсулина в ткани, min-1 Cobelli
    double di=12e-2; // КС деградации инсулина, min-1
    double ka=11.3e-3; // КС всасывания инсулина в плазму, min-1

    // Endogenous glucose production
    // kp1=2.70; %mg/kg/min
    double kp2=0.0021; //min^-1
    double kp3=0.009; //mg/kg/min per pmol/l
    double kp4=0.0618; //mg/kg/min per pmol/kg
    double ki=0.0079; //min^-1

    // Glucose appearance
    double fract=0.9; // Assimilated glucose fraction
    double kgabs=0.057; // RC glucose absorption, min-1,  Dalla Man
    double kgri=0.056; // RC grinding, min-1,  Dalla Man
    double kmin=0.008; // min rate of gut empying
    double kmax=0.056; // max rate of gut empying
    // c=0.115; // percentage, kgut=(kmin+kmax)/2 left обеспечивает 2 волны
    // b=0.72; // percentage, kgut=(kmin+kmax)/2 right
    Vg=1.8; // plasma per BW, dl/kg

    // Glucose compartments
    double k1gg=0.065;
    double k2gg=0.079;

    // Glucose Utilization
    // Fcns=0.8; %mg/kg/min
    double Vm0=2.0; //mg/kg/min
    // Vmxx=0.047; %mg/kg/min per pmol/l
    double Km0=205.59; //mg/kg
    double p2U=0.0731; //min^-1

    // Renal scretion
    double ke1=0.0005; // Rate of renal excretion, min-1, DM
    double ke2=339/Vg; // Threshold of renal excretion, mg/dl

    // Secretion
    double gamma=0.5;  // min^-1
    double alpha=0.050; // min^-1
    double betha=0.11; // pmol/kg per (mg/dl)

    // Fluctuations
    double vbas;
    double Abol;
    double bol1;
    double bol2;
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

    double gp;
    gp = start.gp;
    double Id;
    Id = start.Id;
    double Ipo;
    Ipo = start.Ipo;
    double fgut;
    fgut = start.fgut;
    double fsol;
    fsol = start.fsol;
    double fliq;
    fliq = start.fliq;
    double gt;
    gt = start.gt;
    double Xt;
    Xt = start.Xt;
    double I1;
    I1 = start.I1;
    double Ip;
    Ip = start.Ip;
    double It;
    It = start.It;
    double Il;
    Il = start.Il;
    double Yt;
    Yt = start.Yt;
    double Ii;
    Ii = start.Ii;
    /*************************************************                      Конец дано)                              ********************************************************/

    for (int j = 1; j <= stepCount; j++)
    {

        /* */
        t = currentTime;

        // Fluctuations
        vbas=(Vbas+0.00001)*6000/60; // pmol/min
        Abol=(Dbol+0.001)*6000;
        //точки из матлаб убрал в последующих двух уравнениях после второй скобки множителя
        bol1=1/Ti1*(1-del)*(1./(1+exp(-3*(t-tm1+10-ti1))))*(1./(1+exp(-3*(tm1-10+ti1-t+Ti1))));
        bol2=1/Ti2*del*(1./(1+exp(-3*(t-tm1+10-ti2))))*(1./(1+exp(-3*(tm1-10+ti2-t+Ti2))));
        vbol=Abol*(bol1+bol2);

        double Heavi1;
        if( (t-tm1-tms) >= 0){
             Heavi1 = 1;}
        if( (t-tm1-tms) < 0){
             Heavi1 = 0;}

        double Heavi2;
        if( (-(t-tm1-Tm-tms)) >= 0){
             Heavi2 = 1;}
        if( (-(t-tm1-Tm-tms)) < 0){
             Heavi2 = 0;}
        vm=Dig/Tm*Heavi1*Heavi2; // mg/min

        //
        double alp=5/(2*Dig*(1-b));
        double bet=5/(2*Dig*c);

        //



        double Heavi3;
        if( (kp1-kp2*gp-kp3*Id - kp4*Ipo)  >= 0){
             Heavi3 = 1;}
        if( (kp1-kp2*gp-kp3*Id-kp4*Ipo) < 0){
             Heavi3 = 0;}
        double EGP=(kp1-kp2*gp-kp3*Id-kp4*Ipo)*Heavi3;


        double Vmx=Vm0+Vmxx*Xt;


        double Uid=(Vmx*gt)/(Km0+gt);
        Gtb=(Fcns-EGPb+k1gg*Gpb)/k2gg;


        double kgut=kmin+(kmax-kmin)/2*(tanh(alp*(fsol+fliq-b*Dig))-tanh(bet*(fsol+fliq-c*Dig))+2);
        double G=gp/Vg;


        double fmeal=fract/mt*kgabs*fgut; // Meal

        double Heavi4;
        if( (gp-ke2)  >= 0){
             Heavi4 = 1;}
        if( (gp-ke2) < 0){
             Heavi4 = 0;}
        double eren=ke1*(gp-ke2)*Heavi4; // Renal glucose excretion

        /* ДУ */

        // h - шаг, h=5

        // Решаем dgp через Рунге-Кутта                                                                  // или оставить прирощение только для главной величины???
        // dgp = EGP+fmeal-Fcns-eren-k1gg*gp+k2gg*gt; // plasma glucose mg/dl
        // dgp = f(gp,gt)?
        K[0] = EGP+fmeal-Fcns-eren-k1gg*gp+k2gg*gt;
        K[1] = EGP+fmeal-Fcns-eren-k1gg*(gp+ h*K[0]/2.0)+k2gg*(gt+ h/2.0);
        K[2] = EGP+fmeal-Fcns-eren-k1gg*(gp+ h*K[1]/2.0)+k2gg*(gt+ h/2.0);
        K[3] = EGP+fmeal-Fcns-eren-k1gg*(gp+ h*K[2])+k2gg*gt;
        dgp   = gp + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dgp через Дормана-Принса
/*
        K[0] = EGP+fmeal-Fcns-eren-k1gg*gp+k2gg*gt;
        K[1] = EGP+fmeal-Fcns-eren-k1gg*(gp+ h*K[0]/5.0)+k2gg*(gt+ h/5.0);
        K[2] = EGP+fmeal-Fcns-eren-k1gg*(gp+ (3*h*K[0]/40.0) + (9*h*K[1]/40.0) )+k2gg*(gt+ 3*h/10.0);
        K[3] = EGP+fmeal-Fcns-eren-k1gg*(gp+ (44*h*K[0]/45.0) + (-56*h*K[1]/15.0) + (32*h*K[2]/9.0) )+k2gg*(gt+ 4*h/5.0);
        K[4] = EGP+fmeal-Fcns-eren-k1gg*(gp+ (19372*h*K[0]/6561.0) + (-25360*h*K[1]/2187.0) + (64448*h*K[2]/6561.0) + (-212*h*K[3]/729.0) )+k2gg*(gt+ 8*h/9.0);
        K[5] = EGP+fmeal-Fcns-eren-k1gg*(gp+ (9017*h*K[0]/3168.0) + (-355*h*K[1]/33.0) + (46732*h*K[2]/5247.0) + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )+k2gg*(gt+ h);
        K[6] = EGP+fmeal-Fcns-eren-k1gg*(gp+ (35*h*K[0]/384.0) + (500*h*K[2]/1113.0) + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + (11*h*K[5]/84.0) )+k2gg*(gt+ h);

        dgp   = gp + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = gp + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dgt через Рунге-Кутта
        // dgt=-Uid+k1gg*gp-k2gg*gt;
        // dgt = f(gp,gt)?
        K[0] = -Uid+k1gg*gp-k2gg*gt;
        K[1] = -Uid+k1gg*(gp+ h/2.0)-k2gg*(gt+ h*K[0]/2.0);
        K[2] = -Uid+k1gg*(gp+ h/2.0)-k2gg*(gt+ h*K[1]/2.0);
        K[3] = -Uid+k1gg*(gp+ h)-k2gg*(gt+ h*K[2]);
        dgt   = gt + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dgt через Дормана-Принса
/*
        K[0] = -Uid+k1gg*gp-k2gg*gt;
        K[1] = -Uid+k1gg*(gp+ h/5.0)-k2gg*(gt+ h*K[0]/5.0);
        K[2] = -Uid+k1gg*(gp+ 3*h/10.0)-k2gg*(gt+ 3*h*K[0]/40.0 +9*h*K[1]/40.0);
        K[3] = -Uid+k1gg*(gp+ 4*h/5.0)-k2gg*(gt+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0);
        K[4] = -Uid+k1gg*(gp+ 8*h/9.0)-k2gg*(gt+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) );
        K[5] = -Uid+k1gg*(gp+ h)-k2gg*(gt+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) );
        K[6] = -Uid+k1gg*(gp+ h)-k2gg*(gt+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0 );

        dgt   = gt + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = gt + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dI1 через Рунге-Кутта
        // dI1=-ki*(I1-Ip/V1);
        // dI1 = f(Id,I1)?
        K[0] = -ki*(Id-I1);
        K[1] = -ki*((Id+ h/2.0)-(I1+ h*K[0]/2.0));
        K[2] = -ki*((Id+ h/2.0)-(I1+ h*K[1]/2.0));
        K[3] = -ki*((Id+ h)-(I1+ h*K[2]));
        dI1   = I1 + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dI1 через Дормана-Принса
/*
        K[0] = -ki*(Id-I1);
        K[1] = -ki*((Id+ h/5.0)-(I1+ h*K[0]/5.0));
        K[2] = -ki*((Id+ 3*h/10.0)-(I1+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0));
        K[3] = -ki*((Id+ 4*h/5.0)-(I1+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0));
        K[4] = -ki*((Id+ 8*h/9.0)-(I1+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) ));
        K[5] = -ki*((Id+ h)-(I1+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) ));
        K[6] = -ki*((Id+ h)-(I1+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0));

        dI1   = I1 + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = I1 + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dId через Рунге-Кутта
        // dId=-ki*(Id-I1);
        // dId = f(Id,I1)?
        K[0] = -ki*(Id-I1);
        K[1] = -ki*((Id+ h*K[0]/2.0)-(I1+ h/2.0));
        K[2] = -ki*((Id+ h*K[1]/2.0)-(I1+ h/2.0));
        K[3] = -ki*((Id+ h*K[2])-(I1+ h));
        dId   = Id + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dId через Дормана-Принса
/*
        K[0] = -ki*(Id-I1);
        K[1] = -ki*((Id+ h*K[0]/5.0)-(I1+ h/5.0));
        K[2] = -ki*((Id+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)-(I1+ 3*h/10.0));
        K[3] = -ki*((Id+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0 )-(I1+ 4*h/5.0));
        K[4] = -ki*((Id+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )-(I1+ 8*h/9.0));
        K[5] = -ki*((Id+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (-49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )-(I1+ h));
        K[6] = -ki*((Id+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0)-(I1+ h));

        dId   = Id + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = Id + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dXt через Рунге-Кутта
        // dXt=-p2U*Xt+p2U*((Ip/V1)-Ib)*heaviside((Ip/V1)-Ib);
        // dXt = f(Xt,Ip)?
        double Heavi5;
        if( ((Ip/V1)-Ib)  >= 0){
             Heavi5 = 1;}
        if( ((Ip/V1)-Ib) < 0){
             Heavi5 = 0;}
        K[0] = -p2U*Xt+p2U*((Ip/V1)-Ib)*Heavi5;

        double Heavi6;
        if( (((Ip+ h/2.0)/V1)-Ib)  >= 0){
             Heavi6 = 1;}
        if( (((Ip+ h/2.0)/V1)-Ib) < 0){
             Heavi6 = 0;}
        K[1] = -p2U*(Xt+ h*K[0]/2.0)+p2U*(((Ip+ h/2.0)/V1)-Ib)*Heavi6;

        double Heavi7;
        if( (((Ip+ h/2.0)/V1)-Ib)  >= 0){
             Heavi7 = 1;}
        if( (((Ip+ h/2.0)/V1)-Ib) < 0){
             Heavi7 = 0;}
        K[2] = -p2U*(Xt+ h*K[1]/2.0)+p2U*(((Ip+ h/2.0)/V1)-Ib)*Heavi7;

        double Heavi8;
        if( (((Ip+ h)/V1)-Ib)  >= 0){
             Heavi8 = 1;}
        if( (((Ip+ h)/V1)-Ib) < 0){
             Heavi8 = 0;}
        K[3] = -p2U*(Xt+ h*K[2])+p2U*(((Ip+ h)/V1)-Ib)*Heavi8;
        dXt   = Xt + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dXt через Дормана-Принса
/*
        double Heavi13;
        if( ((Ip/V1)-Ib)  >= 0){
             Heavi13 = 1;}
        if( ((Ip/V1)-Ib) < 0){
             Heavi13 = 0;}
        K[0] = -p2U*Xt+p2U*((Ip/V1)-Ib)*Heavi13;
        double Heavi14;
        if( (((Ip+ h/5.0)/V1)-Ib)  >= 0){
             Heavi14 = 1;}
        if( (((Ip+ h/5.0)/V1)-Ib) < 0){
             Heavi14 = 0;}
        K[1] = -p2U*(Xt+ h*K[0]/5.0)+p2U*(((Ip+ h/5.0)/V1)-Ib)*Heavi14;
        double Heavi15;
        if( (((Ip+ 3*h/10.0)/V1)-Ib)  >= 0){
             Heavi15 = 1;}
        if( (((Ip+ 3*h/10.0)/V1)-Ib) < 0){
             Heavi15 = 0;}
        K[2] = -p2U*(Xt+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)+p2U*(((Ip+ 3*h/10.0)/V1)-Ib)*Heavi15;
        double Heavi16;
        if( (((Ip+ 4*h/5.0)/V1)-Ib)  >= 0){
             Heavi16 = 1;}
        if( (((Ip+ 4*h/5.0)/V1)-Ib) < 0){
             Heavi16 = 0;}
        K[3] = -p2U*(Xt+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)+p2U*(((Ip+ 4*h/5.0)/V1)-Ib)*Heavi16;
        double Heavi17;
        if( (((Ip+ 8*h/9.0)/V1)-Ib)  >= 0){
             Heavi17 = 1;}
        if( (((Ip+ 8*h/9.0)/V1)-Ib) < 0){
             Heavi17 = 0;}
        K[4] = -p2U*(Xt+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )+p2U*(((Ip+ 8*h/9.0)/V1)-Ib)*Heavi17;
        double Heavi18;
        if( (((Ip+ h)/V1)-Ib)  >= 0){
             Heavi18 = 1;}
        if( (((Ip+ h)/V1)-Ib) < 0){
             Heavi18 = 0;}
        K[5] = -p2U*(Xt+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )+p2U*(((Ip+ h)/V1)-Ib)*Heavi18;
        double Heavi19;
        if( (((Ip+ h)/V1)-Ib)  >= 0){
             Heavi19 = 1;}
        if( (((Ip+ h)/V1)-Ib) < 0){
             Heavi19 = 0;}
        K[6] = -p2U*(Xt+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0)+p2U*(((Ip+ h)/V1)-Ib)*Heavi19;

        dXt   = Xt + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = Xt + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dIl через Рунге-Кутта
        // dIl=-(m1+m3t)*Il+m2*Ip; // Liver insulin pmol/kg
        // dIl = f(Il,Ip)?
        K[0] = -(m1+m3t)*Il+m2*Ip;
        K[1] = -(m1+m3t)*(Il+ h*K[0]/2.0)+m2*(Ip+ h/2.0);
        K[2] = -(m1+m3t)*(Il+ h*K[1]/2.0)+m2*(Ip+ h/2.0);
        K[3] = -(m1+m3t)*(Il+ h*K[2])+m2*(Ip+ h);
        dIl   = Il + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

         // Решение dIl через Дормана-Принса
/*
        K[0] = -(m1+m3t)*Il+m2*Ip;
        K[1] = -(m1+m3t)*(Il+ h*K[0]/5.0)+m2*(Ip+ h/5.0);
        K[2] = -(m1+m3t)*(Il+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)+m2*(Ip+ 3*h/10.0);
        K[3] = -(m1+m3t)*(Il+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)+m2*(Ip+ 4*h/5.0);
        K[4] = -(m1+m3t)*(Il+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )+m2*(Ip+ 8*h/9.0);
        K[5] = -(m1+m3t)*(Il+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )+m2*(Ip+ h);
        K[6] = -(m1+m3t)*(Il+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0)+m2*(Ip+ h);

        dIl   = Il + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = Il + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/


        // Решаем dIp через Рунге-Кутта
        // dIp=-m2*Ip+m1*Il+ka/mt*It-di*Ip; // Plasma insulin pmol/kg
        // dIp = f(Ip,It)? // f(Ip,Il)? // f(Ip,Il,It)? ВОПРОС ПРО ПЕРЕМЕННЫЕ!
        K[0] = -m2*Ip+m1*Il+ka/mt*It-di*Ip;
        K[1] = -m2*(Ip+ h*K[0]/2.0)+m1*(Il+ h/2.0)+ka/mt*(It+ h/2.0)-di*(Ip+ h*K[0]/2.0);
        K[2] = -m2*(Ip+ h*K[1]/2.0)+m1*(Il+ h/2.0)+ka/mt*(It+ h/2.0)-di*(Ip+ h*K[1]/2.0);
        K[3] = -m2*(Ip+ h*K[2])+m1*(Il+ h)+ka/mt*(It+ h)-di*(Ip+ h*K[2]);
        dIp   = Ip + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIp через Дормана-Принса
/*
        K[0] = -m2*Ip+m1*Il+ka/mt*It-di*Ip;
        K[1] = -m2*(Ip+ h*K[0]/5.0)+m1*(Il+ h/5.0)+ka/mt*(It+ h/5.0)-di*(Ip+ h*K[0]/5.0);
        K[2] = -m2*(Ip+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)+m1*(Il+ 3*h/10.0)+ka/mt*(It+ 3*h/10.0)-di*(Ip+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0);
        K[3] = -m2*(Ip+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)+m1*(Il+ 4*h/5.0)+ka/mt*(It+ 4*h/5.0)-di*(Ip+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0);
        K[4] = -m2*(Ip+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )+m1*(Il+ 8*h/9.0)+ka/mt*(It+ 8*h/9.0)-di*(Ip+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) );
        K[5] = -m2*(Ip+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )+m1*(Il+ h)+ka/mt*(It+ h)-di*(Ip+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) );
        K[6] = -m2*(Ip+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0 )+m1*(Il+ h)+ka/mt*(It+ h)-di*(Ip+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0 );

        dIp   = Ip + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = Ip + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dfgut через Рунге-Кутта
        // dfgut=-kgabs*fgut+kgut*fliq; // Gut glucose
        // dfgut = f(fgut,fliq)?
        K[0] = -kgabs*fgut+kgut*fliq;
        K[1] = -kgabs*(fgut+ h*K[0]/2.0)+kgut*(fliq+ h/2.0);
        K[2] = -kgabs*(fgut+ h*K[1]/2.0)+kgut*(fliq+ h/2.0);
        K[3] = -kgabs*(fgut+ h*K[2])+kgut*(fliq+ h);
        dfgut   = fgut + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfgut через Дормана-Принса
/*
        K[0] = -kgabs*fgut+kgut*fliq;
        K[1] = -kgabs*(fgut+ h*K[0]/5.0)+kgut*(fliq+ h/5.0);
        K[2] = -kgabs*(fgut+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)+kgut*(fliq+ 3*h/10.0);
        K[3] = -kgabs*(fgut+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)+kgut*(fliq+ 4*h/5.0);
        K[4] = -kgabs*(fgut+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )+kgut*(fliq+ 8*h/9.0);
        K[5] = -kgabs*(fgut+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )+kgut*(fliq+ h);
        K[6] = -kgabs*(fgut+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0)+kgut*(fliq+ h);

        dfgut   = fgut + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = fgut + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dfliq через Рунге-Кутта
        // dfliq=-kgut*fliq+kgri*fsol; // Glucose in liquid phase
        // dfliq = f(fliq,fsol)?
        K[0] = -kgut*fliq+kgri*fsol;
        K[1] = -kgut*(fliq+ h*K[0]/2.0)+kgri*(fsol+ h/2.0);
        K[2] = -kgut*(fliq+ h*K[1]/2.0)+kgri*(fsol+ h/2.0);
        K[3] = -kgut*(fliq+ h*K[2])+kgri*(fsol+ h);
        dfliq   = fliq + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

         // Решение dfliq через Дормана-Принса
/*
        K[0] = -kgut*fliq+kgri*fsol;
        K[1] = -kgut*(fliq+ h*K[0]/5.0)+kgri*(fsol+ h/5.0);
        K[2] = -kgut*(fliq+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)+kgri*(fsol+ 3*h/10.0);
        K[3] = -kgut*(fliq+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)+kgri*(fsol+ 4*h/5.0);
        K[4] = -kgut*(fliq+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )+kgri*(fsol+ 8*h/9.0);
        K[5] = -kgut*(fliq+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )+kgri*(fsol+ h);
        K[6] = -kgut*(fliq+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0)+kgri*(fsol+ h);

        dfliq   = fliq + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = fliq + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dfsol через Рунге-Кутта
        // dfsol=-kgri*fsol+vm; // Glucose in solid phase
        // dfsol = f(fsol,vm)?
        K[0] = -kgri*fsol+vm;
        K[1] = -kgri*(fsol+ h*K[0]/2.0)+(vm+ h/2.0);
        K[2] = -kgri*(fsol+ h*K[1]/2.0)+(vm+ h/2.0);
        K[3] = -kgri*(fsol+ h*K[2])+(vm+ h);
        dfsol   = fsol + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dfsol через Дормана-Принса
/*
        K[0] = -kgri*fsol+vm;
        K[1] = -kgri*(fsol+ h*K[0]/5.0)+(vm+ h/5.0);
        K[2] = -kgri*(fsol+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)+(vm+ 3*h/10.0);
        K[3] = -kgri*(fsol+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)+(vm+ 4*h/5.0);
        K[4] = -kgri*(fsol+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )+(vm+ 8*h/9.0);
        K[5] = -kgri*(fsol+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) +(-5103*h*K[4]/18656.0) )+(vm+ h);
        K[6] = -kgri*(fsol+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) +(-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0)+(vm+ h);

        dfsol   = fsol + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = fsol + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dIpo через Рунге-Кутта
        // dIpo=-gamma*Ipo+(Yt+Sb)*heaviside(dgp/Vg) + (Yt+Sb)*(heaviside(-dgp/Vg));
        // dIpo = f(Ipo,Yt)?
        double Heavi9;
        if( (dgp/Vg)  >= 0){
             Heavi9 = 1;}
        if( (dgp/Vg) < 0){
             Heavi9 = 0;}
        double Heavi10;
        if( (-dgp/Vg)  >= 0){
             Heavi10 = 1;}
        if( (-dgp/Vg) < 0){
             Heavi10 = 0;}
        K[0] = -gamma*Ipo+(Yt+Sb)*Heavi9 + (Yt+Sb)*(Heavi10);
        K[1] = -gamma*(Ipo+ h*K[0]/2.0)+((Yt+ h/2.0)+Sb)*Heavi9 + ((Yt+ h/2.0)+Sb)*(Heavi10);
        K[2] = -gamma*(Ipo+ h*K[1]/2.0)+((Yt+ h/2.0)+Sb)*Heavi9 + ((Yt+ h/2.0)+Sb)*(Heavi10);
        K[3] = -gamma*(Ipo+ h*K[2])+((Yt+ h)+Sb)*Heavi9 + ((Yt+ h)+Sb)*(Heavi10);
        dIpo   = Ipo + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIpo через Дормана-Принса
/*
        double Heavi20;
        if( (dgp/Vg)  >= 0){
             Heavi20 = 1;}
        if( (dgp/Vg) < 0){
             Heavi20 = 0;}
        double Heavi21;
        if( (-dgp/Vg)  >= 0){
             Heavi21 = 1;}
        if( (-dgp/Vg) < 0){
             Heavi21 = 0;}
        K[0] = -gamma*Ipo+(Yt+Sb)*Heavi20 + (Yt+Sb)*(Heavi21);
        K[1] = -gamma*(Ipo+ h*K[0]/5.0)+((Yt+ h/5.0)+Sb)*Heavi20 + ((Yt+ h/5.0)+Sb)*(Heavi21);
        K[2] = -gamma*(Ipo+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)+((Yt+ 3*h/10.0)+Sb)*Heavi20 + ((Yt+ 3*h/10.0)+Sb)*(Heavi21);
        K[3] = -gamma*(Ipo+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)+((Yt+ 4*h/5.0)+Sb)*Heavi20 + ((Yt+ 4*h/5.0)+Sb)*(Heavi21);
        K[4] = -gamma*(Ipo+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )+((Yt+ 8*h/9.0)+Sb)*Heavi20 + ((Yt+ 8*h/9.0)+Sb)*(Heavi21);
        K[5] = -gamma*(Ipo+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )+((Yt+ h)+Sb)*Heavi20 + ((Yt+ h)+Sb)*(Heavi21);
        K[6] = -gamma*(Ipo+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0)+((Yt+ h)+Sb)*Heavi20 + ((Yt+ h)+Sb)*(Heavi21);

        dIpo   = Ipo + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = Ipo + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/


        // Решаем dYt через Рунге-Кутта
        // dYt=-alpha*(Yt-betha*(G-Gb))*heaviside(betha*(G-Gb)+Sb)+(-alpha*Yt-alpha*Sb)*(heaviside(-Sb-betha*(G-Gb)));
        // dYt = f(Yt,G)?
        double Heavi11;
        if( (betha*(G-Gb)+Sb)  >= 0){
             Heavi11 = 1;}
        if( (betha*(G-Gb)+Sb) < 0){
             Heavi11 = 0;}
        double Heavi12;
        if( (-Sb-betha*(G-Gb))  >= 0){
             Heavi12 = 1;}
        if( (-Sb-betha*(G-Gb)) < 0){
             Heavi12 = 0;}
        K[0] = -alpha*(Yt-betha*(G-Gb))*Heavi11+(-alpha*Yt-alpha*Sb)*(Heavi12);
        double Heavi22;
        if( (betha*((G+h/2.0)-Gb)+Sb)  >= 0){
             Heavi22 = 1;}
        if( (betha*((G+h/2.0)-Gb)+Sb) < 0){
             Heavi22 = 0;}
        double Heavi23;
        if( (-Sb-betha*((G+h/2.0)-Gb))  >= 0){
             Heavi23 = 1;}
        if( (-Sb-betha*((G+h/2.0)-Gb)) < 0){
             Heavi23 = 0;}
        K[1] = -alpha*((Yt+ h*K[0]/2.0)-betha*((G+h/2.0)-Gb))*Heavi22+(-alpha*(Yt+ h*K[0]/2.0)-alpha*Sb)*(Heavi23);
        double Heavi24;
        if( (betha*((G+h/2.0)-Gb)+Sb)  >= 0){
             Heavi24 = 1;}
        if( (betha*((G+h/2.0)-Gb)+Sb) < 0){
             Heavi24 = 0;}
        double Heavi25;
        if( (-Sb-betha*((G+h/2.0)-Gb))  >= 0){
             Heavi25 = 1;}
        if( (-Sb-betha*((G+h/2.0)-Gb)) < 0){
             Heavi25 = 0;}
        K[2] = -alpha*((Yt+ h*K[1]/2.0)-betha*((G+ h/2.0)-Gb))*Heavi24+(-alpha*(Yt+ h*K[1]/2.0)-alpha*Sb)*(Heavi25);
        double Heavi26;
        if( (betha*((G+h)-Gb)+Sb)  >= 0){
             Heavi26 = 1;}
        if( (betha*((G+h)-Gb)+Sb) < 0){
             Heavi26 = 0;}
        double Heavi27;
        if( (-Sb-betha*((G+h)-Gb))  >= 0){
             Heavi27 = 1;}
        if( (-Sb-betha*((G+h)-Gb)) < 0){
             Heavi27 = 0;}
        K[3] = -alpha*((Yt+ h*K[2])-betha*((G+h)-Gb))*Heavi26+(-alpha*(Yt+ h*K[2])-alpha*Sb)*(Heavi27);
        dYt   = Yt + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dYt через Дормана-Принса
/*
        double Heavi28;
        if( (betha*(G-Gb)+Sb)  >= 0){
             Heavi28 = 1;}
        if( (betha*(G-Gb)+Sb) < 0){
             Heavi28 = 0;}
        double Heavi29;
        if( (-Sb-betha*(G-Gb))  >= 0){
             Heavi29 = 1;}
        if( (-Sb-betha*(G-Gb)) < 0){
             Heavi29 = 0;}
        K[0] = -alpha*(Yt-betha*(G-Gb))*Heavi28+(-alpha*Yt-alpha*Sb)*(Heavi29);
        double Heavi30;
        if( (betha*((G+h/5.0)-Gb)+Sb)  >= 0){
             Heavi30 = 1;}
        if( (betha*((G+h/5.0)-Gb)+Sb) < 0){
             Heavi30 = 0;}
        double Heavi31;
        if( (-Sb-betha*((G+h/5.0)-Gb))  >= 0){
             Heavi31 = 1;}
        if( (-Sb-betha*((G+h/5.0)-Gb)) < 0){
             Heavi31 = 0;}
        K[1] = -alpha*((Yt+ h*K[0]/5.0)-betha*((G+h/5.0)-Gb))*Heavi30+(-alpha*(Yt+ h*K[0]/5.0)-alpha*Sb)*(Heavi31);
        double Heavi32;
        if( (betha*((G+3*h/10.0)-Gb)+Sb)  >= 0){
             Heavi32 = 1;}
        if( (betha*((G+3*h/10.0)-Gb)+Sb) < 0){
             Heavi32 = 0;}
        double Heavi33;
        if( (-Sb-betha*((G+3*h/10.0)-Gb))  >= 0){
             Heavi33 = 1;}
        if( (-Sb-betha*((G+3*h/10.0)-Gb)) < 0){
             Heavi33 = 0;}
        K[2] = -alpha*((Yt+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)-betha*((G+3*h/10.0)-Gb))*Heavi32+(-alpha*(Yt+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)-alpha*Sb)*(Heavi33);
        double Heavi34;
        if( (betha*((G+4*h/5.0)-Gb)+Sb)  >= 0){
             Heavi34 = 1;}
        if( (betha*((G+4*h/5.0)-Gb)+Sb) < 0){
             Heavi34 = 0;}
        double Heavi35;
        if( (-Sb-betha*((G+4*h/5.0)-Gb))  >= 0){
             Heavi35 = 1;}
        if( (-Sb-betha*((G+4*h/5.0)-Gb)) < 0){
             Heavi35 = 0;}
        K[3] = -alpha*((Yt+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)-betha*((G+4*h/5.0)-Gb))*Heavi34+(-alpha*(Yt+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)-alpha*Sb)*(Heavi35);
        double Heavi36;
        if( (betha*((G+8*h/9.0)-Gb)+Sb)  >= 0){
             Heavi36 = 1;}
        if( (betha*((G+8*h/9.0)-Gb)+Sb) < 0){
             Heavi36 = 0;}
        double Heavi37;
        if( (-Sb-betha*((G+8*h/9.0)-Gb))  >= 0){
             Heavi37 = 1;}
        if( (-Sb-betha*((G+8*h/9.0)-Gb)) < 0){
             Heavi37 = 0;}
        K[4] = -alpha*((Yt+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )-betha*((G+8*h/9.0)-Gb))*Heavi36+(-alpha*(Yt+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )-alpha*Sb)*(Heavi37);
        double Heavi38;
        if( (betha*((G+h)-Gb)+Sb)  >= 0){
             Heavi38 = 1;}
        if( (betha*((G+h)-Gb)+Sb) < 0){
             Heavi38 = 0;}
        double Heavi39;
        if( (-Sb-betha*((G+h)-Gb))  >= 0){
             Heavi39 = 1;}
        if( (-Sb-betha*((G+h)-Gb)) < 0){
             Heavi39 = 0;}
        K[5] = -alpha*((Yt+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )-betha*((G+h)-Gb))*Heavi38+(-alpha*(Yt+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )-alpha*Sb)*(Heavi39);
        double Heavi40;
        if( (betha*((G+h)-Gb)+Sb)  >= 0){
             Heavi40 = 1;}
        if( (betha*((G+h)-Gb)+Sb) < 0){
             Heavi40 = 0;}
        double Heavi41;
        if( (-Sb-betha*((G+h)-Gb))  >= 0){
             Heavi41 = 1;}
        if( (-Sb-betha*((G+h)-Gb)) < 0){
             Heavi41 = 0;}
        K[5] = -alpha*((Yt+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) +11*h*K[5]/84.0)-betha*((G+h)-Gb))*Heavi40+(-alpha*(Yt+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) +11*h*K[5]/84.0)-alpha*Sb)*(Heavi41);

        dYt   = Yt + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = Yt + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        // Решаем dIt через Рунге-Кутта
        // dIt=k21*Ii-ka*It; // tissue insulin, pmol
        // dIt = f(It,Ii)?
        K[0] = k21*Ii-ka*It;
        K[1] = k21*(Ii+ h/2.0)-ka*(It+ h*K[0]/2.0);
        K[2] = k21*(Ii+ h/2.0)-ka*(It+ h*K[1]/2.0);
        K[3] = k21*(Ii+ h)-ka*(It+ h*K[2]);
        dIt   = It + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIt через Дормана-Принса
/*
        K[0] = k21*Ii-ka*It;
        K[1] = k21*(Ii+ h/5.0)-ka*(It+ h*K[0]/5.0);
        K[2] = k21*(Ii+ 3*h/10.0)-ka*(It+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0);
        K[3] = k21*(Ii+ 4*h/5.0)-ka*(It+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0);
        K[4] = k21*(Ii+ 8*h/9.0)-ka*(It+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) );
        K[5] = k21*(Ii+ h)-ka*(It+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) );
        K[6] = k21*(Ii+ h)-ka*(It+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0);

        dIt   = It + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = It + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/


        // Решаем dIi через Рунге-Кутта
        // dIi=-k21*Ii+vbas+vbol;   // interstitial insulin, pmol
        // dIi = f(It,vbol)?
        K[0] = -k21*Ii+vbas+vbol;
        K[1] = -k21*(Ii+ h*K[0]/2.0)+vbas+(vbol+ h/2.0);
        K[2] = -k21*(Ii+ h*K[1]/2.0)+vbas+(vbol+ h/2.0);
        K[3] = -k21*(Ii+ h*K[2])+vbas+(vbol+ h);
        dIi   = Ii + (K[0] + 2.0*K[1] + 2.0*K[2] + K[3])/6.0 *h;

        // Решение dIi через Дормана-Принса
/*
        K[0] = -k21*Ii+vbas+vbol;
        K[1] = -k21*(Ii+ h*K[0]/5.0)+vbas+(vbol+ h/5.0);
        K[2] = -k21*(Ii+ 3*h*K[0]/40.0 + 9*h*K[1]/40.0)+vbas+(vbol+ 3*h/10.0);
        K[3] = -k21*(Ii+ 44*h*K[0]/45.0 + (-56*h*K[1]/15.0) + 32*h*K[2]/9.0)+vbas+(vbol+ 4*h/5.0);
        K[4] = -k21*(Ii+ 19372*h*K[0]/6561.0 + (-25360*h*K[1]/2187.0) + 64448*h*K[2]/6561.0 + (-212*h*K[3]/729.0) )+vbas+(vbol+ 8*h/9.0);
        K[5] = -k21*(Ii+ 9017*h*K[0]/3168.0 + (-355*h*K[1]/33.0) + 46732*h*K[2]/5247.0 + (49*h*K[3]/176.0) + (-5103*h*K[4]/18656.0) )+vbas+(vbol+ h);
        K[5] = -k21*(Ii+ 35*h*K[0]/384.0 + 500*h*K[2]/1113.0 + (125*h*K[3]/192.0) + (-2187*h*K[4]/6784.0) + 11*h*K[5]/84.0)+vbas+(vbol+ h);

        dIi   = Ii + (35*K[0]/384.0 + 500*K[2]/1113.0 + 125*K[3]/192.0 - 2187*K[4]/6784.0 + 11*K[5]/84.0)*h; // во многих источниках домножается на h, но в одном нету такого.
        // fx = Ii + (5179*K[0]/57600.0 + 7571*K[2]/16695.0 + 393*K[3]/640.0 - 92097*K[4]/339200.0 + 187*K[5]/2100.0 + K[6]/40.0)*h; // Альтернативное решение(при вычитании из первого решения даёт оценку ошибки) Мб надо для расчёта шага? Пока не знаю
*/

        /* */

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

        currentTime += h;
    }

    QVector<double> tick;
    QVector<double> CGB;
    QVector<double> Ipg;
    for (DataZ dt : data)       // чтение и выгрузка нуэных результов из структуры
    {
        cout << dt.t << "\t" <<  dt.gp/1.8 << "\t" <<  dt.Ip*20 << "\t" << endl;
        tick.append(dt.t);
        CGB.append(dt.gp/1.8);
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

