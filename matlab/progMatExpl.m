clc 
clear all
close all
%% Global constants  

load('PatConst.mat')

global Vbas tm Tm Dbol Ti ti k bas M ttt

%% Input model parameters
Ginit=122.00;       % +
y0(1)=Ginit*1.8;    %g0  +
y0(2)=2.61;         %Il0 +
y0(3)=5.2040;       %Ip0 +
y0(4)=0;            % fgut mg +
y0(5)=0;            % fliq0 mg +
y0(6)=0;            % fsol0 mg +
y0(7)=169.72;       %gt0 +
y0(8)=104.08;       %I10 +
y0(9)=104.08;       %Id0 +
y0(10)=0;           %Xt0 +
y0(11)=3.08;        %Ipo0 +
y0(12)=0;           %Yt0 +
y0(13)=4120;        %Ii0 +
y0(14)=10830;       %It0 +

%% BasalCalc
bas(1)=1.2238; %Vbas0 +
bas(2)=140; %Gxx ?
bas(3)=110; %Gnn ?
bas(4)=70; %Gt   ? понимаю что это но не понимаю куда ставить и где учитывается)
bas(5)=1.5; %kv ??
bas(6)=bas(1)*(1+bas(5)); %Vmax ??
bas(7)=180; %Gcrit ??

%% Meal & Insulin
M=90;  % + ?
tm=60; % +
Tm=20; % +
ti=[30 60 10];
Ti=[10 10 10];
del = 0.4; %result of bolus calculation
OB = 6.63; %result of bolus calculation
Dbol=[del*OB (1-del)*OB 0];
Vb=bas(1);
Vbas=Vb;
t1=30;
t2=30+120;

%% Time 
progT=720;
time=0:5:progT;

%% ODE solving and Criteria calculating

GG=[];
IP=[];
for i=1:(length(time)-1)

    time1=i*5-5;
    time2=i*5;
    tspan = time1:time2; 
    [t,G] = ode45(@ProgODE, tspan, y0);
    y0=G(end,:);
    GG =[GG ; G(1:end-1,1)/k(24)];
    IP =[IP ; G(1:end-1,3)*20];
    if (time1<t1) || (time1>t2)
        Vbas=basalCalc2(G(end,1)/k(24));
    else
        Vbas=Vb;
    end
end


%% Fitting
figure
hold on
Gpr =[GG(2:end) ; G(end,1)/k(24)];
Ipr=[IP(2:end) ; IP(end,1)*20];
plot(1:720,Gpr,'g');
plot(1:720,IP,'b');
grid on

figure
plot(ttt)
csvwrite('progMatExplDate1.csv',GG);
csvwrite('progMatExplDate2.csv',IP);
