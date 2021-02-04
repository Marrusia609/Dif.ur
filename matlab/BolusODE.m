function dZ = odeDM(t,Z)


global Dig ti1 Dbol Vbas mt kp1 g0 gt0 Ip0 b c Fcns m2 Vmxx del tm Ti1 Ti2 ti2 tms Tm Vg

%%
dZ=zeros(14,1);
gp=Z(1); % plasma glucose, mg/kg
Il=Z(2); % liver insulin, pmol/kg 
Ip=Z(3); % plasma insulin, pmol/kg
fgut=Z(4); % gut carbs, mg
fliq=Z(5); % liquid stomach carbs, mg
fsol=Z(6); % solid stomach carbs, mg
gt=Z(7); % tissue glucose, mg/kg
I1=Z(8); % 1st compartment insulin signal, pmol/l
Id=Z(9); % delayed insulin signal, pmol/l
Xt=Z(10); %  interstitial insulin, pmol/l
Ipo=Z(11); % portal vein insulin, pmol/kg
Yt=Z(12); % ??? , pmol/kg/min
Ii=Z(13); % pmol
It=Z(14); % pmol

%% Initial
Gpb=g0; % mg/kg
Gtb=gt0; % mg/kg
Ipb=Ip0; % pmol/kg
Ib=Ipb*20; % pmol/l
EGPb=1.92; % mg/kg/min
Gb=Gpb/1.8; % mg/dl
Sb=1.55; % pmol/kg/min
%% Insulin
V1=0.05; %l/kg
m1=0.190; %min^-1
m6=0.6471; %dimensionless
m3t=(m6*m1)/(1- m6);
%% Insulin Kraegen
k21=2.97e-2; % КС всасывания инсулина в ткани, min-1 Cobelli
di=12e-2; % КС деградации инсулина, min-1
ka=11.3e-3; % КС всасывания инсулина в плазму, min-1
%% Endogenous glucose production
% kp1=2.70; %mg/kg/min
kp2=0.0021; %min^-1
kp3=0.009; %mg/kg/min per pmol/l
kp4=0.0618; %mg/kg/min per pmol/kg
ki=0.0079; %min^-1
%% Glucose appearance
fract=0.9; % Assimilated glucose fraction
kgabs=0.057; % RC glucose absorption, min-1,  Dalla Man
kgri=0.056; % RC grinding, min-1,  Dalla Man
kmin=0.008; % min rate of gut empying
kmax=0.056; % max rate of gut empying
% c=0.115; % percentage, kgut=(kmin+kmax)/2 left обеспечивает 2 волны
% b=0.72; % percentage, kgut=(kmin+kmax)/2 right
Vg=1.8; % plasma per BW, dl/kg
%% Glucose compartments
k1gg=0.065;
k2gg=0.079;
%% Glucose Utilization
% Fcns=0.8; %mg/kg/min
Vm0=2.0; %mg/kg/min
% Vmxx=0.047; %mg/kg/min per pmol/l
Km0=205.59; %mg/kg
p2U=0.0731; %min^-1
%% Renal scretion
ke1=0.0005; % Rate of renal excretion, min-1, DM
ke2=339/Vg; % Threshold of renal excretion, mg/dl
%% Secretion
gamma=0.5;  % min^-1
alpha=0.050; % min^-1
betha=0.11; % pmol/kg per (mg/dl)
%% Fluctuations

vbas=(Vbas+0.00001)*6000/60; % pmol/min
Abol=(Dbol+0.001)*6000;
bol1=1/Ti1*(1-del)*(1./(1+exp(-3*(t-tm+10-ti1)))).*(1./(1+exp(-3*(tm-10+ti1-t+Ti1))));
bol2=1/Ti2*del*(1./(1+exp(-3*(t-tm+10-ti2)))).*(1./(1+exp(-3*(tm-10+ti2-t+Ti2))));
vbol=Abol*(bol1+bol2);
vm=Dig/Tm*heaviside(t-tm-tms)*heaviside(-(t-tm-Tm-tms)); % mg/min

%% 
alp=5/(2*Dig*(1-b));
bet=5/(2*Dig*c);

EGP=(kp1-kp2*gp-kp3*Id-kp4*Ipo)*heaviside(kp1-kp2*gp-kp3*Id-kp4*Ipo);
Vmx=Vm0+Vmxx*Xt;
Uid=(Vmx*gt)/(Km0+gt);
Gtb=(Fcns-EGPb+k1gg*Gpb)/k2gg;
kgut=kmin+(kmax-kmin)/2*(tanh(alp*(fsol+fliq-b*Dig))-tanh(bet*(fsol+fliq-c*Dig))+2);
G=gp/Vg;
fmeal=fract/mt*kgabs*fgut; % Meal
eren=ke1*(gp-ke2)*heaviside(gp-ke2); % Renal glucose excretion
    
dgp=EGP+fmeal-Fcns-eren-k1gg*gp+k2gg*gt; % plasma glucose mg/dl
dgt=-Uid+k1gg*gp-k2gg*gt;
dI1=-ki*(I1-Ip/V1);
dId=-ki*(Id-I1);
dXt=-p2U*Xt+p2U*((Ip/V1)-Ib)*heaviside((Ip/V1)-Ib);
dIl=-(m1+m3t)*Il+m2*Ip; % Liver insulin pmol/kg
dIp=-m2*Ip+m1*Il+ka/mt*It-di*Ip; % Plasma insulin pmol/kg
dfgut=-kgabs*fgut+kgut*fliq; % Gut glucose
dfliq=-kgut*fliq+kgri*fsol; % Glucose in liquid phase
dfsol=-kgri*fsol+vm; % Glucose in solid phase
dIpo=-gamma*Ipo+(Yt+Sb)*heaviside(dgp/Vg) + (Yt+Sb)*(heaviside(-dgp/Vg));
dYt=-alpha*(Yt-betha*(G-Gb))*heaviside(betha*(G-Gb)+Sb)+(-alpha*Yt-alpha*Sb)*(heaviside(-Sb-betha*(G-Gb)));
dIt=k21*Ii-ka*It; % tissue insulin, pmol
dIi=-k21*Ii+vbas+vbol;   % interstitial insulin, pmol

dZ(1)=dgp;  %g0=Gpb/Vg
dZ(2)=dIl;
dZ(3)=dIp;
dZ(4)=dfgut;
dZ(5)=dfliq;
dZ(6)=dfsol;
dZ(7)=dgt; % gt0=Gtb/Vg
dZ(8)=dI1;
dZ(9)=dId;
dZ(10)=dXt;
dZ(11)=dIpo;
dZ(12)=dYt;
dZ(13)=dIi;
dZ(14)=dIt;

end

