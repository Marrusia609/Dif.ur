function dZ = odeDM(t,Z) 

global Vbas tm Tm Dbol Ti ti k M

dZ=zeros(14,1);
gp=Z(1); Il=Z(2); Ip=Z(3); fgut=Z(4); fliq=Z(5); fsol=Z(6); gt=Z(7); 
I1=Z(8); Id=Z(9); Xt=Z(10); Ipo=Z(11); Yt=Z(12); Ii=Z(13); It=Z(14);

%% Basal
vbas=Vbas*100; % pmol/min
%% Bolus
bol1=1/Ti(1)*Dbol(1)*(1./(1+exp(-3*(t+10-ti(1))))).*(1./(1+exp(-3*(-10+ti(1)-t+Ti(1)))));
bol2=1/Ti(2)*Dbol(2)*(1./(1+exp(-3*(t+10-ti(2))))).*(1./(1+exp(-3*(-10+ti(2)-t+Ti(2)))));
bol3=1/Ti(3)*Dbol(3)*(1./(1+exp(-3*(t+10-ti(3))))).*(1./(1+exp(-3*(-10+ti(3)-t+Ti(3)))));
vbol=6000*(bol1+bol2+bol3);
%% Meal
Dig=1176*M+1;
vm=Dig/Tm*heaviside(t-tm)*heaviside(-(t-tm-Tm));

%% ODE
EGP=(k(12)-k(13)*gp-k(14)*Id-k(15)*Ipo)*heaviside(k(12)-k(13)*gp-k(14)*Id-k(15)*Ipo);
kgut=k(20)+(k(21)-k(20))/2*(tanh((5/(2*Dig*(1-k(22))))*(fsol+fliq-k(22)*Dig))-tanh((5/(2*Dig*k(23)))*(fsol+fliq-k(23)*Dig))+2);
    
dgp=EGP+k(17)/k(1)*k(18)*fgut-k(27)-k(32)*(gp-k(33))*heaviside(gp-k(33))-k(25)*gp+k(26)*gt; % plasma glucose mg/dl
dgt=-((k(28)+k(29)*Xt)*gt)/(k(30)+gt)+k(25)*gp-k(26)*gt;
dI1=-k(16)*(I1-Ip/k(5));
dId=-k(16)*(Id-I1);
dXt=-k(31)*Xt+k(31)*((Ip/k(5))-k(2))*heaviside((Ip/k(5))-k(2));
dIl=-(k(6)+k(8))*Il+k(7)*Ip;
dIp=-k(7)*Ip+k(6)*Il+k(11)/k(1)*It-k(10)*Ip;
dfgut=-k(18)*fgut+kgut*fliq;
dfliq=-kgut*fliq+k(19)*fsol;
dfsol=-k(19)*fsol+vm;
dIpo=-k(34)*Ipo+(Yt+k(4))*heaviside(dgp)+(Yt+k(4))*(heaviside(-dgp));
dYt=-k(35)*(Yt-k(36)*(gp/k(24)-k(3)))*heaviside(k(36)*(gp/k(24)-k(3))+k(4))+(-k(35)*Yt-k(35)*k(4))*(heaviside(-k(4)-k(36)*(gp/k(24)-k(3))));
dIt=k(9)*Ii-k(11)*It;
dIi=-k(9)*Ii+vbas+vbol;

dZ=[dgp dIl dIp dfgut dfliq dfsol dgt dI1 dId dXt dIpo dYt dIi dIt]';

end