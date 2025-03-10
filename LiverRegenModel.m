function dy=LiverRegenModel(t,y)
% This fuction contains the parameters values of LeAnne paper
% This model is generated to study the sensitive parameters of the Liver
% Regeneration Model

% parameters metabolic load M, and apoptosis parameters kapop, thetaa,
% betaa are different/modified from 2009 paper

% Parameter values 
%% metabolic load
% M=16.8; % first paper ref for rat
M = 5.8; %human donor. optimized from 5.8 - 5.9

%% TNF
ktnf = 1.5;%TNF production
kaptnf = 0.9;%TNF decay+

%% JAK
vjak = 2e4;%JAK production
kmjak = 1e4;
kapjak = 4e-1;%JAK decay

%% STAT3
kprost3 = 2;%[proSTAT3]
vst3 = 750;%STAT3 production
kmst3 = 0.4;
kapst3 = 0.1;%STAT3 decay

%% SOCS3
vsoc = 2.4e4;%SOCS3 production
kmsoc = 7e-4;
kisoc = 15e-3;%dissociation constant for SOCS3 inhibition
kapsoc = 0.4;%SOCS3 decay

%% IE
vie = 2.5e2;%IE production
kmie = 18;
kapie = 5;%IE decay

%% GF
kup = 6e-2; %GF uptake by ECM
kgf = 1.125e-1;%GF production
kapgf = 0.23;%GF decay

%% ECM
kapecm = 33;%ECM decay
kdeg = 7; %ECM degradation by TNF
kecm = kapecm + kdeg; %constant ECM production

%% steady state production rates
prtnf = vjak/(kmjak+1) + kaptnf - ktnf*M;% all recipients original
prst3 = - vst3*(kprost3^2)/(kprost3^2 + kmst3*(1 + 1/kisoc)) + vie/(1+kmie) + vsoc/(1+kmsoc)+ kapst3;
prsoc =  - vsoc/(1+kmsoc) + kapsoc;
prgf = kup + kapgf - kgf*M;
prie = - vie/(1 + kmie) + kapie;
prjak = - vjak/(kmjak+1) + kapjak;
    
%% cell cycle rates
kq = (7e-3); %Q->P %rat
kr = (5.4e-2);%R->Q %rat
kp = (4.4e-3);%P->R %rat
kprol = 2.0e-2;%rat; proliferation rate

%% requiescence (P->Q)
kreq = (1e-1); %rat; requiescence rate for primed cells
thetar = 8; %sigmoid parameters for requiescence
betar = 3;
   
%% apoptosis
% all the apoptosis parameters are modified for human data
% these are different from 2009 paper
kapop = (1e-1)/24;%apoptosis. adjusted for human
% kapop = 1e-1;%rat apoptosis rate
thetaa = 0.15;%sigmoid terms for apoptosis
betaa = 0.075; %M from the old derivs was absorbed into this 

%% ODE function
dy = zeros(10,1);
dy(1) = ktnf*M/((y(8)+y(9)+y(10))) - vjak*y(1)/(kmjak+y(1)) - kaptnf*y(1) + prtnf;%TNF  
dy(2) = vjak*y(1)/(kmjak+y(1)) - kapjak*y(2) + prjak;%JAK
dy(3) = vst3*(kprost3^2)*y(2)/(kprost3^2 + kmst3*(1 + y(4)/kisoc)) - vie*y(3)/(y(3)+kmie) - vsoc*y(3)/(y(3)+kmsoc)- kapst3*y(3) + prst3;%STAT3
dy(4) = vsoc*y(3)/(y(3)+kmsoc) - kapsoc*y(4) + prsoc;%SOCS3 
dy(5) = -kdeg*y(1)*y(5) - kapecm*y(5) + kecm;%ECM
dy(6) = vie*y(3)/(y(3)+kmie) - kapie*y(6) + prie;%IE 
dy(7) = kgf*M/((y(8)+y(9)+y(10))) - kup*y(7)*y(5) - kapgf*y(7) + prgf;%GF
dy(8) = -kq*(y(6)-1)*y(8) + kr*y(5)*y(10) + 0.5*(1 + tanh((thetar-y(7))/betar))*kreq*y(9) - kapop*0.5*(1 + tanh((thetaa - (y(8)+y(9)+y(10)))/betaa))*y(8); %Q  - kapop2*0.5/2*(1 + tanh((y(8)+y(9)+y(10) - 1)/betaa2))*y(8)                
dy(9) = -kp*(y(7)-1)*y(9) + kq*(y(6)-1)*y(8) - 0.5*(1 + tanh((thetar-y(7))/betar))*kreq*y(9) - kapop*0.5*(1 + tanh((thetaa - (y(8)+y(9)+y(10)))/betaa))*y(9);%P
dy(10) = -kr*y(5)*y(10) + kp*(y(7)-1)*y(9) + kprol*y(10) - kapop*0.5*(1 + tanh((thetaa - (y(8)+y(9)+y(10)))/betaa))*y(10);%P

end
