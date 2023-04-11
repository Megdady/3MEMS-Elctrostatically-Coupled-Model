clear all, close all, u = 10^-6;     % u = micro, clc
% "Hamed Design"    : Upper: 28.40, 23.10V, Middle: 31.60, 19.30V, Lower: 27.80, 22.70V
% "Pull-in/Pull-out": Upper: 28.97, 23.05V, Middle: 31.67, 19.13V, Lower: 28.28, 22.58V

Vm1  = 8.5; % 8.5
Vm2  = 27.5; % 27.5
Vm3  = 8.6; % 8.6

Vcd1 = Vm1; InputRatio = 0.814; %OFF 
Vcd3 = Vm3;  % InputVoltage %OFF
Vpp1 = Vm1; %OFF
Vpp2 = 0;
Vpp3 = Vm3; %OFF

% % %%%
% 
% % Vm1  = 0;
% % Vm2  = 0;
% % Vm3  = 0; % GND
% % 
% % Vcd1 = 0; InputRatio = 0.814;
% % Vcd3 = 0;
% % Vpp1 = 0; % InputVoltage
% % Vpp2 = 0;
% % Vpp3 = 0;
% 
% % %%%
% % 
% % Vm1  = 36;
% % Vm2  = 0;
% % Vm3  = 26;
% % 
% % Vcd1 = 24; InputRatio = 0.814; % InputRatio*InputVoltage
% % Vcd3 = 24; % InputVoltage
% % Vpp1 = 12;
% % Vpp2 = 28;
% % Vpp3 = 9;
% 
% % %%%

g = 9.81;
FAmp = -6*2*18*g; AAmp = 0; % The input parameters, Force Amplitude (Shock), Acceleration Amplitude (Base Excitation)

% wn = 3400;
u = 10^-6; % u = micro

rho = 2330; Depth = 50*u; % Density
Mfactor = 2; M1 = 0.9*rho*Depth*58910*u^2*Mfactor; M2 = rho*Depth*75513.66*u^2*Mfactor; M3 = rho*Depth*58910*u^2*Mfactor; MPCB = 0.5*Mfactor; %kg
g0up = (7)*1e-6; g0dn = (14)*1e-6; S = (3.5)*1e-6; % As = b*215*1e-6; Cpar = (1.05 +0.0)*1e-12;
% CapSensMax_I = eps*As*(8/(S) +7/(g0dn+g0up-S)) +Cpar;    CapSensMin_I = eps*As*(8/(g0up) +7/(g0dn)) +Cpar;
Mfactor = 36*0.80; M1 = rho*Depth*58910*u^2*Mfactor; M2 = rho*Depth*75513.66*u^2*Mfactor; M3 = M1; MPCB = 0.5*Mfactor; %kg
Kfactor = 1 ; K1 = 3*Kfactor*.9; K2 = 3*Kfactor*0.9; K3 = 3*Kfactor*.9;
Wn1 = sqrt(K1/M1); Wn2 = sqrt(K2/M2); Wn3 = sqrt(K3/M3); WnPCB = 0.1*Wn2; KPCB = MPCB*WnPCB^2;
zeta = 0.1;
Cfactor = 1; C1 = 2*M1*Wn1*zeta; C2 = 2*M2*Wn2*zeta; C3 = 2*M3*Wn3*zeta; CPCB = 0;%2*MPCB*WnPCB*zeta; %Ns/m

VCD1  = (Vm1 -Vcd1);
VCD3  = (Vm3 -Vcd3);
VPP1  = (Vm1 -Vpp1);
VPP2  = (Vm2 -Vpp2);
VPP3  = (Vm3 -Vpp3);
VPP12 = (Vm1 -Vm2);
VPP32 = (Vm3 -Vm2);

Time = 1 ;% Triangle


%%%
inaccuracy = 0.5;

epsilon = 8.854*u^2; t1 = 50*u; t2 = t1; t3 = t1;
S1 = [-inf 4.4]*u; S2 = [-(4.6) (4.6)]*u; S3 = [-3.7 inf]*u; % Stoppers values, where Si = [xLowerLimit xUpperLimit]
S2 = [-(4.9+inaccuracy) (4.9-inaccuracy)]*u; % S2 upper gap must remain as 4.6*u to match the Pull-in/Pull-out charactatestics

NCD1 = 25  ; wCD1 = 3.5*u; WCD1 = 427.5*u; gCD1 = 5*u ; dCD11 = 6.5*u; dCD12 = 7.5*u; dCD13 = 17*u;
NCD3 = NCD1; wCD3 = wCD1 ; WCD3 = WCD1   ; gCD3 = gCD1; dCD31 = dCD11; dCD32 = dCD12; dCD33 = dCD13;

N1 = 8  ; gppU1 = (7+.4-1.1+.2)*u; gppL1 = (14+.4)*u; L1 = (215-.4)*u; % gppU1 dominant
N2 = 2*5; gppU2 = (7+.4-inaccuracy)*u; gppL2 = (7+.4+inaccuracy)*u; L2 = (215-.4)*u;
N3 = N1 ; gppU3 = (14+.4)*u; gppL3 = (7+.4-1.1+.1)*u; L3 = (215-.4)*u;

N12 = 6; gppU12 = (11.5+.4-0.7-inaccuracy)*u; gppL12 = (20.5+.4+0.7+inaccuracy)*u; L12 = (310-.4)*u;
N32 = 6; gppU32 = (11.5+.4-0.7-inaccuracy)*u; gppL32 = (20.5+.4+0.7+inaccuracy)*u; L32 = (310-.4)*u;




