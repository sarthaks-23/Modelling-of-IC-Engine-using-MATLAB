% Initial conditions
R = 8.314; %Universal Gas Constant in J/mol-K or Pa-m^3/mol-K.
Cv = 1.5*R;
Cp = 2.5*R;
gamma = Cp/Cv;
T1 = 70+273.15; %Converting Temperature into K
P1 = 1e5; % 1 bar = 100 kPa = 10^5 Pa
%Taking number of moles in this problem to be 1
n = 1;
%Uisng ideal gas relation PV = nRT
V1 = n*R*T1/P1;
eta = 0.75; %Irreversible work efficiency

%Process 1->2 is Adiabatic Compression from 70 C to 150 C
T2 = 150+273.15;
P2 = P1*(T1/T2)^(gamma/(1-gamma));
V2 = V1*(P1/P2)^(1/gamma);
W1_2 = (P2*V2 - P1*V1)/(1-gamma);
Q1_2 = 0;
U1_2 = Cv*(T2-T1);
H1_2 = Cp*(T2-T1);

%Process 2->3 is Constant Pressure Process from 150C to 70C
T3 = T1;
P3 = P2;
V3 = n*R*T3/P3;
W2_3 = P3*(V3-V2);
U2_3 = Cv*(T3-T2);
Q2_3 = Cp*(T3-T2);
H2_3 = Q2_3;

%Process 3->1 is Isothermal Expansion to original state
U3_1 = 0;
H3_1 = 0;
W3_1 = n*R*T1*log(V1/V3);
Q3_1 = W3_1;

%Overall cycle 1->2->3->1:
U_total = U1_2 + U2_3 + U3_1;
H_total = H1_2 + H2_3 + H3_1;
Q_total = Q1_2 + Q2_3 + Q3_1;
W_total = W1_2 + W2_3 + W3_1;

%% Plotting PV,PT,TV curves for the process
%Process 1->2: Adiabatic
V_ad = linspace(V1,V2,100);
P_ad = P1*(V1./V_ad).^(gamma);
T_ad = T1*(V1./V_ad).^(gamma-1);
%Process 2->3: Isobaric
T_ib = linspace(T2,T3,100);
V_ib = V2.*(T_ib/T2);
P_ib = linspace(P2,P2,100);
%Process 3->1: Isothermal
V_iso = linspace(V3,V1,100);
P_iso = P3*(V3./V_iso);
T_iso = linspace(T3,T1,100);

%Plotting PV curve for the process
figure;
hold on;
plot(V_ad,P_ad./1e5,'r-','LineWidth',1.5); %V in m^3, P in Pa
plot(V_ib,P_ib./1e5,'b-','LineWidth',1.5); %V in m^3, P in Pa
plot(V_iso,P_iso./1e5,'g-','LineWidth',1.5); %V in m^3, P in Pa
title('PV Diagram');
xlabel('Volume (in m^3)');
ylabel('Pressure in bar');
legend('Adiabatic Compression','Isobaric Cooling','Isothermal Expansion','Location','best');
grid on;

%Plotting PT curve for the process
figure;
hold on;
plot(T_ad,P_ad./1e5,'r-','LineWidth',1.5);
plot(T_ib,P_ib./1e5,'b-','LineWidth',1.5);
plot(T_iso,P_iso./1e5,'g-','LineWidth',1.5);
title('PT diagram');
xlabel('Temperature in K');
ylabel('Pressure in bar');
legend('Adiabatic compression','Isobaric Cooling','Isothermal Expansion','Location','best');
grid on;
hold off;
%Plotting TV curve for the process
figure;
hold on;
plot(V_ad,T_ad,'r-','LineWidth',1.5);
plot(V_ib,T_ib,'b-','LineWidth',1.5);
plot(V_iso,T_iso,'g-','LineWidth',1.5);
title('TV diagram');
xlabel('Volume (in m^3)');
ylabel('Temperature (in K)');
legend('Adiabatic compression','isobaric Cooling','Isothermal Expansion','Location','best');
grid on;

%% Area under PV curve using integral function
f1 = @(V_ad) P1*(V1./V_ad).^(gamma); %adiabatic process
f2 = P2*(V3-V2);
f3 = @(V_iso) P3*(V3./V_iso); %isothermal process

area = integral(f1,V1,V2)+f2+integral(f3,V3,V1); %Area under PV Curve = total work done

%% Irreversible Work Calculations and plotting its curve
%irreversible work
W1_2_irr = eta*W1_2;
W2_3_irr = eta*W2_3;
W3_1_irr = eta*W3_1;

%irreversible heat transfer
Q1_2_irr = U1_2 + W1_2_irr;
Q2_3_irr = U2_3 + W2_3_irr;
Q3_1_irr = U3_1 + W3_1_irr;

%Total work done (irreversible)
W_irr = W1_2_irr+W2_3_irr+W3_1_irr;

%Total heat transferred in the irreversible process
Q_irr = Q1_2_irr+Q2_3_irr+Q3_1_irr;

%For process 1->2
P1_irr = 1e5;
T1_irr = 273.15+70;
V1_irr = n*R*T1_irr/P1_irr;
P2_irr = P1_irr + eta*(P2-P1_irr);
%For process 2->3
T2_irr = T2;
V2_irr = n*R*T2_irr/P2_irr;
%For process 3->1
T3_irr = T3;
P3_irr = P2_irr;
V3_irr = n*R*T3_irr/P3_irr;

%plotting PV curve
V_ad_irr = linspace(V1_irr,V2_irr,100);
T_ad_irr = linspace(T1_irr,T2_irr,100);
P_ad_irr = P1_irr + eta.*(P_ad-P1_irr);
T_ib_irr = linspace(T2_irr,T3_irr,100);
P_ib_irr = linspace(P2_irr,P3_irr,100);
V_ib_irr = n*R*T_ib_irr./P_ib_irr;
P_iso_irr = linspace(P3_irr, P1_irr, 100);
T_iso_irr = linspace(T3_irr,T1_irr,100);
V_iso_irr = n * R * T3 ./ P_iso_irr; % Volume during isothermal expansion

figure;
hold on;
plot(V_ad_irr,P_ad_irr./1e5,'r--','LineWidth',1.5);
plot(V_ib_irr,P_ib_irr./1e5,'b--','LineWidth',1.5);
plot(V_iso_irr,P_iso_irr./1e5,'g--','LineWidth',1.5);
title('PV Diagram for irreversible process');
xlabel('V(in m^3)');
ylabel('P(in Pa)');
legend('Adiabatic Compression','Constant Pressure Compression','Isothermal expansion');
grid on;

figure;
hold on;
plot(T_ad_irr,P_ad_irr./1e5,'r--','LineWidth',1.5);
plot(T_ib_irr,P_ib_irr./1e5,'b--','LineWidth',1.5);
plot(T_iso_irr,P_iso_irr./1e5,'g--','LineWidth',1.5);
title('PT Diagram for irreversible process');
xlabel('T(in K)');
ylabel('P(in Pa)');
legend('Adiabatic Compression','Constant Pressure Compression','Isothermal expansion');
grid on;

figure;
hold on;
plot(V_ad_irr,T_ad_irr,'r--','LineWidth',1.5);
plot(V_ib_irr,T_ib_irr,'b--','LineWidth',1.5);
plot(V_iso_irr,T_iso_irr,'g--','LineWidth',1.5);
title('TV Diagram for irreversible process');
xlabel('V(in m^3)');
ylabel('T (in K)');
legend('Adiabatic Compression','Constant Pressure Compression','Isothermal expansion');
grid on;

%% Difference in work can also be observed in PV diagram

figure;
hold on;
plot(V_ad_irr,P_ad_irr./1e5,'r--','LineWidth',1.5);
plot(V_ib_irr,P_ib_irr./1e5,'b--','LineWidth',1.5);
plot(V_iso_irr,P_iso_irr./1e5,'g--','LineWidth',1.5);
title('PV Diagram for irreversible process');
xlabel('V(in m^3)');
ylabel('P(in Pa)');
plot(V_ad,P_ad./1e5,'r-','LineWidth',1.5); %V in m^3, P in Pa
plot(V_ib,P_ib./1e5,'b-','LineWidth',1.5); %V in m^3, P in Pa
plot(V_iso,P_iso./1e5,'g-','LineWidth',1.5); %V in m^3, P in Pa
title('PV Diagram (Revresible vs Irreversible)');
xlabel('Volume (in m^3)');
ylabel('Pressure in bar');
legend('Adiabatic Compression (irreversible)','Isobaric Cooling (irreversible)','Isothermal Expansion (irreversible)','Adiabatic Compression (reversible)','Isobaric Cooling (reversible)','Isothermal Expansion (reversible)','Location','best');
grid on;


