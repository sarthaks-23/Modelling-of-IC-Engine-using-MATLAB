% Initial Conditions
P1 = 100;
T1 = 273.15 + 30;
V1 = 0.0038;
T3 = 273.15+1200; %peak temperature
r = 8; %compression ratio
gamma = 1.4;
Cv = 0.718;
R = 8.314/28.7; %R in kPa-m^3/kg-K
m = (P1*V1)/(R*T1);

% 1 -> 2 is isentropic compression
V2 = V1/r;
T2 = T1*r^(gamma-1);
P2 = P1*r^(gamma);

% 2 -> 3 constant volume heat addition
qin = m*Cv*(T3-T2);
P3 = P2*T3/T2;
V3 = V2;

% 3 -> 4 isentropic expansion
V4 = V1;
P4 = P3/(r^(gamma));
T4 = T3/(r^(gamma-1));

% 4 -> 1 constant volume heat rejection
qout = m*Cv*(T4-T1);

% First law of thermodynamics: Q-W = del(U). For a complete cycle, del(U) = 0.
Wnet = qin - qout;

% Thermal efficiency:
eta = 1-(1/r^(gamma-1));

% Mean Effective pressure:
MEP = Wnet/(V1-V2);

%Plotting P-V Diagram
%Process 1->2 Isentropic Compression
V1_2 = linspace(V1, V2, 100);
T1_2 = linspace(T1,T2,100);
P1_2 = P1*(V1./V1_2).^(gamma);

%Process 2->3 Constant volume heat addition
V2_3 = linspace(V2,V3,100);
T2_3 = linspace(T2,T3,100);
P2_3 = P2.*T2_3/T2;

%Process 3->4 Isentropic Expansion
V3_4 = linspace(V3,V4,100);
P3_4 = P3*(V3./V3_4).^(gamma);
T3_4 = linspace(T3,T4,100);

%Process 4->1 Constant Volume heat rejection
V4_1 = linspace(V4,V1,100);
T4_1 = linspace(T4,T1,100);
P4_1 = (P4/T4).*T4_1;


figure;
hold on;
grid on;
plot(V1_2,P1_2,'r-','LineWidth',2,'DisplayName','Isentropic Compression');
plot(V2_3,P2_3,'b-','LineWidth',2,'DisplayName','Constant Volume Heat Addition');
plot(V3_4,P3_4,'g-','LineWidth',2,'DisplayName','Isentropic Expansion');
plot(V4_1,P4_1,'k-','LineWidth',2,'DisplayName','Constant Volume Heat Rejection');
legend;
xlabel('Volume (in m^3)');
ylabel('Pressure (in kPa)');

%Calculating area under PV curve
% Process 1 -> 2
f1 = @(V) P1*(V1./V).^(gamma);
val1 = abs(integral(f1,V1,V2));

f2 = @(V) P3*(V3./V).^(gamma);
val2 = abs(integral(f2,V3,V4));

W_curve = val2 - val1;

fprintf('Heat Rejection (Q_out): %.3f kJ\n',qout);
fprintf('Net Work Output (W_net): %.3f kJ\n',Wnet);
fprintf('Thermal Efficiency(η): %.3f %% \n',eta*100);
fprintf('Mean Effective Pressure (MEP): %.3fkPa\n',MEP);





