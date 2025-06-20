clc;
clear;

%% Initial Conditions of ideal otto cycle
disp('Provide the following initial conditions to simulate an ideal otto cycle:');
T1 = input('Enter Initial Temperature (in K): '); % Initial Temperature in K
P1 = input('Enter Initial Pressure (in Pa): '); % Initial Pressure in Pa
T3 = input('Enter Maximum Temperature (in K):'); % Maximum Temperature in K
rod = input('Enter connecting rod length (in m):'); % Connecting rod length in m
stroke = input('Enter stroke (in m): '); % Stroke in m
bore = input('Enter bore (in m): '); % Bore (cylinder diameter) in m
cr = input('Enter compression ratio: '); % Compression ratio

% Known constants
R_u = 8.314; % Universal gas constant in J/mol-K

%% Calculating Swept and Clearance Volume for the piston

% Swept and clearance volume calculations
V_swept = (pi/4)*bore^2*stroke; % Swept volume in m^3
V_clearance = V_swept/(cr-1); % Clearance volume in m^3
V1 = V_swept + V_clearance; % Total volume at BDC
V2 = V_clearance; % Volume at TDC
R = R_u/(28.7); % Specific gas constant for air in J/kg-K
m = P1*V1/(R*T1); % Mass of air in kg


%% Defining the variable specific heats for the system

% Coefficients for Cp (J/mol.K): Cp = a + b*T + c*T^2 + d*T^3
a = 28.11; % Coefficient a
b = 0.1967e-2; % Coefficient b
c = 0.4802e-5; % Coefficient c
d = -1.966e-9; % Coefficient d

% Function for Cp as a function of temperature
Cp = @(T) a + b.*T + c.*T.^2 + d.*T.^3; % Specific heat capacity (J/mol.K)
Cv = @(T) Cp(T) - R_u; % Cv from Cp and R
gamma = @(T) Cp(T) ./ Cv(T); % Specific heat ratio


%% Thermodynamics for the given parameters

% Process 1->2 Isentropic Compression
T2 = T1*(V1/V2)^(gamma(T1)-1);
P2 = P1*(V1/V2)*(T2/T1);

% Process 2->3 Constant Volume Heat Addition
P3 = P2*T3/T2;
V3 = V2;
cv_in = @(T) m*Cv(T);
Q_in = integral(cv_in,T2,T3);

% Process 3->4 Isentropic Expansion
V4 = V1;
T4 = T3*(V3/V4)^(gamma(T3)-1);
P4 = (T4/T3)*(P3)*(V3/V4);

% Process 4->1 Constant Volume Heat Rejection
cv_out = @(T) m*Cv(T);
Q_out = integral(cv_out,T1,T4);

% Thermal Efficiency
eta = @(T)1-(1/cr)^(gamma(T)-1);
eta_process = eta((T1+T3)/2);

% Work done in the process: Q-W=d(U) and d(U) = 0 for a complete cycle
W = Q_in - Q_out;

%% Plotting PV plot with intake and exhaust strokes

figure;

% Intake stroke (Isochoric process from TDC to BDC)
V_intake = linspace(V2, V1, 1000); % Volume increases during intake
P_intake = P1*ones(size(V_intake)); % Constant pressure (for intake)

% Exhaust stroke (Isochoric process from BDC to TDC)
V_exhaust = linspace(V1, V2, 1000); % Volume decreases during exhaust
P_exhaust = P1*ones(size(V_exhaust)); % Constant pressure (for exhaust)

% Process 1->2 Isentropic Compression
V1_2 = linspace(V1, V2, 1000);
T1_2 = T1*(V1./V1_2).^(gamma(T1)-1);
P1_2 = (V1./V1_2).*(T1_2./T1)*P1;

plot(V_intake, P_intake, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Intake Stroke');
hold on;

plot(V1_2, P1_2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Isentropic Compression');

% Process 2->3 Constant Volume Heat Addition
T2_3 = linspace(T2, T3, 1000);
P2_3 = P2.*(T2_3./T2);
V2_3 = linspace(V2, V3, 1000);

plot(V2_3, P2_3, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Constant Volume Heat Addition');

% Process 3->4 Isentropic Expansion
V3_4 = linspace(V3, V4, 1000);
T3_4 = T3*(V3./V3_4).^(gamma(T3)-1);
P3_4 = (V3./V3_4).*(T3_4./T3)*P3;

plot(V3_4, P3_4, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Isentropic Expansion');

% Process 4->1 Constant Volume Heat Rejection
T4_1 = linspace(T4, T1, 1000);
P4_1 = P4.*(T4_1./T4);
V4_1 = linspace(V4, V1, 1000);

plot(V_exhaust, P_exhaust, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Exhaust Stroke');
plot(V4_1, P4_1, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Constant Volume Heat Rejection');

xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('Otto Cycle PV Diagram with Intake and Exhaust Strokes');
legend('Location','best');
grid on;

%% Adiabatic Parameters vs without using adiabatic parameters
figure;
hold on;
plot([V1 V2], [P1 P2], 'r--', 'LineWidth', 2,'DisplayName','Compression without using adiabatic parameters');
plot(V1_2, P1_2, 'b', 'LineWidth', 1.5,'DisplayName','Compression using adiabatic parameters');
plot([V2 V1], [P3 P4], 'r--', 'LineWidth', 2,'DisplayName','Expansion without using adiabatic paramters');
plot(V3_4, P3_4, 'color','g', 'LineWidth', 1.5,'DisplayName','Expansion using adiabatic parameters');
xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('Otto Cycle P-V Diagram with Geometric Inputs');
grid on;
legend('Location', 'best');

%% Plotting efficiency with respect to compression ratio for a certain gamma at a certain temperature (Tavg = (T1+T3)/2)
Tavg = (T1+T3)/2;
eta = @(r) 1 - (1./r).^(gamma(Tavg)-1);
eta_T1 = @(r) 1-(1./r).^(gamma(T1)-1);
eta_T2 = @(r) 1-(1./r).^(gamma(T2)-1);
eta_T3 = @(r) 1-(1./r).^(gamma(T3)-1);
eta_T4 = @(r) 1-(1./r).^(gamma(T4)-1); 
T = linspace(T1,T3,5);
figure;
r = linspace(1,12,100);
plot(r,eta_T1(r),'b-','LineWidth',1.5,'DisplayName',['gamma =' num2str(gamma(T1))]);
hold on
plot(r,eta_T2(r),'g-','LineWidth',1.5,'DisplayName',['gamma =' num2str(gamma(T2))]);
plot(r,eta(r),'r-','LineWidth',1.5,'DisplayName',['gamma =' num2str(gamma(Tavg))]);
plot(r,eta_T3(r),'c-','LineWidth',1.5,'DisplayName',['gamma =' num2str(gamma(T3))]);
plot(r,eta_T4(r),'y-','LineWidth',1.5,'DisplayName',['gamma =' num2str(gamma(T4))]);
grid on;
plot(cr,eta(cr),'ko','LineWidth',1.5,'DisplayName','System efficiency we are using right now');
legend ('Location','best');
xlabel('compression ratio (r)');
ylabel('efficiency');
title('efficiency(eta) vs compression ratio(r)');

%% Plotting log curves of PV

figure;hold on;
plot(log10(V_intake),log10(P_intake),'r--','LineWidth',1.5); 
plot(log10(V1_2),log10(P1_2),'r-','LineWidth',1.5); 
plot(log10(V2_3),log10(P2_3),'b-','LineWidth',1.5);
plot(log10(V3_4),log10(P3_4),'g-','LineWidth',1.5);
plot(log10(V4_1),log10(P4_1),'k-','LineWidth',1.5);
plot(log10(V_exhaust),log10(P_exhaust),'b--','LineWidth',1.5); 
title('log PV curve for the cycle');
xlabel('log_1_0V');
ylabel('log_1_0P');
grid on;
legend('intake stroke','isentropic compression','constant volume heat addition','isentropic compression','constant volume heat rejection','exhaust stroke','location','best');

%% Plotting change in pressure with respect to crank angle

% Crank Angle
theta1 = linspace(0, pi, 1000); % Intake
theta2 = linspace(pi, 2*pi, 1000); % Compression
theta3 = linspace(2*pi, 3*pi, 1000); % Heat addition
theta4 = linspace(3*pi, 4*pi, 1000); % Expansion
theta5 = linspace(4*pi, 5*pi, 1000); % Heat Rejection
theta6 = linspace(5*pi,6*pi,1000); %Exhaust
% Volume function (adjust based on engine geometry)
V_crank = @(theta) V_clearance + (V1 - V_clearance) * (1 - cos(theta))/2 + sqrt((rod^2) - (V1^2) * sin(theta).^2);

% Pressure calculations for each phase
P_crank1 = P_intake; % Intake
P_crank2 = P1_2;      % Compression
P_crank3 = P2_3;      % Heat Addition
P_crank4 = P3_4;      % Expansion
P_crank5 = P4_1;      % Heat Rejection
P_crank6 = P_exhaust;  % Exhaust

% Plotting Pressure as a function of crank angle
figure;
hold on;
plot(theta1, P_crank1, 'g-', 'LineWidth', 1.5);
plot(theta2, P_crank2, 'r-', 'LineWidth', 1.5);
plot(theta3, P_crank3, 'b-', 'LineWidth', 1.5);
plot(theta4, P_crank4, 'r-', 'LineWidth', 1.5);
plot(theta5, P_crank5, 'b-', 'LineWidth', 1.5);
plot(theta6, P_crank6, 'g-', 'LineWidth', 1.5);
xlabel('Crank Angle (\theta) [radians]');
ylabel('Pressure (Pa)');
title('Pressure vs Crank Angle for the Otto Cycle');
grid on;
legend('Intake', 'Compression', 'Heat Addition', 'Expansion', 'Heat Rejection', 'Exhaust');

%% Displaying Results
disp('Results for the simulation');
fprintf('Theoretical Efficiency of the system: %f\n',eta_process);
fprintf('State variables [P1, V1, T1]: [%f kPa, %f m^3, %f K]\n',P1/1000,V1,T1);
fprintf('State variables [P2, V2, T2]: [%f kPa, %f m^3, %f K]\n',P2/1000,V2,T2);
fprintf('State variables [P3, V3, T3]: [%f kPa, %f m^3, %f K]\n',P3/1000,V3,T3);
fprintf('State variables [P4, V4, T4]: [%f kPa, %f m^3, %f K]\n',P4/1000,V4,T4);
fprintf('Heat (input) of cycle: %f kJ\n',Q_in/1000);
fprintf('Heat (output) of cycle: %f kJ\n',Q_out/1000);
fprintf('Work done of the cycle: %f kJ\n',W/1000);
fprintf('Mean Effective Pressure (MEP): %f kPa\n',W/(V_swept*1000));
