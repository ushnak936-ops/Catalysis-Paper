%% T_Prediction_Model.m
% Prediction of catalytic reaction at given temperature using Arrhenius parameters

clear; clc; close all;

%% Reaction temperature (K)
T = 473.15;

%% Initial concentrations
CA0 = 0.042;       % Initial fatty acid concentration [mmol/L]
CC0 = 0.000;       % Initial alkane (RP-I) concentration [mmol/L]
CD0 = 0.000;       % Initial alcohol concentration [mmol/L]
CE0 = 0.000;       % Initial alkane (RP-II) concentration [mmol/L]
C0 = [CA0 CC0 CD0 CE0];   % Initial concentration vector

%% Arrhenius parameters
A1  = 2.95465E-08;
A2  = 16.49240546;
A3  = 19.59550102;
A4  = 119.0329089;
Aeq = 107370926.5;

p1  = -2.2226e+03;
p2  = -1.9000e+03;
p3  = -1.9407e+03;
p4  = -2.7725e+03;
peq = -5.7496e+03;

% Calculate rate constants (at temperature T)
k1  = A1  * exp(p1/T);
k2  = A2  * exp(p2/T);
k3  = A3  * exp(p3/T);
k4  = A4  * exp(p4/T);
keq = Aeq * exp(peq/T);

%% Solve the ODE system
tspan = 0:1:1440;  % time in minutes
[t, C_model] = ode45(@(t, C) kinetic_model_prediction(t, C, k1, k2, k3, k4, keq), tspan, C0);

%% Post-processing — calculate yield (mol%)
total_feed = CA0;     % [mmoles/litre]
Yield_A = C_model(:,1) ./ total_feed * 100;     % [%]
Yield_C = C_model(:,2) ./ total_feed * 100;
Yield_D = C_model(:,3) ./ total_feed * 100;
Yield_E = C_model(:,4) ./ total_feed * 100;

%% Experimental data for validation (at T = 473.15 K)
t_exp = [0, 120, 360, 720, 1080, 1440];   % [min]

CA_exp = [0.042, 0.025497754, 0.017307278, 0.016656885, 0.008691079, 0.006988958];    % Fatty Acid
CC_exp = [0.000, 0.001368438, 0.002169306, 0.002643177, 0.003761127, 0.003962345];    % Alkane (RP-I)
CD_exp = [0.000, 0.005030847, 0.004406842, 0.002434665, 0.001237993, 0.000955983];    % Alcohol
CE_exp = [0.000, 0.005323307, 0.011587307, 0.014774524, 0.021575801, 0.023234188];    % Alkane (RP-II)

Yield_A_exp = CA_exp ./ total_feed * 100;
Yield_C_exp = CC_exp ./ total_feed * 100;
Yield_D_exp = CD_exp ./ total_feed * 100;
Yield_E_exp = CE_exp ./ total_feed * 100;

%% Plot comparison of total alkane yield (model vs experimental)
sum_yield_model = Yield_C + Yield_E;
sum_yield_exp = Yield_C_exp + Yield_E_exp;

figure;
plot(t, sum_yield_model, 'm--', 'LineWidth', 2);
hold on;
plot(t_exp, sum_yield_exp, 'k*', 'MarkerSize', 8);
xlabel('Time (min)');
ylabel('Alkanes Yield (mol%)');
legend('Model Prediction (473.15 K)', 'Experimental Data', 'Location', 'best');
grid on;

%% Local function — reaction kinetics
function dCdt = kinetic_model_prediction(t, C, k1, k2, k3, k4, keq)
    % Concentrations
    CA = C(1);
    CC = C(2);
    CD = C(3);
    CE = C(4);

    % Reaction conditions
    C0_TG = 0.042;      % Initial TG concentration (mmol/L)
    v = 40;             % Liquid volume (mL)
    V = 60;             % Gas volume (mL)
    R = 8.314;          % Gas constant (J/mol·K)
    T = 473.15;         % Temperature (K)
    p0_H2 = 3e6;        % Initial H2 partial pressure (Pa) (30 bar)

    % Calculate partial pressure of H2 at time t
    pH2 = p0_H2 - ((v*(3*C0_TG + CC + 2*CD + 3*CE)*R*T)/V)/1e3;  % [Pa]

    % Intermediate species (CB)
    CB = (k1*CA*pH2 + (k3/keq)*CD) / (k2 + k3*pH2);

    % Rate equations
    dCA_dt = -k1*CA*pH2;
    dCC_dt = k2*CB;
    dCD_dt = (k3*CB*pH2) - ((k3/keq)*CD) - (k4*CD);
    dCE_dt = k4*CD;

    dCdt = [dCA_dt; dCC_dt; dCD_dt; dCE_dt];
end
