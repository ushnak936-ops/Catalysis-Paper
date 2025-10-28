clear, clc, close all

%% Constants and Initial Conditions
Cata_Conc = 0.07;  % Catalyst concentration (g) or 0.175 w/v%

CA0 = 0.042;  % Initial fatty acid concentration 
CC0 = 0.000;  % Initial alkane (RP-I)
CD0 = 0.000;  % Initial alcohol
CE0 = 0.000;  % Initial alkane (RP-II)

C0 = [CA0 CC0 CD0 CE0];  % Initial concentration vector

%% Kinetic parameters (as functions of catalyst concentration)
a1 = 2.1200e-09;  b1 = 2.1800e-10;
a2 = 0.5951;      b2 = 0.3102;
a3 = 1.9953;      b3 = 0.2749;
a4 = 0.8786;      b4 = 0.3861;
a_eq = -111.8000; beq = 511.1900;

K1  = (a1 * Cata_Conc) + b1;
K2  = (a2 * Cata_Conc) + b2;
K3  = (a3 * Cata_Conc) + b3;
K4  = (a4 * Cata_Conc) + b4;
Keq = (a_eq * Cata_Conc) + beq;

%% Solve ODE system
tspan = 0:1:1440;  % Time in minutes
[t, C_model] = ode45(@(t,C) C_Pred(t, C, K1, K2, K3, K4, Keq), tspan, C0);

%% Post-Processing
total_feed = CA0;  
Yield_A = C_model(:,1) ./ total_feed * 100;
Yield_C = C_model(:,2) ./ total_feed * 100;
Yield_D = C_model(:,3) ./ total_feed * 100;
Yield_E = C_model(:,4) ./ total_feed * 100;

sum_yield_model = Yield_C + Yield_E;


%% Plotting
sum_yield_model = Yield_C + Yield_E;

plot(t,sum_yield_model, 'm--', 'LineWidth',2)
hold on

xlabel('Time (mins)'); 
ylabel('Alkanes Yield (mol%)');
legend('Catalyst concentration=0.07g');



%% ============================================================
% Function file
% ============================================================
function dCdt = C_Pred(t, C, k1, k2, k3, k4, keq)
    % Species concentrations
    CA = C(1);  
    CC = C(2); 
    CD = C(3); 
    CE = C(4); 

    % Reaction conditions
    C0_TG = 0.042;     % Initial TG concentration (mmol/L)
    v = 40;            % Liquid volume (mL)
    V = 60;            % Gas volume (mL)
    R = 8.314;         % Gas constant (J/mol.K)
    T = 473.15;        % Temperature (K)
    p0_H2 = 3.5e6;       % Initial H2 partial pressure (Pa) [35 bar]
    
    % Calculate partial pressure of H2 at time t
    pH2 = p0_H2 - ((v * (3*C0_TG + CC + 2*CD + 3*CE) * R * T) / V) / 1e3;  % [Pa]
    
    % Calculate CB based on quasi-steady-state assumption
    CB = (k1 * CA * pH2 + (k3/keq) * CD) / (k2 + k3 * pH2);
    
    % Differential equations
    dCA_dt = -k1 * CA * pH2;
    dCC_dt = k2 * CB;
    dCD_dt = (k3 * CB * pH2) - ((k3/keq) * CD) - (k4 * CD);
    dCE_dt = k4 * CD;
    
    % Output derivatives
    dCdt = [dCA_dt; dCC_dt; dCD_dt; dCE_dt];
end
