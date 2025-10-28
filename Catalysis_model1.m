%% Catalysis_model.m
% Main script to model HDO kinetics and optimize rate constants

clear; clc; close all;
global k1 k2 k3 k4 keq

% === Initial concentrations ===
CA = 0.042;      % Initial fatty acid concentration (mmol/L)
CC = 0.000;      % Initial alkane (RP-I)
CD = 0.000;      % Initial alcohol
CE = 0.000;      % Initial alkane (RP-II)
C0 = [CA CC CD CE];  

% === Time span ===
tspan = [0:1:1440];  %[min]

% === Example initial kinetic parameters (optional) ===
k1 = 3e-10;
k2 = 0.02;
k3 = 0.4;
k4 = 0.5;
keq = 500;

% === Solve ODEs for base case ===
[t,C_model] = ode45(@kinetic_model_B, tspan, C0);

% === Plot results ===
% defining the Yield (mole percentage)
total_feed = CA(1);     % [mmoles/litre]
Yield_A = C_model(:,1) ./ total_feed * 100;     % [%]
Yield_C = C_model(:,2) ./ total_feed * 100;
Yield_D  = C_model(:,3) ./ total_feed * 100;
Yield_E  = C_model(:,4) ./ total_feed * 100;

 %At 35 bar
t_exp = [0, 120, 360, 720, 1080, 1440];     % [min]

CA_exp = [0.042, 0.022800751, 0.015548694, 0.013424514, 0.007499138, 0.004944984];    % Fatty Acid;
CC_exp = [0.00, 0.001453939, 0.002273772, 0.003194277, 0.003941226, 0.004213656];         % Alkane (RP-I)
CD_exp = [0.000, 0.005323241, 0.004512345,0.001450064, 0.001226393, 0.000912183];    % Alcohol
CE_exp = [0.000, 0.00560754, 0.012775642, 0.018721192, 0.024009201, 0.026868251];      % Alkane (RP-II)

   
  total_feed = CA(1);
% at temperature of 473.15K
%total_C_EXP = CA +CC +CD+ CE;
Yield_A_exp = CA_exp ./total_feed * 100;
Yield_C_exp = CC_exp./total_feed * 100;
Yield_D_exp = CD_exp./total_feed * 100;
Yield_E_exp = CE_exp./total_feed * 100;

% figure
sum_yield_model = Yield_C + Yield_E;
sum_yield_exp = Yield_C_exp + Yield_E_exp;

plot(t,sum_yield_model, 'm--', 'LineWidth',2)
hold on
plot(t_exp,sum_yield_exp, '*')
hold on

xlabel('Time (mins)'); 
ylabel('Alkanes Yield (mol%)');
legend('P=35 bar');


% === Genetic Algorithm optimization ===
 lb = [3e-10, 2.3e-1, 3.7e-1, 4.3e-1, 500];
 ub = [4.5e-10, 3.4e-1, 6e-1, 7.6e-1, 550]; 
options = optimoptions('ga', 'Display', 'iter');

x_best = ga(@kinetic_optimizer_B, 5, [], [], [], [], lb, ub, [], options);
disp(x_best)

%% === Local function definitions (must come after main script) ===

% ---------------------------------------------------------------
function dCdt = kinetic_model_B(t, C)
    % Kinetic Model for HDO of Triglycerides 
    global k1 k2 k3 k4 keq;
    close all

    % State variables
    CA = C(1);  
    CC = C(2); 
    CD = C(3); 
    CE = C(4); 

    % Reaction conditions
    C0_TG = 0.042;   
    v = 40;           
    V = 60;           
    R = 8.314;        
    T = 473.15;       
    p0_H2 = 3.5e6;      
    
    % Calculate partial pressure of H2
    pH2 = p0_H2 - ((v*(3*C0_TG + CC + 2*CD + 3*CE)*R*T)/V)/1e3;
    CB = (k1*CA*pH2 + (k3/keq)*CD) / (k2 + k3*pH2);
    
    % Rate equations
    dCA_dt = -k1*CA*pH2;
    dCC_dt = k2*CB;
    dCD_dt = (k3*CB*pH2) - ((k3/keq)*CD) - (k4*CD);
    dCE_dt = k4*CD;
    
    dCdt = [dCA_dt; dCC_dt; dCD_dt; dCE_dt];
end
% ---------------------------------------------------------------

% Optimizer function for kinetic parameter fitting
function MSE_tot = kinetic_optimizer_B(param)
    global k1 k2 k3 k4 keq
        close all

    k1 = param(1);              
    k2 = param(2);     
    k3 = param(3);                 
    k4 = param(4);
    keq = param(5);


    % Initial concentrations 
    C0 = [0.042, 0, 0, 0];  

    % Normalization values
max_CA = 0.042;
max_CC = 0.004213656;
max_CD = 0.005323241;
max_CE = 0.026868251;


    % Time span
    tspan = 0:1:1440;

    % Solve ODEs
    [t,C_model] = ode45(@kinetic_model_B, tspan, C0);

    % Experimental data
t_exp = [0, 120, 360, 720, 1080, 1440];     % [min]

CA_exp = [0.042, 0.022800751, 0.015548694, 0.013424514, 0.007499138, 0.004944984];    % Fatty Acid;
CC_exp= [0.00, 0.001453939, 0.002273772, 0.003194277, 0.003941226, 0.004213656];         % Alkane (RP-I)
CD_exp = [0.000, 0.005323241, 0.004512345,0.001450064, 0.001226393, 0.000912183];    % Alcohol
CE_exp = [0.000, 0.00560754, 0.012775642, 0.018721192, 0.024009201, 0.026868251];      % Alkane (RP-II)

%measure the error (MSE_tot)
MSE_CA =((C_model(121,1)-CA_exp(1,2))^2)+((C_model(361,1)-CA_exp(1,3))^2)+((C_model(721,1)-CA_exp(1,4))^2)+((C_model(1081,1)-CA_exp(1,5))^2)+((C_model(1441,1)-CA_exp(1,6))^2);
MSE_CC= ((C_model(121,2)-CC_exp(1,2))^2)+((C_model(361,2)-CC_exp(1,3))^2)+((C_model(721,2)-CC_exp(1,4))^2)+((C_model(1081,2)-CC_exp(1,5))^2)+((C_model(1441,2)-CC_exp(1,6))^2);
MSE_CD = ((C_model(121,3)-CD_exp(1,2))^2)+((C_model(361,3)-CD_exp(1,3))^2)+((C_model(721,3)-CD_exp(1,4))^2)+((C_model(1081,3)-CD_exp(1,5))^2)+((C_model(1441,3)-CD_exp(1,6))^2);
MSE_CE = ((C_model(121,4)-CE_exp(1,2))^2)+((C_model(361,4)-CE_exp(1,3))^2)+((C_model(721,4)-CE_exp(1,4))^2)+((C_model(1081,4)-CE_exp(1,5))^2)+((C_model(1441,4)-CE_exp(1,6))^2);


MSE_CA_Nor = MSE_CA/max_CA;
MSE_CC_Nor = MSE_CC/max_CC;
MSE_CD_Nor = MSE_CD/max_CD;
MSE_CE_Nor = MSE_CE/max_CE;

MSE_tot = MSE_CA_Nor + MSE_CC_Nor +MSE_CE_Nor + MSE_CD_Nor;


end
