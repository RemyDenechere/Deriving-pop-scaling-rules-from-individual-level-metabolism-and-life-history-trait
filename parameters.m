function param = parameters()

%% load the data: 
load('Egg_size_data_Benthos.mat')
load('Egg_size_data_Elasmobranch.mat')
load('Egg_size_data_Teleost.mat')
load('Egg_size_data_Copepods.mat')
load('Egg_size_data_Mammal.mat')
load('Growth_data_Copepod.mat')
load('Growth_data_Bivalve.mat')
load('Growth_data_Elasmobranch.mat') 
load('Growth_data_Teleost.mat')
load('r_max_Teleost_Hutchings2012.mat')
load('r_max_Elasmobranch_Zhou2011.mat')
load('r_max_Mammalia_Hutchings2012.mat')
load('Maturity_data.mat')

%% Basic parameter for the plots:

param.A = 5.35; % Fish individual growth rate g^(1-n) year^{-1)
param.n = 3/4; % metabolic exponent
param.M = Grid(10^(-3), 10^(7)); % Asymptotic size range g
param.a = 0.42; % physiological mortality year^(-1)

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 Individual level strategies from data                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.B_c = 3.2*10^(-2); % Benthos length to weight coefficient g/cm^3 
param.B_eta = 16*10^(-2); % g, ratio between asymptotic size and size at maturation
param.F_c = 0.01; % Fish Length 2 weight g/cm^3
param.F_eta = 0.28; % g, ratio between asymptotic size and size at maturation

% Bivalves-----------------------------------------------------------------
Linf_B = Growth_data_Bivalve.Linf; % data
param.Winf_B = param.B_c*Linf_B.^(3); % conversion length to wet weigth
K_B = Growth_data_Bivalve.K; % von Bertalanffy growth coef K 
% Conversion from K to A 
param.A_B = 3*K_B.*Linf_B.^(3/4)*param.B_c^(1/4)*param.B_eta^(-1/12); %g^(-1/4) year^{-1)
% fit:
param.FitAB = exp(-0.01941)*param.Winf_B.^(0.2149);

% Teleost parameters: ----------------------------------------------------
Linf_F = Growth_data_Teleost.Linf;
param.Winf_F = param.F_c*Linf_F.^(3);
param.A_F = Growth_data_Teleost.A;
% fit:
param.FitAF = exp(1.292)*param.Winf_F.^(0.04331); 

% Elasmobranch parameters: -----------------------------------------------
Linf_E = Growth_data_Elasmobranch.Linf;
param.Winf_E = param.F_c*Linf_E.^(3);
K_E = Growth_data_Elasmobranch.K;
% Conversion from K to A 
param.A_E = 3*K_E.*Linf_E.^(3/4)*param.F_c^(1/4)*param.F_eta^(-1/12);
% fit:
param.FitAE = exp(1.706)*param.Winf_E.^(0.02404);

% Copepods ----------------------------------------------------------------
param.ixAct = Growth_data_Copepod.Feeding == 'Active feeders' | Growth_data_Copepod.Feeding == 'Mixed feeders'; % selecte active and mixed feeders. 
param.ixPas = param.ixAct == 0;

param.A_C = Growth_data_Copepod.A; 
param.Winf_C = Growth_data_Copepod.Winf;
% fit model: 
param.FitAC_Act = param.Winf_C(~isnan(param.Winf_C(param.ixAct))).^(0)*6.8346; 
param.FitAC_Pass = param.Winf_C(~isnan(param.Winf_C(param.ixPas))).^(0)* exp(0.3652);

% Get Linf: 
Linf_B = Egg_size_data_Benthos.L_i;
Linf_E = Egg_size_data_Elasmobranch.Linf;
Linf_T = Egg_size_data_Teleost.Linf;

% Get Winf: [g]
param.Winf_B2 = param.B_c*Linf_B.^(3);
param.Winf_E2 = param.F_c*Linf_E.^(3); % ! elasmobranch same c as for Telost 
param.Winf_T2 = param.F_c*Linf_T.^(3);
param.Winf_C2 = Egg_size_data_Copepods.AdultSize*10^(-3); % mg to g
param.Winf_M2 = Egg_size_data_Mammal.AdultSize; % already in g. 

% get w0: [g]
param.w0_B = Egg_size_data_Benthos.Ww_0;
param.w0_E = Egg_size_data_Elasmobranch.w0;
param.w0_T = Egg_size_data_Teleost.w0;
param.w0_C = Egg_size_data_Copepods.ProgenySize*10^(-3); % mg to g
param.w0_M = Egg_size_data_Mammal.ProgenySize; % already in g

% Fit: 
param.Fit_B = exp(14.28)*Winf_B2.^(1.07); 
param.Fit_E = exp(8.718)*Winf_E2.^(-0.2819); % 362.5060*Winf_E.^(0); % 
param.Fit_T = exp(7.30)*Winf_T2.^(0.9085);
param.Fit_C = exp(8.871)*Winf_C2.^(0.3256); % 176.6654*Winf_C.^(0); 
param.Fit_M = exp(1.29)*Winf_M2.^(0.1294); 



