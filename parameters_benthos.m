function param = parameters_benthos(Winf)
% Parameters is a function producing a structured variable containing all
% the physiological and ecological parameters for Benthos

%% Body size
param.c = 3.2*10^(-2); % g/cm^3, Length to weight coefficient 
param.Linf = (Winf./param.c).^(1/3); % g, Asymptotic length
param.eta = 16*10^(-2); % g, ratio between asymptotic size and size at maturation
param.W0 = 5.3*10^(-7); 0.001; %Fish;   % g, egg size benthos
param.Ws = 2.2*10^(-5); % g, size at settelment
param.Wm = Winf*param.eta; % g, size at maturation
% param.k = 3*0.26; % 1/yr (averaged). 3xK
% param.A = param.k.*B_Winf.^(1/3); % Reproduction coefficient

%% consumption
param.n = 3/4; 
param.epsa = 0.72; % Assimilation efficiency
param.q = 2/3; % Exponent for clearance rate
param.gamma = 3.3540*365; % [g^(-q)L/yr] (4.788019 + 1.9200)/2
param.gamma0 = 0.003962261*365; % [g^(-q)L/yr]
param.qinf = 1/3; % asymptotic size scaling coefficient for clearance rate. 
param.fc = 0.2; % critical feeding level
param.f0 = 0.6; % expected average feeding level for fish
param.h = 7.7745; % 16.3794; % or   g^(1-n)/yr. 

%% Von bertalanffy parameters: 
    % based on the formula : A = A0 Winf^b
param.A0 = exp(-0.04469); % [g^{1/3}/yr] 
param.b = 0.2149; % exponent for growth rate scaling with asymptotic size. 
param.h0 = param.A0/(param.epsa*(param.f0 - param.fc));

%% mortality:
param.a = 0.42; % coefficient of physiological mortality for fish
param.mu0 = param.a*param.epsa*(param.f0 - param.fc)*param.h; % Background mortality coefficient

%% Reproducion
param.qr = 1; % Asymptotic size scalling exponent for maximum reproduction
param.kr = 3*10^(-4)*365; % g^(-qr)/yr, Coefficient for max reproduction rate
param.epsEeg = 0.22; % 0.9;
param.epsR = 0.03;

end

