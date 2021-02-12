function param = parameters_fish(Winf)
% This function contain all the physiological and ecological parameter used
% for Fish (Andersen, book, 2019).

%% Body size
param.c = 0.01; % g/cm^3, Length to weight coefficient 
param.Linf = (Winf/param.c).^(1/3); % g, Asymptotic length
param.eta = 0.28; % g, ratio between asymptotic size and size at maturation
param.W0 = 0.001; % g, egg size
% param.Ws = 2.2*10^(-5); % g, size at settelment
param.Wm = Winf*param.eta; % g, size at maturation

%% consumption:
param.n = 3/4; 
param.epsa = 0.6; % Assimilation efficiency
param.q = 0.8; % Exponent for clearance rate
param.gamma = 1.9753e+03; % [g^(-q)L/yr] 
param.fc = 0.2; % critical feeding level
param.f0 = 0.6; % expected average feeding level for fish

% param.A = 5.35; % VonBertalanffy growth A   g^(1-n)/yr.
% param.h = param.A/(param.epsa*(param.f0-param.fc)); % Maximum consumption parameter g^(1-n)/yr.

%% Von bertalanffy parameters: 
    % based on the formula : A = A0 Winf^b
param.A0 = exp(1.249); % [g^{1/3}/yr] 
param.b = 0.04195; % exponent for growth rate scaling with asymptotic size. 
param.h0 = param.A0/(param.epsa*(param.f0 - param.fc));

%% mortality:
param.a = 0.42; % coefficient of physiological mortality for fish
% param.mu0 = param.a*param.epsa*(param.f0 - param.fc)*param.h; % Background mortality coefficient

%% reproduction: 
param.epsEeg = 0.22;      % reproductive efficiency from Gunderson (1997).
% param.epsR =  0.03;% 0.015; %          %  Recruitment efficiency

end

