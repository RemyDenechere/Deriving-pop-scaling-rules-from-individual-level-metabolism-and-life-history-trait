%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           Figure PhD paper 1                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("C:/Users/rden/OneDrive - Danmarks Tekniske Universitet/Ph.D/Project_1/Modelling/data")
addpath("C:/Users/rden/OneDrive - Danmarks Tekniske Universitet/Ph.D/Project_1/Modelling/matlab_function")
addpath("C:/Users/rden/OneDrive - Danmarks Tekniske Universitet/Ph.D/Project_1/Modelling/Symbols")
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
load("Maturity_data.mat")

%% Figure 1: Theoretical population growth rate
% Assumtions: 
%   - A is constant with adult size M. 
%   - Resource is constant (implies by A = constant). 
%   - 2 trade-offs: M0 is either constant or proportional to M. 
clf
figure1 = gcf;
% Pop groth rate: 
A = 5.35; 
n = 3/4;
Winf = Grid(10^(-3), 10^(7));
% Winf = linspace(10^(-3), 10^(7), 10^(4));

a = 0.42;

subplot1 = subplot(1,2,1)
ratio = 10^3;
epsR = [0.02, 0.05, 1 ]';
r_max = A.*(1 - n) .* Winf.^(n - 1) .*((1 - a).*log(ratio) + log(epsR)); 
loglog(Winf, r_max, 'k', 'LineWidth', 1.25)
xlabel('Adult size, M')
set(gca, 'XTick',[0.01 100 10^(6)])
ylabel('Maximum population growth rate, r_{max}')
title('A')


subplot2 = subplot(1,2,2)
ratio = Winf/0.01;
epsR = [0.0005, 0.001, 0.01, 1]';
r_max = A.*(1 - n) .* Winf.^(n - 1) .*((1 - a).*log(ratio) + log(epsR)); 
loglog(Winf, r_max, 'k', 'LineWidth', 1.25)
hold on 
plot(Winf, Winf.^(-1/4),'k--', 'LineWidth', 1.25)
xlabel('Adult size, M')
set(gca, 'XTick',[0.01 100 10^(6)])
title('B')

xlim(subplot2,[0.5*10^(-2) 10^(6)]);
xlim(subplot1,[10^(-2) 10^(6)]);
ylim(subplot1,[0.0528783494877261 3.88425212252619]);

% Create textbox
annotation(figure1,'textbox',...
    [0.672838113449621 0.406200066137574 0.0882857142857146 0.0682860449735366],...
    'String','\epsilon_R=0.01',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.237238095511539 0.805286596119937 0.0626517561820777 0.0507936507936516],...
    'VerticalAlignment','middle',...
    'String','\epsilon_R=1',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.148577200751209 0.304653439153445 0.0882857142857148 0.0622354497354442],...
    'VerticalAlignment','middle',...
    'String','\epsilon_R=0.02',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.195981714640099 0.613833333333333 0.0882857142857147 0.0740833333333333],...
    'VerticalAlignment','middle',...
    'String','\epsilon_R=0.05',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.70260373844962 0.194534060846569 0.0882857142857146 0.0507936507936517],...
    'VerticalAlignment','middle',...
    'String','\epsilon_R=0.001',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.798521408117632 0.161013558201066 0.0882857142857149 0.0507936507936515],...
    'VerticalAlignment','middle',...
    'String','\epsilon_R=5x10^{-4}',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.558187020101585 0.401791335978843 0.0882857142857146 0.0507936507936517],...
    'VerticalAlignment','middle',...
    'String','\epsilon_R=1',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.586493055555552 0.59707638888889 0.0363802083333331 0.0730833333333338],...
    'String','3',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create ellipse
annotation(figure1,'ellipse',...
    [0.64040972222222 0.859895833333334 0.0287656250000003 0.045861111111111]);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.588595486111108 0.612069444444445 0.0287656250000003 0.0458611111111111]);

% Create textbox
annotation(figure1,'textbox',...
    [0.638858506944441 0.844020833333335 0.036380208333333 0.0730833333333333],...
    'String','2',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create ellipse
annotation(figure1,'ellipse',...
    [0.86365190972222 0.740833333333334 0.0287656250000003 0.045861111111111]);

% Create textbox
annotation(figure1,'textbox',...
    [0.859895833333331 0.72584027777778 0.0363802083333331 0.0730833333333333],...
    'String',{'1'},...
    'FitBoxToText','off',...
    'EdgeColor','none');

x0=0;
y0=0;
width=16;
height=10; 
set(figure1, 'units', 'centimeters', 'position',[x0,y0,width,height])

set(figure1,'Units','centimeters');
screenposition = get(gcf,'Position'); % get the figure size
set(figure1,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); % make the print paper size fits the figure size
print -dpdf -painters Fig_0_Theoretical_rmax

%% Figure 2: Calculation of R* over adult size M, for various growth strategies

epsa = 0.6; 
f0 = 0.6; 
fc = 0.2;
gamma = 1.9753*10^3; 
q = 0.8; 
M = linspace(10^(-3), 10^(7), 10^(4));

figure1 = figure()
% Varying b (grstrategy): 
subplot1 = subplot(1,2,1)
A_0 = [5, 10] ; b = [0, q-n, 0.2]'; 
LineStyle = ['-', ':']; 
for i = 1:2
    h = A_0(i) * M.^(b)/(epsa*(f0 - fc));
    R_star = (fc * h).* M.^(n - q) / (gamma * (1 - fc)); 
    
    plot(M, R_star, 'Color', 'k', 'LineStyle', LineStyle(i), 'LineWidth', 1.25)
    hold on 
end
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
title('A' ,'position',[0.020574769119154,0.03169512472463,0])
ylabel('Minimum sustainable resource R^*')
xlabel('Adult size, M')


% Varying q parameter: 
subplot2 = subplot(1,2,2)
A_0 = 5 ; b = q-n; 
q = [0.5, 0.75, 1]';

h = A_0 * M.^(b)/(epsa*(f0 - fc));
R_star = (fc * h).* M.^(n - q) / (gamma * (1 - fc)); 
loglog(M, R_star, 'Color', 'k', 'LineWidth', 1.25)
xlabel('Adult size, M')
title('B', 'position', [0.081777313638694,0.065712872215395,0])

xlim(subplot1,[0.0049141264204139 40625778.5455678]);
ylim(subplot1,[0.000660894157161765 0.0314005397084759]);

xlim(subplot2,[0.0287860005734426 813274.595671752]);
ylim(subplot2,[0.000215455548253829 0.0656210196600505]);


set(subplot2,'XMinorTick','on','XScale','log','YMinorTick','on','YScale',...
    'log');

% Create textboxs:---------------------------------------------------------
annotation(figure1,'textbox',...
    [0.403425578831311 0.424162257495587 0.0475115766262406 0.0458553791887129],...
    'VerticalAlignment','middle',...
    'String','b=n-q',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.403425578831311 0.161375661375658 0.0475115766262407 0.0458553791887127],...
    'VerticalAlignment','middle',...
    'String','b<n-q',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.143880926130099 0.650793650793649 0.156560088202867 0.0458553791887127],...
    'VerticalAlignment','middle',...
    'String','Increasing A_0',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

annotation(figure1,'textbox',...
    [0.403425578831311 0.75132098765432 0.0475115766262407 0.0458553791887126],...
    'VerticalAlignment','middle',...
    'String','b>n-q',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.751826901874309 0.583774250440917 0.143432194046308 0.0705467372134018],...
    'VerticalAlignment','middle',...
    'String','Increasing q',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create arrows:-----------------------------------------------------------
annotation(figure1,'arrow',[0.305402425578831 0.305402425578831],...
    [0.619811618165782 0.756614087301586]);

annotation(figure1,'arrow',[0.751929437706725 0.751929437706725],...
    [0.310405643738977 0.714285714285714]);

x0=0;
y0=0;
width=16;
height=10; 
set(gcf, 'units', 'centimeters', 'position',[x0,y0,width,height])

set(gcf,'Units','centimeters');
screenposition = get(gcf,'Position'); % get the figure size
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); % make the print paper size fits the figure size
print -dpdf -painters Fig_0b_Theoretical_R_star

%% Figure 3: Growth and egg size strategies - DATA

col = my_color();

% Benthos parameters: ----------------------------------------------------
paramB = parameters_benthos(500);
Linf_B = Growth_data_Bivalve.Linf;
Winf_B = paramB.c*Linf_B.^(3);
K_B = Growth_data_Bivalve.K; % [1/yr]
A_B = 3*K_B.*Linf_B.^(3/4)*paramB.c^(1/4)*paramB.eta^(-1/12);

% fit model:
FitAB = exp(-0.01941)*Winf_B.^(0.2149);
% FitKB = exp(-1.29602)*Winf_B.^(-0.03510);

% Teleost parameters: ----------------------------------------------------
paramF = parameters_fish(1);
Linf_F = Growth_data_Teleost.Linf;
Winf_F = paramF.c*Linf_F.^(3);
K_F = Growth_data_Teleost.K;
A_F = Growth_data_Teleost.A;

% fit model:
FitAF = exp(1.292)*Winf_F.^(0.04331); 
% FitKF = exp(0.08684)*Winf_F.^( -0.2067);

% Elasmobranch parameters: -----------------------------------------------
paramF = parameters_fish(1);
Linf_E = Growth_data_Elasmobranch.Linf;
Winf_E = paramF.c*Linf_E.^(3);
K_E = Growth_data_Elasmobranch.K;
A_E = 3*K_E.*Linf_E.^(3/4)*paramF.c^(1/4)*paramF.eta^(-1/12);

%  fit model:
FitAE = exp(1.706)*Winf_E.^(0.02404); 
% FitKF = exp(0.5012)*Winf_E.^(-0.226);

% Copepods ----------------------------------------------------------------
i = Growth_data_Copepod.Feeding == 'Active feeders' | Growth_data_Copepod.Feeding == 'Mixed feeders';
j = i == 0;
A_C = Growth_data_Copepod.A; 
Winf_C = Growth_data_Copepod.Winf;

% fit model: 

FitAC_Act = Winf_C(~isnan(Winf_C)).^(0)*  6.8346; 
FitAC_Act = FitAC_Act(i); 
FitAC_Pass = Winf_C(~isnan(Winf_C)).^(0)* exp(0.3652);
FitAC_Pass = FitAC_Pass(j);
%
% plot Von Bertalanffy coefficient: --------------------------------------
% Growth parameter A:
figure()

subplot(2,1,1)
loglog(Winf_B, A_B, 'o', 'Color', col.yellight, 'MarkerSize', 5, 'MarkerFaceColor', [0.99 0.94 0.67] ) % Benthos data A %
hold on 
plot(Winf_F, A_F, '.', 'Color', [col.bleulight, 0.9], 'MarkerSize', 7) % Teleost data A
hold on 
plot(Winf_E, A_E, 'v', 'Color', col.redlight, 'MarkerSize', 5) % Elasmobranch data A
hold on 
plot(Winf_C(i), A_C(i), 's', 'Color', [0.49,0.18,0.56], 'MarkerSize', 5, 'MarkerFaceColor', [0.78,0.6, 0.82]) % Copepods data A active feeders 
hold on 
plot(Winf_C(j), A_C(j), 'd', 'Color', [0.49,0.18,0.56], 'MarkerSize', 5, 'MarkerFaceColor', [0.78,0.6, 0.82]) % Copepods data A Passive feeders 
hold on 
plot(Winf_B, FitAB, '-', 'Color', col.yel, 'Linewidth', 1.5) % Benthos fit A
hold on 
plot(Winf_F, FitAF, '-', 'Color', col.bleu, 'Linewidth', 1.5) % teleost fit A
hold on 
plot(Winf_E, FitAE, '-', 'Color', col.red, 'Linewidth', 1.5) % Elasmobranch fit A
hold on 
plot(Winf_C(i), FitAC_Act, '-', 'Color', [0.49,0.18,0.56], 'Linewidth', 1.5) % Copepods fit A
hold on 
plot(Winf_C(j), FitAC_Pass, '-', 'Color', [0.49,0.18,0.56], 'Linewidth', 1.5) % Copepods fit A

ylabel('Growth coefficient, A [g^{1/4}/yr]')
title('A')
legend('Bivalves', 'Teleost', 'Elasmobranch', 'Copepods', 'Location', 'Southeast', 'EdgeColor', 'none')
set(gca, 'FontSize', 10)                                                                                                                                                                                              


% Relation between egg size and Linf 
% Get Linf: 
Linf_B = Egg_size_data_Benthos.L_i;
Linf_E = Egg_size_data_Elasmobranch.Linf;
Linf_T = Egg_size_data_Teleost.Linf;

% Get Winf: [g]
Winf_B = paramB.c*Linf_B.^(3);
Winf_E = paramF.c*Linf_E.^(3); % ! elasmobranch same c as for Telost 
Winf_T = paramF.c*Linf_T.^(3);
Winf_C = Egg_size_data_Copepods.AdultSize*10^(-3); % mg to g
Winf_M = Egg_size_data_Mammal.AdultSize; % already in g. 

% get w0: [g]
w0_B = Egg_size_data_Benthos.Ww_0;
w0_E = Egg_size_data_Elasmobranch.w0;
w0_T = Egg_size_data_Teleost.w0;
w0_C = Egg_size_data_Copepods.ProgenySize*10^(-3); % mg to g
w0_M = Egg_size_data_Mammal.ProgenySize; % already in g

x = log(Winf_M); y = log(Winf_M./w0_M); 

% Fit: 
Fit_B = exp(14.28)*Winf_B.^(1.07); 
Fit_E = exp(8.718)*Winf_E.^(-0.2819); % 362.5060*Winf_E.^(0); % 
Fit_T = exp(7.30)*Winf_T.^(0.9085);
Fit_C = exp(8.871)*Winf_C.^(0.3256); % 176.6654*Winf_C.^(0); 
Fit_M = exp(1.29)*Winf_M.^(0.1294); 

subplot(2,1,2)
% Data: 
loglog(Winf_B', Winf_B'./w0_B', 'o', 'Color', col.yellight, 'MarkerSize', 5, 'MarkerFaceColor', [0.99 0.94 0.67] )
hold on 
plot(Winf_T, Winf_T./w0_T, '.', 'Color', col.bleulight, 'MarkerSize', 7)
hold on 
plot(Winf_E, Winf_E./w0_E, 'v', 'Color', col.redlight, 'MarkerSize', 5)
hold on
plot(Winf_C, Winf_C./w0_C, 's', 'Color', [0.49,0.18,0.56], 'MarkerSize', 5, 'MarkerFaceColor', [0.78,0.6, 0.82])
hold on 
plot(Winf_M, Winf_M./w0_M, 's', 'Color', [0.51,0.66,0.31], 'MarkerSize', 5, 'MarkerFaceColor', [0.51,0.66,0.31]+ 0.3)
% Fit: 
plot(Winf_B, Fit_B, '-', 'Color', col.yel, 'LineWidth', 1.5)
hold on 
plot(Winf_T, Fit_T, '-', 'Color', col.bleu, 'LineWidth', 1.5)
hold on 
plot(Winf_E, Fit_E, '-', 'Color', col.red,  'LineWidth', 1.5)
hold on 
plot(Winf_C, Fit_C, '-', 'Color', [0.49,0.18,0.56],  'LineWidth', 1.5)
hold on 
plot(Winf_M, Fit_M, '-', 'Color', [0.51,0.66,0.31],  'LineWidth', 1.5)

xlabel('Asymptotic weight, M_{\infty} [g]')
ylabel('Adult:Offspring size ratio, M_{\infty}/M_0 [g]')
title('B')
legend('Bivalve', 'Teleost', 'Elasmobranch', 'Copepod','Mammal', 'Location', 'Northwest', 'EdgeColor', 'none')
set(gca, 'FontSize', 10)

% Print pdf 
x0=0;
y0=0;
width=16;
height=17; 
set(gcf, 'units', 'centimeters', 'position',[x0,y0,width,height])

set(gcf,'Units','centimeters');
screenposition = get(gcf,'Position'); % get the figure size
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); % make the print paper size fits the figure size
print -dpdf -painters Fig_1_Data_growth_fecondity

%% Figure 4: 
% -------------------------------------------------------------------------
% Make Epsilon_R (recruitment efficiency to be a free parameter
% -------------------------------------------------------------------------
    % Linear regression with imposed slope: 
    % General model:
    %      f(x) = (-1/4)*x + b
    % Coefficients (with 95% confidence bounds):
    %        b =       1.681  (1.576, 1.786)

% Set the parameter to calculat r_max for Elasmobranch:--------------------
param = parameters_fish(1); % set the parameters or a random Winf. 
A0 = exp(1.70586) ; b = 0; ratio = 6.1119*10^3;
    minWinf = 1.7412*10^(3); maxWinf = 4.8799*10^(5); % g
syms Eps

h0 = A0/(param.epsa*(param.f0 - param.fc));

% growth coefficient: -----------------------------------------------------
    R = 1; % Resource
    A = param.epsa.*h0.*Winf.^(b).*(1./(1+ (param.h0.*Winf.^(b))./(param.gamma.*R)) - param.fc);% individual growth 
    exp_b = A.*(1 - param.n).* ((1 - param.a) .* log(ratio) + log(param.epsEeg*Eps)); % constant term in linear regression
    S = vpa(solve(exp_b == exp(1.681), Eps)); % solving eps_R for a value of constant linear term establish from the linear regression. 
    
%% Figure 4: Population growth rate 4 strategies

color = [0.00,0.45,0.74 ; 0.85,0.33,0.10 ; 0.93,0.69,0.13; 0.49,0.18,0.56 ; 0.84,0.63,0.13]; % set the color 
color_light = [0.41,0.76,0.99 ; 1.00,0.60,0.43 ; 1.00,0.82,0.39; 0.84,0.63,0.13];  
% blue / red / Jaune / violet
param = parameters_fish(1); % set the parameters or a random Winf. 

[Winf, ~] = Grid(10^(-7), 10^8); % g 
f = [0.3, 0.6, 1]'; % resource
clear p 

figure()
for i = 1:5
    % Add r_max data: 
if i == 1 % Teleost
    subplot(2,2,1)
    x = r_max_Teleost_Hutchings2012.Maximum_weight_g;
    y = r_max_Teleost_Hutchings2012.rmax;
    plot(x, y, 'o', 'Color', [0.41,0.76,0.99], 'MarkerSize', 4);
    hold on 
elseif i == 2 % Elasmobranch                                                transformation from Linf to Winf with c = 0.01 (Teleost value)
    subplot(2,2,2)    
    x = data_rmax_Elasmobranch_Zhou2011.Winf;
    y = data_rmax_Elasmobranch_Zhou2011.M2;
    plot(x, y, 'o', 'Color', [1.00,0.60,0.43], 'MarkerSize', 4);
    hold on 
    plot(r_max_Mammalia_Hutchings2012.Maximum_weight_g , r_max_Mammalia_Hutchings2012.rmax, ...
        'o', 'Color', [0.51,0.66,0.31], 'MarkerSize', 4);
    hold on
end

if i == 5 
    subplot(2,2,1)
    [~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f, param, i, 1);
    indx = Winf >= minWinf & Winf <= maxWinf; 
else 
    subplot(2,2,i)
    [~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f, param, i, 1);
    indx = Winf >= minWinf & Winf <= maxWinf; 
    % 90% confidence calculation:  
    [~ ,r_max_inf, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 2);   % calculation for lower confidence  
    [~ ,r_max_sup, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 3);   % calculation for higher confidence

    idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0 ;
    ciplot(r_max_inf(idx), r_max_sup(idx), Winf(idx), color(i,:), 0.25, 'none')
    hold on 
end

pi = plot(Winf(indx), r_max(1:3, indx), '-',...
    Winf, r_max', '--', 'Color', color(i,:))
p(i) = pi(1);
hold on 
    for j = 1:3
        pj = pi(j);
        pj.LineWidth = 0.5 +  j*0.5;
    end

    if i == 4 
    else 
        for j = 4:6
        pj = pi(j);
        pj.LineWidth = 0.5 +  (j-3)*0.5;
        end
    end

plot(Winf, 10^(0)*Winf.^(param.n-1), 'k--', 'LineWidth', 1.5);
set(gca, 'XTick',[10^(-7) 10^(-4) 10^(-1) 10^(2) 10^(5) 10^(8)])
set(gca, 'YTick',[10^(-2) 10^(-1) 10^(0) 10^(1) 10^(2) 10^(3)])
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
xlim([Winf(1) Winf(end)])
ylim([min(min(r_max)) max(max(r_max))])

if i == 1
    ylabel({'Maximum population','growth rate r_{max} [yr^{-1}]'})
elseif i == 2 
elseif i == 3
    ylabel({'Maximum population','growth rate r_{max} [yr^{-1}]'})
    xlabel('Adult size, M [g]')
elseif i == 4
    xlabel('Adult size, M [g]')
end

end

subplot(2,2,1)
xlim([0.0445266595310423 91714913.4792109]);
ylim([10^(-3) 2.08485164319904]);
subplot(2,2,2)
xlim([0.000626577522456087 116599536.669696]);
ylim([0.00644069496127611 20.1775349442806]);
subplot(2,2,3)
xlim([3.8534816261347e-05 96061089.3779234]);
ylim([0.00785807778090396 17.4760846210439]);
subplot(2,2,4)
xlim([9.97339375552866e-08 0.0180406199581103]);
ylim([3.81735231876277 2280.59479422943]);


% Create textbox
annotation(gcf,'textbox',...
    [0.573815104166665 0.920292735042721 0.0374826388888888 0.0525598290598356],...
    'String','B',...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.135047743055555 0.918935897435887 0.0374826388888888 0.0525598290598356],...
    'String',{'A'},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.135047743055555 0.446077991452979 0.0374826388888888 0.052559829059835],...
    'String','C',...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.574917534722221 0.446077991452979 0.0374826388888888 0.0525598290598355],...
    'String','D',...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% pdf properties:----------------------------------------------------------
x0=0;
y0=0;
width=16;
height=13;
set(gcf, 'units', 'centimeters', 'position',[x0,y0,width,height])

set(gcf,'Units','centimeters');
screenposition = get(gcf,'Position'); % get the figure size
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); % make the print paper size fits the figure size
print -dpdf -painters Fig_2_Rmax_data_4cases
%% Figure 4 second version:
clf
color = [0.77,0.45,0.51; 0.00,0.45,0.74 ; 0.85,0.33,0.10 ; 0.93,0.69,0.13; 0.49,0.18,0.56 ; 0.49,0.18,0.56; ]; % set the color 
color_light = [0.77,0.45,0.51; 0.41,0.76,0.99 ; 1.00,0.60,0.43 ; 1.00,0.82,0.39; 1.00,0.82,0.39];  
% blue / red / Jaune / violet
param = parameters_fish(1); % set the parameters or a random Winf. 

[Winf, ~] = Grid(10^(-7), 10^8); % g 
f = [0.3, 0.6, 1]'; % resource
clear p 

subplot(2,2,[1,2])

for i = 1:6
    
[~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f(2), param, i, 1);
indx = Winf >= minWinf & Winf <= maxWinf; 

% 90% confidence calculation:  
[~ ,r_max_inf, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 2); % calculation conf. interval low boundary
[~ ,r_max_sup, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 3); % calculation conf. interval upper boundary

% Limitation for copepods and elasmobranchs
if i == 3 % Elasmobranchs
    idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 1;
elseif i == 5 | i == 6 % Copepods
    idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0 & Winf < 0.1;
else % Others
    idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0;
end 

if i ==6
  i = i;   
end

if i == 6 
    plot(Winf(indx), r_max(indx), '',  Winf(idx), r_max(idx)', '--',...
    'Color', color(i,:), 'Linewidth', 1.5 )
else 
% Confident interval: based on variation of A from data collected:
ciplot(r_max_inf(idx), r_max_sup(idx), Winf(idx), color(i,:), 0.25, 'none')
hold on 
plot(Winf(indx), r_max(indx), '-',  Winf(idx), r_max(idx)', '--',...
    'Color', color(i,:), 'Linewidth', 1.5 )

end 

end

plot(Winf, 10^(1)*Winf.^(param.n-1), 'k--', 'LineWidth', 1.5);
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
xlim([9.59252353156652e-08 101912836.96271]);
ylim([0.00185088935383643 336.189436874531]);
ylabel('Maximum population growth rate r_{max} [yr^{-1}]')
xlabel('Adult size, M [g]')
title('A')


% Plot fish panel: 

subplot(2,2,3)

x = r_max_Teleost_Hutchings2012.Maximum_weight_g;
y = r_max_Teleost_Hutchings2012.rmax;
plot(x, y, 'o', 'Color', [0, 0.26,0.99], 'MarkerSize', 3, 'LineWidth', 1);
hold on

i = 2; 

[~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f(2), param, i, 1);
indx = Winf >= minWinf & Winf <= maxWinf; 
% 90% confidence calculation:  only for feeding 0.6
[~ ,r_max_inf, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 2);   % calculation for lower confidence  
[~ ,r_max_sup, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 3);   % calculation for higher confidence

idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0 ;
ciplot(r_max_inf(idx), r_max_sup(idx), Winf(idx), color(i,:), 0.25, 'none')
hold on 

plot(Winf(indx), r_max(indx), '-', Winf, r_max', '--',...
    'Color', color(i,:), 'Linewidth', 1.5)

hold on 
plot(Winf(idx), 10^(0)*Winf(idx).^(param.n-1), 'k--', 'LineWidth', 1.5);
ylim([0.01 13])
xlim([min(Winf(idx)) Winf(end)])
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
title('B')

% Plot elasmobranch panel
subplot(2,2,4)

x = data_rmax_Elasmobranch_Zhou2011.Winf;
y = data_rmax_Elasmobranch_Zhou2011.M2;
plot(x, y, 'o', 'Color', [1.00,0.3,0], 'MarkerSize', 3, 'LineWidth', 1);
hold on

Marine = r_max_Mammalia_Hutchings2012(r_max_Mammalia_Hutchings2012.Type == 'Marine', :);
Terrest = r_max_Mammalia_Hutchings2012(r_max_Mammalia_Hutchings2012.Type == 'Terrestrial', :);

plot(Marine.Maximum_weight_g , Marine.rmax, ...
    'x', 'Color', [0,0.66,0], 'MarkerSize', 4, 'LineWidth', 1, 'MarkerFaceColor', [0.51,0.66,0.31]+ 0.3);
hold on
plot(Terrest.Maximum_weight_g , Terrest.rmax, ...
    'o', 'Color', [0,0.66, 0], 'MarkerSize', 3, 'LineWidth', 1, 'MarkerFaceColor', [0.51,0.66,0.31]+ 0.3)

i = 3; 
[~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f(2), param, i, 1);
indx = Winf >= minWinf & Winf <= maxWinf; 
% 90% confidence calculation:  only for feeding 0.6
[~ ,r_max_inf, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 2);   % calculation for lower confidence  
[~ ,r_max_sup, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 3);   % calculation for higher confidence

idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0 ;
ciplot(r_max_inf(idx), r_max_sup(idx), Winf(idx), color(i,:), 0.25, 'none')
hold on 

plot(Winf(indx), r_max(indx), '-', Winf, r_max', '--',...
    'Color', color(i,:), 'Linewidth', 1.5)
    
plot(Winf(idx), 10^(0)*Winf(idx).^(param.n-1), 'k--', 'LineWidth', 1.5);
ylim([0.005 0.5*10^2])
xlim([10^(-1) Winf(end)])
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
title('C')

%% plot: -----------------------------------------------------
x0=0;
y0=0;
width=16;
height=18;
set(gcf, 'units', 'centimeters', 'position',[x0,y0,width,height])

set(gcf,'Units','centimeters');
screenposition = get(gcf,'Position'); % get the figure size
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); % make the print paper size fits the figure size
print -dpdf -painters Fig_2_Rmax_3_panel

%% Figure 5: R^* vs Ma 4 strategies
figure

[Winf_val, ~] = Grid(10^(-7), 10^7);
for i = 1:4

[R_star , ~, ~, ~] = Pop_growth_rate(Winf_val, 0, param , i, 1);

loglog(Winf_val, R_star, 'Color', color(i,:), 'Linewidth', 2)
hold on 
end

% plot proprieties: -------------------------------------------------------
xlim([Winf_val(1) Winf_val(end)])
ylabel('R^* [g.L^{-1}]')
xlabel('Adult size, M [g]')
set(gca, 'FontSize', 10)
legend('Teleost', 'Elasmobranch', 'Bivalves', 'Copepods', 'Location', 'Northwest', 'EdgeColor', 'none')

% pdf properties:----------------------------------------------------------
x0=0;
y0=0;
width=16;
height=10;
set(gcf, 'units', 'centimeters', 'position',[x0,y0,width,height])

set(gcf,'Units','centimeters');
screenposition = get(gcf,'Position'); % get the figure size
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); % make the print paper size fits the figure size
print -dpdf -painters Fig_3_Rstar_data_4cases

%% suplementary: plot Winf/Wm ratio

Winf =  datamaturity.Winf; x = log10(Winf);
Wm = datamaturity.Wm; y = log10(Wm);
phylum = datamaturity.Var4;
Col = my_color();
clr = [Col.bleulight; Col.redlight; Col.yellight]; 
sym = '.'; 
Fit1 =  exp(-3.61528)*(Winf.^1.19749); %  0.4095*Winf  -48.34; % 

subplot(2,1,1)
gscatter(Winf, Wm ,phylum, clr, sym, 20)
ylabel('Maturation weight, M_m [g]')
xlabel(' ')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
hold on 
plot(Winf, Fit,  'k', 'LineWidth', 1.5)
leg1 = legend('Teleosts', 'Elasmobranches', 'Bivalves', 'Fit');
leg1.Location = 'NorthWest'; leg1.EdgeColor = 'none';

subplot(2,1,2)
gscatter(Winf, Winf./Wm ,phylum, clr, sym, 20)
xlabel('Asymptotic weight, M_{\infty} [g]')
ylabel('Maturation ratio, M_{\infty}/M_m')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
hold on 
plot([min(Winf) max(Winf)], [3 3],  'k', 'LineWidth', 1.5)
leg2 = legend('Teleosts', 'Elasmobranches', 'Bivalves');
leg2.Location = 'NorthEast'; leg2.EdgeColor = 'none';
ylim([10^(-2) 10^(7)])
xlim([min(Winf) max(Winf)])

% pdf properties:----------------------------------------------------------
x0=0;
y0=0;
width=16;
height=12.5;
set(gcf, 'units', 'centimeters', 'position',[x0,y0,width,height])

set(gcf,'Units','centimeters');
screenposition = get(gcf,'Position'); % get the figure size
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]); % make the print paper size fits the figure size
print -dpdf -painters Fig_4_Supp_Maturation_size



