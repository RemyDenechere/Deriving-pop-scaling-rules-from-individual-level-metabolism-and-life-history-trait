%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           Figure PhD paper 1                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
addpath('data', 'Suppl_code')

% set tha parameters
param = parameters();
col = my_color();

%% Figure 1: Theoretical population growth rate
% Assumtions: 
%   - A is constant with adult size M. 
%   - Resource is constant (implies that A = constant). 
%   - 2 strategies: M0 is either constant or proportional to M. 

figure()
% Constant strategy: ------------------------------------------------------
subplot1 = subplot(1,2,1)
ratio = 10^3; % adult:offsprping size ratio M/M_0
epsR = [0.02, 0.05, 1]'; % recruitment efficiency 
r_max = param.A.*(1 - param.n) .* param.M.^(param.n - 1) .*((1 - param.a).*log(ratio) + log(epsR)); 

loglog(param.M, r_max, 'k', 'LineWidth', 1.25)
set(gca, 'XTick',[0.01 100 10^(6)])
xlabel('Adult size, M')
ylabel('Maximum population growth rate, r_{max}')
title('A')

% Proportional strategy:---------------------------------------------------
ratio = param.M/0.01; % adult:offsprping size ratio M/M_0
epsR = [0.0005, 0.001, 0.01, 1]'; % recruitment efficiency 
r_max = param.A.*(1 - param.n) .* param.M.^(param.n - 1) .*((1 - param.a).*log(ratio) + log(epsR)); 

subplot2 = subplot(1,2,2)
loglog(param.M, r_max, 'k', 'LineWidth', 1.25)
hold on 
plot(param.M, param.M.^(-1/4),'k--', 'LineWidth', 1.25)
xlabel('Adult size, M')
set(gca, 'XTick',[0.01 100 10^(6)])
title('B')

% limits for the axis: 
xlim(subplot2,[0.5*10^(-2) 10^(6)]);
xlim(subplot1,[10^(-2) 10^(6)]);
ylim(subplot1,[0.0528783494877261 3.88425212252619]);

% Textboxes: 
Extra_code_F1(gcf)

% Save the figure as a pdf: 
save_pdf(gcf, 'Fig1_Theoretical_rmax',16, 10)

%% Figure 2: Growth and egg size strategies - DATA
figure()

% Plot Growth parameter A: ------------------------------------------------
subplot(2,1,1)
loglog(param.Winf_B, param.A_B, 'o', 'Color', col.yellight, 'MarkerSize', 5, 'MarkerFaceColor', [0.99 0.94 0.67] ) % Benthos data A %
hold on 
plot(param.Winf_F, param.A_F, '.', 'Color', [0, 0.26,0.99], 'MarkerSize', 7) % Teleost data A
hold on 
plot(param.Winf_E, param.A_E, 'o', 'Color', col.redlight, 'MarkerSize', 5) % Elasmobranch data A
hold on 
plot(param.Winf_C(param.ixAct), param.A_C(param.ixAct), 'o', 'Color', [0.49,0.18,0.56], 'MarkerSize', 5, 'MarkerFaceColor', [0.78,0.6, 0.82]) % Copepods data A active feeders 
hold on 
plot(param.Winf_C(param.ixPas), param.A_C(param.ixPas), 'd', 'Color', [0.49,0.18,0.56], 'MarkerSize', 5, 'MarkerFaceColor', [0.78,0.6, 0.82]) % Copepods data A Passive feeders 
hold on 
plot(param.Winf_B, param.FitAB, '-', 'Color', col.yel, 'Linewidth', 1.5) % Benthos fit A
hold on 
plot(param.Winf_F, param.FitAF, '-', 'Color', col.bleu, 'Linewidth', 1.5) % teleost fit A
hold on 
plot(param.Winf_E, param.FitAE, '-', 'Color', col.red, 'Linewidth', 1.5) % Elasmobranch fit A
hold on 
plot(param.Winf_C(param.ixAct), param.FitAC_Act, '-', 'Color', [0.49,0.18,0.56], 'Linewidth', 1.5) % Copepods fit A
hold on 
plot(param.Winf_C(param.ixPas), param.FitAC_Pass, '-', 'Color', [0.49,0.18,0.56], 'Linewidth', 1.5) % Copepods fit A
ylabel('Growth coefficient, A [g^{1/4}/yr]')
title('A')
legend('Bivalves', 'Teleost', 'Elasmobranch', 'Copepod A.F. ', 'Copepod A.F.', 'Location', 'Southeast', 'EdgeColor', 'none')
set(gca, 'FontSize', 10)                                                                                                                                                                                              

% Relation between egg size and Linf 
% 
subplot(2,1,2)
% Data: 
loglog(param.Winf_B2, param.Winf_B2./param.w0_B, 'o', 'Color', col.yellight, 'MarkerSize', 5, 'MarkerFaceColor', [0.99 0.94 0.67] )
hold on 
plot(param.Winf_T2, param.Winf_T2./param.w0_T, '.', 'Color', [0, 0.26,0.99], 'MarkerSize', 7)
hold on 
plot(param.Winf_E2, param.Winf_E2./param.w0_E, 'o', 'Color', col.redlight, 'MarkerSize', 5)
hold on
plot(param.Winf_C2, param.Winf_C2./param.w0_C, 'o', 'Color', [0.49,0.18,0.56], 'MarkerSize', 5, 'MarkerFaceColor', [0.78,0.6, 0.82])
hold on 
plot(param.Winf_M2, param.Winf_M2./param.w0_M, 'o', 'Color', [0.51,0.66,0.31], 'MarkerSize', 5, 'MarkerFaceColor', [0.51,0.66,0.31]+ 0.3)
% Fit: 
plot(param.Winf_B2, param.Fit_B, '-', 'Color', col.yel, 'LineWidth', 1.5)
hold on 
plot(param.Winf_T2, param.Fit_T, '-', 'Color', col.bleu, 'LineWidth', 1.5)
hold on 
plot(param.Winf_E2, param.Fit_E, '-', 'Color', col.red,  'LineWidth', 1.5)
hold on 
plot(param.Winf_C2, param.Fit_C, '-', 'Color', [0.49,0.18,0.56],  'LineWidth', 1.5)
hold on 
plot(param.Winf_M2, param.Fit_M, '-', 'Color', [0.51,0.66,0.31],  'LineWidth', 1.5)

xlabel('Asymptotic weight, M_{\infty} [g]')
ylabel('Adult:Offspring size ratio, M_{\infty}/M_0 [g]')
title('B')
legend('Bivalve', 'Teleost', 'Elasmobranch', 'Copepod', 'Mammal', 'Location', 'Northwest', 'EdgeColor', 'none')
set(gca, 'FontSize', 10)

% Save the figure as a pdf: 

save_pdf(gcf, 'Fig2_Data_growth_fecondity',16, 17)
saveas(gcf, 'Fig2_Data_growth_fecondity.emf')
    
%% Figure 3: Rmax simulation and data
figure
color = [0.77,0.45,0.51; 0.00,0.45,0.74 ; 0.85,0.33,0.10 ; 0.93,0.69,0.13; 0.49,0.18,0.56 ; 0.49,0.18,0.56; ]; % set the color 
color_light = [0.77,0.45,0.51; 0.41,0.76,0.99 ; 1.00,0.60,0.43 ; 1.00,0.82,0.39; 1.00,0.82,0.39];  
% blue / red / Jaune / violet

[Winf, ~] = Grid(10^(-7), 10^8); % asymptotic size range g 
f = [0.3, 0.6, 1]'; % feeding level

subplot(2,2,[1,2])
for i = 1:6
    
[~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f(2), param, i, 1);
if i == 1 
    indx = [];
else 
    indx = Winf >= minWinf & Winf <= maxWinf; 
end

% 90% confidence calculation:  
[~ ,r_max_inf, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 2); % calculation conf. interval low boundary
[~ ,r_max_sup, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 3); % calculation conf. interval upper boundary

% Limitation for copepods and elasmobranchs
if i == 3 % Elasmobranchs
    idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 1;
elseif i == 5 | i == 6 % Copepods
    idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0 & Winf < 0.1;
elseif i == 1 
    idx = Winf == Winf;
else % Others
    idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0;
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

plot(Winf, 10^(0)*Winf.^(param.n-1), 'k--', 'LineWidth', 1.5);
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
xlim([9.59252353156652e-08 101912836.96271]);
ylim([0.00185088935383643 336.189436874531]);
ylabel('Maximum population growth rate r_{max} [yr^{-1}]')
xlabel('Adult size, M [g]')
title('A')


% Teleost panel: ----------------------------------------------------------
% data from hutchings 2012:
subplot(2,2,3)
plot(param.Winf_Hutch, param.rmax_Hutch, '.', 'Color', [0, 0.26,0.99], 'MarkerSize', 9, 'LineWidth', 1);
hold on

% rmax simulation from our model:
[~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f(2), param, 2, 1);
idx = Winf >= minWinf & Winf <= maxWinf; 
plot(Winf(idx), r_max(idx), '-', Winf, r_max', '--',...
    'Color', color(2,:), 'Linewidth', 1.5)
hold on 

% 90% confidence calculation:  only for feeding 0.6
[~ ,r_max_inf, ~, ~] = Pop_growth_rate(Winf, f(2), param, 2, 2);   % calculation for lower confidence  
[~ ,r_max_sup, ~, ~] = Pop_growth_rate(Winf, f(2), param, 2, 3);   % calculation for higher confidence
idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0 ; % positive values only
% area of confidence:
ciplot(r_max_inf(idx), r_max_sup(idx), Winf(idx), color(2,:), 0.25, 'none')
hold on 

% prediction from Metabolic Theo of Ecology: 
plot(Winf(idx), 10^(0.75)*Winf(idx).^(param.n-1), 'k--', 'LineWidth', 1.5);
ylim([0.01 13])
xlim([min(Winf(idx)) Winf(end)])
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
title('B')

% Elasmobranch panel-------------------------------------------------------
subplot(2,2,4)
% Mammals hutchings 2012:
plot(param.Winf_Marine , param.rmax_Marine, ...
    'x', 'Color', [0.51,0.66,0.31]-0.3, 'MarkerSize', 5)
    % 'x', 'Color', [0,0.66,0], 'MarkerSize', 4, 'LineWidth', 1, 'MarkerFaceColor', [0.51,0.66,0.31]+ 0.3);
hold on 
plot(param.Winf_Ter , param.rmax_Ter, ...
    'o', 'Color', [0.51,0.66,0.31], 'MarkerSize', 4, 'MarkerFaceColor', [0.51,0.66,0.31]+ 0.3)
    %'o', 'Color', [0,0.66, 0], 'MarkerSize', 3, 'LineWidth', 1, 'MarkerFaceColor', [0.51,0.66,0.31]+ 0.3)
hold on
% rmax elasmobranchs Zhou 2011:
plot(param.Winf_zhou, param.rmax_zhou, 'o', 'Color', col.red, 'MarkerSize', 5);

% rmax simulations: 
i = 3; 
[~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f(2), param, i, 1);
idx = Winf >= minWinf & Winf <= maxWinf; 
plot(Winf(idx), r_max(idx), '-', Winf, r_max', '--',...
    'Color', color(i,:), 'Linewidth', 1.5)

% 90% confidence calculation:  only for feeding 0.6
[~ ,r_max_inf, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 2);   % calculation for lower confidence  
[~ ,r_max_sup, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 3);   % calculation for higher confidence
idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0 ;
ciplot(r_max_inf(idx), r_max_sup(idx), Winf(idx), color(i,:), 0.25, 'none')

% prediction from MTE: 
plot(Winf(idx), 10^(0)*Winf(idx).^(param.n-1), 'k--', 'LineWidth', 1.5);
ylim([0.005 0.5*10^2])
xlim([10^(-1) Winf(end)])
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
title('C')
hold off
legend('Mammal mar.', 'Mammal ter.', 'Elasmobranch', 'EdgeColor','none', 'Color','none', 'Location', 'southwest')

%% Save the figure as a pdf: 
save_graph(gcf, 'pdf', 'Fig3_Rmax', 16, 18)

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

%% supplementary: Calculation of R* over adult size M, for various growth strategies

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


