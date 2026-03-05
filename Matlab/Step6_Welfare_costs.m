%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: Welfare costs simulations
%                   
%               
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
%
% Date: March 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping

clearvars
clear global
close all

% Compute welfare for recalibrated economy:

recalibrate = false;

% Load parameters

load('Param.mat')
fprintf("Parameters loaded... \n")

% Load databases

load('FICE.mat')
load('FIP.mat')
load('IIPCC.mat')
load('IICECC.mat')

fprintf("Databases loaded... \n")

% Deep parameters

r = Param.r;
b = Param.b;
g =  Param.g;  % Mean growth rate of permanent component
bn = Param.bn;
eta = Param.eta;
rho = Param.rho;
beta = Param.beta;
kappa = Param.kappa;
omega = Param.omega;

% Simulation parameters

sigmag = Param.sigmag; % Std. dev growth rate of permanent component
Tsim = Param.Tsim;   % Simulation points
burn = Param.burn; % Burn-in period for simulation
nstd = Param.nstd;
init_debt = Param.init_debt;
horizn = Param.horizn;

%% Simulate the Value Functions

% Initialize the vectors

FIP.VSim = zeros(1,Tsim - burn);
FICE.VSim = zeros(1,Tsim - burn);
IIPCC.VSim = zeros(1,Tsim - burn);
IICECC.VSim = zeros(1,Tsim - burn);

for i = 2:Tsim-burn
    FIP.VSim(i) = FIP.V(FICE.Index(i), FICE.SimB(i-1));
    FICE.VSim(i) = FICE.V(FICE.Index(i), FICE.SimB(i-1));
    IIPCC.VSim(i) = IIPCC.V(IICECC.Index(i), IICECC.SimB(i-1));
    IICECC.VSim(i) = IICECC.V(IICECC.Index(i), IICECC.SimB(i-1)) ;    

    if mod(i,(Tsim-burn)/10)==0
        disp(i)
    end
end

FICE.Welfcost_fi =     (((FIP.VSim(1, 2:end)) ./ (FICE.VSim(1, 2:end))).^ (1/(1 - rho)) - 1)*100;
FICE.Welfcost_crossiipsfice = (((IIPCC.VSim(1, 2:end)) ./ (FICE.VSim(1, 2:end))).^ (1/(1 - rho)) - 1)*100;
IICECC.Welfcost_ii =  (((IIPCC.VSim(1, 2:end)) ./ (IICECC.VSim(1, 2:end))).^ (1/(1 - rho)) - 1)*100;
IICECC.Welfcost_crossiicevsfip = (((FIP.VSim(1, 2:end)) ./ (IICECC.VSim(1, 2:end))).^ (1/(1 - rho)) - 1)*100;
IIPCC.Welfcost_pp = (((FIP.VSim(1, 2:end)) ./ (IIPCC.VSim(1, 2:end))).^ (1/(1 - rho)) - 1)*100;
IICECC.Welfcost_ce = (((FICE.VSim(1, 2:end)) ./ (IICECC.VSim(1, 2:end))).^ (1/(1 - rho)) - 1)*100;

% Descriptive statistics

FICE.mWelfcost_fi = mean(FICE.Welfcost_fi);
FICE.stdWelfcost_fi = std(FICE.Welfcost_fi);
IICECC.mWelfcost_ii = mean(IICECC.Welfcost_ii)
IICECC.stdWelfost_ii = std(IICECC.Welfcost_ii);

% Remove extreme outliers to compute mean and std

p = prctile(FICE.Welfcost_crossiipsfice, [.25, 99.75]);
FICE.Welfcost_crossiipsfice = FICE.Welfcost_crossiipsfice(FICE.Welfcost_crossiipsfice>=p(1) & FICE.Welfcost_crossiipsfice<=p(2));
p = prctile(IICECC.Welfcost_crossiicevsfip, [.25, 99.75]);
IICECC.Welfcost_crossiicevsfip = IICECC.Welfcost_crossiicevsfip(IICECC.Welfcost_crossiicevsfip>=p(1) & IICECC.Welfcost_crossiicevsfip<=p(2));
p = prctile(IIPCC.Welfcost_pp, [.25, 99.75]);
IIPCC.Welfcost_pp = IIPCC.Welfcost_pp(IIPCC.Welfcost_pp>=p(1) & IIPCC.Welfcost_pp<=p(2));
p = prctile(IICECC.Welfcost_ce, [.25, 99.75]);
IICECC.Welfcost_ce = IICECC.Welfcost_ce(IICECC.Welfcost_ce>=p(1) & IICECC.Welfcost_ce<=p(2));


FICE.mWelfcost_crossiipsfice = mean(FICE.Welfcost_crossiipsfice)
FICE.stdWelfcost_crossiipsfice = std(FICE.Welfcost_crossiipsfice);
IICECC.mWelfcost_crossiicevsfip = mean(IICECC.Welfcost_crossiicevsfip)
IICECC.stdWelfcost_crossiicevsfip = std(IICECC.Welfcost_crossiicevsfip);
IICECC.mWelfcost_ce = mean(IICECC.Welfcost_ce)
IICECC.stdWelfcost_ce = std(IICECC.Welfcost_ce);
IIPCC.mWelfcost_pp = mean(IIPCC.Welfcost_pp)
IIPCC.stdWelfcost_pp = std(IIPCC.Welfcost_pp);

%% Save Files

save('IIPCC.mat','IIPCC', '-v7.3');
save('IICECC.mat','IICECC', '-v7.3');
save('FICE.mat', 'FICE', '-v7.3');

