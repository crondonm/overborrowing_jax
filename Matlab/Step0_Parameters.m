%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: 
%               1. Create Parameter Matrix
%               2. Replicate Table 2 in the Paper.
%               
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
%
% Date: March, 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping

clearvars
clear global
close all


%% Parameters to control Simulation

fprintf("Starting to create parameters ... \n")
%19
%501
Param.nodes = 10;           % Nodes state-space exogenous variables
Param.grid_points = 101;  % Grid for debt
Param.burn = 1000000;    % Burning period       
Param.Tsim = 1000000+Param.burn; % Horizon of simulation

%% Parameters to control IRFs

Param.horizn    = 21;
Param.init_debt = floor(Param.grid_points/2);  % IRFs start must start at Steady State for Debt


%% Parameters to control degree of information

Param.rho_g = 0.496779121757470;    % Autocorrelation permanent growth 
Param.sigmag = 0.00327510105056510; % Std. dev. growth rate of permanent component 
Param.g = log(1.0131);  % Mean growth rate: Garcia-Cicco et al (2010). 

%% Parameters for AR(1) of endowments

% Autocorrelation Matrix

Param.A  = [0.734679129925539,  -0.255320869154147 ,0;
                  0.0336525070039930, 0.417047371095781, 0;
                  0, 0, Param.rho_g];

% Variance-Covariance Matrix

Param.Sigma = [0.00461817366335938, 0.000379607811872954, 0;
                        0.000379607811872954, 0.00137117483970497, 0;
                        0, 0, Param.sigmag];
        
%% Deep Parameters

% Calibrated to match data

Param.beta = 0.83; % Argentina Average NFA-GDP: -29%
Param.kappa = 0.335; % Argentin: Frequency of crises: 5.5%

% Following Bianchi (2011)

Param.r = 0.04; % Annual Interest Rate
Param.rho = 2;    % Risk Aversion Coef
Param.omega = 0.31; % Weight of Tradables 
Param.eta   = 0.83; % Intratemporal Elasticity of Substitution
Param.ytn = Param.nodes;  % Transitory Component of tradable endowment
Param.ynn = Param.nodes; % Transitory Component of Non-tradable endowment
Param.gn = Param.nodes;   % Permanent Component
Param.gTn = Param.nodes; % Growth Rate Tradable Output
Param.gNn = Param.nodes; % Growth Rate Nontradable Output

Param.bmin = -1.4 ;
Param.bmax = -0.2 ;
Param.bn = Param.grid_points;
Param.b = linspace(Param.bmin,Param.bmax,Param.bn) ;

%% ----> Other Parameters

Param.nstd = 1; % Number of standard deviations to consider Financial Crises
Param.window = 5; % Number of years around crises

%% Save

fprintf("Done \n")

save('Data/Param.mat', '-v7.3')

%% Table 2

names = ["Parameter", "Meaning", "Value", "Source/Target" ];
param = ["r", "\sigma", "\varepsilon", "\omega", "\beta", "\kappa", "\mu_g"]';
value  = [Param.r, Param.rho, Param.eta, Param.omega, Param.beta, Param.kappa, Param.g]';
def = ["Interest Rate", "Inverse of Intertemporal elasticity of substitution", "Elasticity between tradable and non-tradable", "Weight between tradable and nontradable goods", "Discount factor", "Borrowing Constraint Constant", "Avg growth rate"]';
    
T = table( param, value, def);

dispAsciiTable(T)

function dispAsciiTable(T)
    % Print header
    fprintf('+------------+------------+-------------------------+\n');
    fprintf('| %-10s | %-10s | %-23s |\n',...
            T.Properties.VariableNames{:});
    fprintf('+------------+------------+-------------------------+\n');

    % Print rows
    for r = 1:height(T)
        fprintf('| %-10s | %-10.3f | %-23s |\n',...
            string(T{r,1}), T{r,2}, string(T{r,3}));
    end
    fprintf('+------------+------------+-------------------------+\n');
end