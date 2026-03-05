%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: 
%               1. Produce table 1 in the paper.
%               2. Replicates algorithm to calibrate series. Note: Std
%               Errors might change due to approximation issues in Pattern
%               Search. Which happens approximately at the 5th decimal.
%          
%               
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
% Date:      March 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
warning('off', 'all');
clc


load("../Data/Calibration.mat")

disp("************************************************")
disp("********TABLE 1*************")
disp("************************************************")

names = ["rho_zt_zt", "rho_zt_zn", "rho_zn_zt", "rho_zn_zn", "rho_g", "sigma_zt", "covar_zt_zn", "sigma_zn", "sigma_g"  ]
param = [xopt(1), xopt(2), xopt(3), xopt(4), xopt(5), sqrt(xopt(6)), xopt(7), sqrt(xopt(8)), sqrt(xopt(9))  ]
std_dev = sqrt(abs(diag(inv(hessiancsd(opt_LL,xopt)))))

T = table( names', param', std_dev, 'VariableNames', {'Parameter', 'Values', 'Std. Dev'});
disp(T)


%% Replication

replicate = 0;

if replicate == 1

    % Import Data

    opts = spreadsheetImportOptions("NumVariables", 9);
    % Specify sheet and rangedelta
    opts.Sheet = "Data";
    opts.DataRange = "A30:I145";
    % Specify column names and types
    opts.VariableNames = ["Year", "Agriculture", "Manufacturing", "Services", "Tradable", "SignalYT", "SignalYN", "Total", "Ratio"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];
    filename = "Data.xlsx";
    
    % Import the data
    Data = readtable(filename, opts, "UseExcel", false);
    lnsignalYT = (Data.SignalYT) - mean(Data.SignalYT);
    lnsignalYN = (Data.SignalYN) - mean(Data.SignalYN);
    T = size(lnsignalYN,1);
    
    % Fit a quadratic polynomial (order = 2) to the data
    order = 2; 
    t = 1:T;
    p = polyfit(t, lnsignalYN, order);
    quadratic_trend = polyval(p, t');
    lnsignalYN = lnsignalYN - quadratic_trend; % Remove the quadratic trend from the time series
    
    lntotalY = (Data.Ratio);
    num=10;
    std_yt = num*std(lnsignalYT)^2;
    std_yn = num*std(lnsignalYN)^2;
    std_y = num*std(lntotalY)^2;
    cov_yTN = num*cov(lnsignalYT', lnsignalYN');
    grY = 100*1.01/100;
    
    neg_LL = @(x) LL_klm(x, [lnsignalYT';lnsignalYN']);
    opt_LL = @(x) LL_klm_opt(x, [lnsignalYT';lnsignalYN']);
    
    lb   = [ -0.99; -0.99; -0.99;-0.99;  -0.99;  1e-6; -cov_yTN(1,2); 1e-6; 1e-6]';
    ub  = [ 0.99;  0.99;  0.99;  0.99;  0.99;  std_yt; cov_yTN(1,2); std_yn;  std_y]';
    
    x0 = lb;
    
    %% Pattern Search + Simulated Annealing
    
    rng(10,'twister') % for reproducibility
    
    options = optimoptions('fmincon','Display','iter','MaxIterations', 6e4, 'MaxFunctionEvaluations', 10e10, 'HessianApproximation', 'bfgs', 'PlotFcn',{@optimplotfval,@optimplotx,@optimplotfirstorderopt});
    [sol1, fval0,~,~,~,~,~] = fmincon(neg_LL,x0,[],[],[],[],lb,ub,[],options);
    options = optimoptions('simulannealbnd', 'InitialTemperature', 100, 'TemperatureFcn','temperaturefast', 'MaxFunctionEvaluations', 1e100, 'MaxIterations', 10e10);
    [sol2, fval1] = simulannealbnd(neg_LL, sol1, lb, ub, options);
    options = optimoptions('patternsearch','Display','iter','PlotFcn',{@psplotbestf, @psplotbestx}, 'MaxIterations', 6e4, 'MaxFunctionEvaluations', 10e10);
    [sol3, fval2] = patternsearch(neg_LL, sol2, [],[],[],[],lb,ub, options);
    options = optimoptions('fminunc','Display','iter','MaxIterations', 6e4, 'MaxFunctionEvaluations', 10e10, 'HessianApproximation', 'bfgs', 'PlotFcn',{@optimplotfval,@optimplotx,@optimplotfirstorderopt});
    [xopt, fval3,~,~,~,~] = fminunc(neg_LL, sol3, options);
    
    names = ["rho_zt_zt", "rho_zt_zn", "rho_zn_zt", "rho_zn_zn", "rho_g", "sigma_zt", "covar_zt_zn", "sigma_zn", "sigma_g"  ]
    param = [xopt(1), xopt(2), xopt(3), xopt(4), xopt(5), sqrt(xopt(6)), xopt(7), sqrt(xopt(8)), sqrt(xopt(9))  ]
    std_dev = sqrt(abs(diag(inv(hessiancsd(opt_LL,xopt)))))
    
    T = table( names', param', std_dev, 'VariableNames', {'Parameter', 'Values', 'Std. Dev'});
    disp(T)

end

