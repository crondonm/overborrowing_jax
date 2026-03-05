%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code:  Create transition matrix and simulation vectors
%                     Load Parameters.
%                     Discretize Stochastic Processes
%                     Save Transition Matrix and Simulations.
% Ancillary Functions:
%                     tmp.m: Discretize VAR(1) - Schmitt-Ghohé and Uribe (2010)
%                     KF_transition_matrix.m: Compute transition matrix based on simulation and
%                     Kalman Filter.
%                     findClosest2.m: Algorithm to find the closest neighbor.
%
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
%
% Date: March, 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping

clearvars
clear global
close all

maxNumCompThreads(16);

%% Load Parameters

load('Data/Param.mat')

A = Param.A;   % Autocorrelation matrix
Sigma  = chol(Param.Sigma); % Variance - Covariance matrix
ytn = Param.ytn;       % Grid points for tradable component
ynn = Param.ynn;      % Grid points for non-tradable component
gn = Param.gn;         % Grid points for permanent component
gTn = Param.gTn;    % Growth Rate Tradable Output
gNn = Param.gNn;    % Growth Rate Nontradable Output
r = Param.r; % Annual Interest Rate

g =  Param.g;  % Mean growth rate of permanent component
sigmag = Param.sigmag; % Std. dev growth rate of permanent component
Tsim = Param.Tsim;   % Simulation points
Tburn = Param.burn; % Burn-in period for simulation

%% Discretize VAR

fprintf("Discretize VAR ... \n")
[Pi3,S2,Xvec] = tpm(A, Sigma, [ytn;ynn;gn], Tsim, 0,[],0); % Fundamental Transition Matrix, Matrix of states and Simulation
fprintf("Done \n")

%% Transition Matrix for full information exercise

S = S2;
yt = unique(S(:,1)) ; % Transitory component of Tradable Endowment
yn = unique(S(:,2)) ; % Transitory component of Non-tradable Endowment
gt = unique(S(:,3)) ; % Common permanent component
S(:,3) = S(:,3) + g ; % This changes the mean of the growth rate to g

% Reduce the state-space to include only states that are reached in the
% simulation

SumPi=round(sum(Pi3,2)); 
Indd=find(SumPi==1);
Pi2=Pi3(any(Pi3,2),:);
Pi = Pi2(:,any(Pi3,1));
S = S(Indd,:) ;
Yn = length(S) ; % Size of the irreducible state-space

save('Data/TranMatFI.mat', "Pi", "Indd", "Xvec", "S",  "Yn", "yt", "yn", "gt",  '-v7.3')
fprintf("The transition matrix for the full information problem has been created  \n")



%% Create underlying series for observable series
% Xvec contains the simulated underlying states
% Xvec final order is: [yt, yn, gt, ns, gTT, gTN]

clear Pi Indd S yt yn gt Yn Pi2 SumPi   

fprintf("Starting computation of transition matrix for imperfect information ... \n")

Xvec = [[0;0;0] Xvec];
Xvec(4,2:end) = Xvec(1,2:end) - Xvec(1,1:end-1) + Xvec(3,2:end);   % Growth of Tradable Output: implicit equations given by the Kalman Filter 
Xvec(5,2:end) = Xvec(2,2:end) - Xvec(1,1:end-1) + Xvec(3,2:end);   % Implicit equations given by the Kalman Filter 

stdztT = std(Xvec(1,2:end));
stdztN = std(Xvec(2,2:end));
stdgt = std(Xvec(3,2:end));
stdgta = std(Xvec(4,2:end));   % Standard deviation of the growth rate of the tradable endowment
stdgna = std(Xvec(5,2:end));   % Standard deviation of the growth rate of the nontradable endowment

gT = linspace(-log(10)*stdgta,log(10)*stdgta,gTn)';  % Grid for the observable growth rate of tradable endowment
gN = linspace(-log(10)*stdgna,log(10)*stdgna,gNn)';  % Grid for the observable growth rate of Nontradable endowment

yt = unique(S2(:,1));    % Transitory component of the tradable endowment
yn = unique(S2(:,2));    % Transitory component of the nontradable endowment
gt = unique(S2(:,3));    % Common permanent component

S(:,1) = repmat(yt, [ynn*gn*gTn*gNn 1]) ;
S(:,2) = repmat(reshape(repmat(yn', [ytn 1]),[ytn*ynn 1]), [gn*gTn*gNn 1]) ;
S(:,3) = g + repmat(reshape(repmat(gt', [ytn*ynn 1]), [ynn*gn*ytn 1]), [gTn*gNn 1]) ;
S(:,4) = g + repmat(reshape(repmat(gT', [ytn*ynn*gn 1]), [ynn*gn*ytn*gTn 1]), [gNn 1]) ;
S(:,5) = g + repmat(reshape(repmat(gN', [ytn*ynn*gn*gTn 1]), [ynn*gn*ytn*gTn*gNn 1]),[1 1]) ;

% Kalman Filter Matrices

Z = [1 0 1 -1; ... 
     0 1 1 -1]; 
 
T = [A(1,1) A(1,2)    0    0 ; ...
     A(2,1) A(2,2)    0    0 ; ...
        0      0   A(3,3)  0 ; ...
        1      0      0    0 ] ;
    
cc = [0 0 (1-A(3,3))*g 0]' ; % Constant vector including the growth

R = [1 0 0 ; ...
     0 1 0 ; ...
     0 0 1 ; ...
     0 0 0 ]; ...
     
Q = [Sigma(1,1) Sigma(1,2) 0; Sigma(2,1) Sigma(2,2) 0 ; 0 0 Sigma(3, 3)]  ;
P = eye(size(T)) ;
err = 1 ;

fprintf("Kalman Filter Fix-point iteration ... \n")
while err>1e-13
    P1 = P ;
    P = T*P*T' - T*P*Z'*(Z*P*Z')^(-1)*(Z*P*T') + R*Q*R' ; % Ricatti Equation Fix Point Estimation
    err = max(abs(P(:) - P1(:))) ;
end
fprintf("Done \n")

k1 = (eye(size(P)) - P*Z'*(Z*P*Z')^(-1)*Z) ; % Kalman gains
k2 = P*Z'*(Z*P*Z')^(-1) ;                         % Kalman gains

fprintf("Re-estimate transition matrix ... \n")
[Pi,Indd]  = KF_transition_matrix(Xvec',  yt, yn, gt, gT, gN, T, g, cc, k1, k2, ynn, ytn, gn, gTn, gNn) ;  % This is the transition matrix given by the agent's beliefs
S = S(Indd,:);
fprintf("Done \n")

clear Param 
save('Data/TranMatII.mat', "Pi", "Indd", "Xvec", "S", "T", "cc", "yt", "yn", "gt", "gT", "gN", "k1", "k2", '-v7.3')
fprintf("Transition matrix for imperfect information created \n")
