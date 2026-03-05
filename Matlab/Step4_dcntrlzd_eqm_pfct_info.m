%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: Use policy function iteration to solve the decentralized equilibrium
%                    Step 1: Conjecture B' and use budget constraint and price equation.
%                    Step 2: Assume the constraint is binding. Set B'_{k-1} = kappa(yt + P_k yn)
%                    Step 3: Compute Expected value of lambda
%                    Step 4: Simulate
%                    Step 5: Compute IRFs              
%
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
%
% Date: March 2025
%
% Standard Convention:
% Variable named using Upper case: Variable only in terms of the original state 
% Variable named using Lower case: Variable depends on the decision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Housekeeping

clearvars
clear global
close all

maxNumCompThreads(16);

fprintf("Starting to solve the decentralized eqm under perfect information ... \n")

% Load parameters

load('Data/Param.mat')
fprintf("Parameters loaded... \n")

% Load transition matrix and sthocastic states

load('Data/TranMatFI.mat')
fprintf("Transition Matrix loaded... \n")

% Grid and stochastic parameters

A = Param.A;   % Autocorrelation matrix
Yn = length(S);
ytn = Param.ytn;       % Grid points for tradable component
ynn = Param.ynn;      % Grid points for non-tradable component
gn = Param.gn;         % Grid points for permanent component

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

% Defining grids as matrices

yT = repmat(S(:,1),[bn bn]) ;  % Grid for the transitory tradable component
yN = repmat(S(:,2),[bn bn]) ;  % Grid for the transitory Non-tradable component
gT = repmat(S(:,3),[bn bn]) ;  % Grid for the common permanent component
bbb = repmat(reshape(repmat(b,[Yn 1]),[Yn*bn 1]),[1 bn]) ; % Debt today
bbs = repmat(b,[Yn*bn 1]) ; % Decision today of debt Tomorrow 
yInd   = repmat((1:Yn)',[1 bn]) ;
bbsInd = repmat((1:bn),[Yn*bn 1]) ;

% Parameters for policy iteration
Delta  = 1 ;
Conv   = 1e-08 ;
window = Param.window;
Tol = 1e-08 ;
maxiter = 1000;


%% Policy Function Iteration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_conj = reshape(bbb(:,1),[Yn bn]);  % First Step: Conjecture for debt (B' = B)

%  Uncomment the following lines if you use a conjecture that includes
%  points outside the grid

    % IndBk1 = B_conj < b(1);               % Index where the borrowing is to the left of the grid b(1) is the minimum value of debt in the grid
    % B_conj(IndBk1) = b(1);                % Impose the minimum value of debt where the conjecture is out of the grid
    % B_conj_ind = nan(size(B_conj)) ;
    % b1=b;
    % parfor ii = 1:(Yn*bn)                 % This step is to make sure that the conjecture is consistent with the points in the grid, not relevant when the conjecture is the grid itself
    %     B_conj_ind(ii) = find(b1 <= B_conj(ii), 1, 'last' );
    % end

Lambda_conj = ones(Yn,bn) ; % Lambda is the marginal utility of consumption, we start a point where consumption is one
IndLambda = sub2ind(size(Lambda_conj),repmat((1:Yn)',[bn*bn 1]),bbsInd(:)) ; % Changes the subs to indexes, useful to compute the euler equation

% First Order Conditions
fprintf("Creating Matrices ... \n")

CT =  (1+r)*bbb + exp(gT + yT) - bbs.*exp(gT) ;
Pr = (((1-omega))/(omega)).*(exp(yN+gT)./max(CT,1e-10)).^(-1/eta) ;
C  = ((omega)*(max(CT,1e-10).^((eta-1)/eta)) + ((1-omega))*(exp(yN+gT).^((eta-1)/eta))).^((eta/(eta-1))) ;
U  = (C.^(1-rho))/(1-rho) ;
lambda = (omega).*(C.^(-rho+1/eta)).*(CT.^(-(1/eta))) ; % Marginal utility of consumption for each point in the grid
lambda(CT<0) = 1e+10 ; % For each negative value of consumption the marginal utility increasing consumption is infinity
M = - kappa*(exp(yT+gT)+Pr.*exp(yN+gT)) ;

for ii=1:(Yn*bn)  
    Temp1 = find(CT(ii,:)>0,1,'last') ;
    Temp2 = find((bbs(ii,:) -M(ii,:))>=0,1,'first') ;
    if Temp2>Temp1
       M(ii,Temp1) = bbs(ii,Temp1) ;
    end
end

fprintf("Done \n")

% Policy Iteration

fprintf("Starting Policy Iteration ... \n")

iter = 1 ;

while Delta > Conv
    
    % Construct your conjecture:
    
    Lambdalast = Lambda_conj ; % Best guess for Lambda
    Blast = B_conj ;           % Best guess for debt
       
    % Marginal Utility of Consumption:
    Lambda = reshape(Lambda_conj(IndLambda),[Yn*bn bn]) ; % Marginal utility of consumption that I would have tomorrow if I have a given bbs
    
    % Now compute the Euler Equation:
    % mu = lagrange multiplier with respect to the collateral
    mu = 1 - (1+r)*beta*(exp(gT).^(-rho)).*reshape(Pi*reshape(Lambda,[Yn bn*bn]),[Yn*bn bn])./lambda ; 
    Elambda = reshape(Pi*reshape(Lambda,[Yn bn*bn]),[Yn*bn bn]);
    % Step 3: Compute Expected value of Lambda:
    mu2 = mu ; 
    mu2(CT < 0) = 1e+10 ;        % Marginal utility of relaxing the borrowing constraint when consumption is negative: infty   
    mu2((bbs - M) < 0) = 1e+10 ; % If the constraint is violated, marginal utility of relaxing the borrowing constraint is infty
    muabs = abs(mu2);            
    [muopt,muInd] = min(muabs,[],2) ;
    % Index for the optimal debt given muopt
    B_conj_ind = reshape(muInd,[Yn bn]) ;  
    lambda3d = reshape(lambda,[Yn bn bn]) ; % 3-dimensional object with yn, b and bprime
    % I am choosing the index for the optimal value of b given the optimal value for mu
    Ind2 = sub2ind(size(lambda3d),repmat((1:Yn)',[bn 1]),reshape(repmat((1:bn),[Yn 1]),[Yn*bn 1]),muInd) ;  
    Lambda_conj = reshape(lambda(Ind2),[Yn bn]) ; % Assign the new conjecture
    B_conj = b(B_conj_ind) ;
    Delta1 = max(abs(Lambda_conj(:) - Lambdalast(:))) ; 
    Delta2 = max(abs(B_conj(:) - Blast(:))) ; 
    Delta = max(Delta1,Delta2) ;
    disp([iter, Delta, Delta1, Delta2])
    iter=iter+1;
end

fprintf("Done \n")

% Policy Function:

Pol = B_conj_ind ;
Lambda3d = reshape(Elambda, [Yn bn bn]);
BC3d = reshape(M,[Yn bn bn]) ;
P3d  = reshape(Pr,[Yn bn bn]) ;
CT3d = reshape(CT,[Yn bn bn]) ;
C3d  = reshape(C,[Yn bn bn]) ;

BCOptInd = sub2ind(size(BC3d),repmat((1:Yn)',[bn 1]),reshape(repmat((1:bn),[Yn 1]),[bn*Yn 1]),Pol(:)) ;
BCOpt = reshape(M(BCOptInd),[Yn bn]) ;
POpt = reshape(Pr(BCOptInd),[Yn bn]) ;

CTOpt = reshape(CT(BCOptInd),[Yn bn]) ;
COpt = reshape(C(BCOptInd),[Yn bn]) ;
UOpt = reshape(U(BCOptInd),[Yn bn]) ;

%% Computing Value Function

fprintf("Starting computation of value function ... \n")

factor = beta.*exp(gT).^(1-rho) ;

% Initialize  iteration

V = zeros(Yn,bn);

for iter = 1:maxiter
    v = V;
    ExpV = Pi*V;
    for iz = 1:Yn
        for ib = 1:bn
            V(iz,ib) = UOpt(iz,ib)+ factor(iz).*ExpV(iz, Pol(iz,ib));
        end
    end
    Crit = max(abs(V(:) - v(:))) ;
   if mod(iter,10)==0
      disp([iter Crit])
   end
   if Crit<Tol
       disp([iter Crit])
       break
   end
end

fprintf("Done \n")

%% Simulation

SimB    = zeros(1,Tsim);
SimB(1) = floor(bn/2) ;
Temp1 = [findClosest2(S(:,1), yt) findClosest2(S(:,2), yn) findClosest2(S(:,3), gt + g)] ;
Sim = Xvec' ; % Simulation coming from the transition matrix
Simyt = Sim(:,1) ; Simyn = Sim(:,2) ; Simg = Sim(:,3) ; % assigning the values for each endowment and shock
IndApp = [findClosest2(Simyt, yt) findClosest2(Simyn, yn) findClosest2(Simg, gt)]; % Assign the value of sym but in the grid (using the closest one)
Index = zeros(1, Tsim);
Index(1) = find((Temp1(:,1)==findClosest2(Sim(1,1), yt)).*(Temp1(:,2)==findClosest2(Sim(1,2), yn)).*(Temp1(:,3)==findClosest2(Sim(1,3), gt))) ;
ZBBT = std(Simyt) ; % Standard Deviation to consider a Boom Bust
GBBT = std(Simg) ;
 
% Initizalizing variables 
 
CSim = zeros(1, Tsim);
VSim = zeros(1, Tsim);
PSim = zeros(1, Tsim);
BCSim = zeros(1, Tsim);
CTSim = zeros(1, Tsim);
CNSim = zeros(1, Tsim);
SimBhat = zeros(1, Tsim);
grYTSim = zeros(1, Tsim);
grYNSim = zeros(1, Tsim);
LambdaSim = zeros(1, Tsim);
Posterioryt = zeros(1, Tsim);
Posterioryn = zeros(1, Tsim);
Posteriorg = zeros(1, Tsim);

% Starting simulation

fprintf("Starting Simulation ... \n")

for i=2:Tsim

    Posterioryt(i) = findClosest2(Simyt(i), yt); 
    Posterioryn(i) = findClosest2(Simyn(i), yn);
    Posteriorg(i) = findClosest2(Simg(i), gt);
    Index(i) = find((Temp1(:,1)==Posterioryt(i)).*(Temp1(:,2)==Posterioryn(i)).*(Temp1(:,3)==Posteriorg(i))) ;
    SimB(i) = Pol(Index(i),SimB(i - 1)) ;
    SimBhat(i) = b(SimB(i))/exp(Simyt(i - 1)) ; % ReNormalization of b
    BCSim(i) = BC3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1)))/exp(Simyt(i - 1)) ; % Collateral at period i, this is BC=Borrowing Constraint
    CTSim(i) = CT3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1)))/exp(Simyt(i - 1)) ; % Consumption of tradables 
    CSim(i) = C3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1)))/exp(Simyt(i - 1)) ; % Total Consumption
    PSim(i) = P3d(Index(i ), SimB(i - 1), Pol(Index(i), SimB(i - 1))) ; % Price
    LambdaSim(i) = Lambda3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1)))*exp(Simyt(i - 1)) ;    % Price
    grYTSim(i) = exp(Simg(i) + Simyt(i) - Simyt(i - 1) + g);
    grYNSim(i) = exp(Simg(i) + Simyn(i) - Simyt(i - 1) + g);
    CNSim(i) = grYNSim(i).*exp(Simyt(i - 1));

    if mod(i,Tsim/10) == 0
        disp(i);
    end
end



%% Definition of series for aggregate variables

Ytot = ((exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim(2:end).*exp(Simyn(2:end)' + Simg(2:end)' + g)))./exp(Simyt(1:end- 1)');
CtoY = CSim(2:end).*exp(Simyt(1:end- 1)')./(exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim(2:end).*exp(Simyn(2:end)'+Simg(2:end)'+g));
CTtoY = CTSim(2:end).*exp(Simyt(1:end- 1)')./(exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim(2:end).*exp(Simyn(2:end)'+Simg(2:end)'+g));
CNtoY = CNSim(2:end).*exp(Simyt(1:end- 1)')./(exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim(2:end).*exp(Simyn(2:end)'+Simg(2:end)'+g));
DtoY = SimBhat(2:end).*exp(Simyt(1:end-1)')./(exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim(2:end).*exp(Simyn(2:end)'+Simg(2:end)'+g)) ; 
CA = (SimBhat(2:end) - SimBhat(1:end-1)) ;
CAtoY = CA.*exp(Simyt(1:end-1)')./(exp(Simyt(2:end)'+Simg(2:end)'+g) + PSim(2:end).*exp(Simyn(2:end)'+Simg(2:end)' + g)) ; 
mP = mean(PSim);
mCT = mean(CTSim);
mCN = mean(CNSim);

fprintf("Done  \n")

%% Making series of equal length

Simyt = Sim(burn + 1:end, 1);
Simyn = Sim(burn + 1:end, 2);
Simg = Sim(burn + 1:end, 3);

CA = CA(1, burn:end);
CAtoY = CAtoY(1, burn:end);
CtoY = CtoY(1, burn:end);
CTtoY = CTtoY(1, burn:end);
CNtoY = CNtoY(1, burn:end);
DtoY = DtoY(1, burn:end);
Ytot = Ytot(1, burn:end);

SimB = SimB(1, burn + 1:end);
CSim = CSim(1, burn + 1:end);
PSim = PSim(1, burn + 1:end);
BCSim = BCSim(1, burn + 1:end);

CTSim = CTSim(1, burn  + 1:end);
CNSim = CNSim(1, burn + 1:end);
SimBhat = SimBhat(1, burn + 1:end);
grYTSim = grYTSim(1, burn + 1:end);
grYNSim = grYNSim(1, burn + 1:end);

Posteriorg  = Posteriorg(1, burn + 1:end);
Posterioryt = Posterioryt(1, burn + 1:end);
Posterioryn = Posterioryn(1, burn + 1:end);

%% Definition of crises 

fprintf("Starting crises analysis ... \n")

AAA    = SimBhat > (BCSim + (b(2) - b(1))/2) ; 
CCC    = CA; 
CCCT   = std(CCC);
Crisis = (CCC>CCCT).*(1-AAA) ;
Crisis2 = sum(Crisis)/(length(CCC)) ;
%Crisis  = [0 Crisis] ; 
CrInd = find(Crisis==1) ;
CrInd = CrInd(CrInd > window + 1 ) ; 
CrInd = CrInd(CrInd< Tsim - burn - window) ;  

for i=-window : window
    Cycles.IRB(i + window + 1,:) = SimBhat(CrInd + i);
    Cycles.IRCA(i + window + 1,:) = CA(CrInd + i) ;
    Cycles.IRBC(i + window + 1,:) = BCSim(CrInd + i) ;
    Cycles.IRC(i + window + 1,:) = CSim(CrInd + i) ;
    Cycles.IRCT(i + window + 1,:) = CTSim(CrInd + i) ;
    Cycles.IRCN(i + window + 1,:) = CNSim(CrInd + i) ;
    Cycles.IRP(i + window + 1,:) = PSim(CrInd + i) ;
    Cycles.IRZ(i + window + 1,:) = yt(Posterioryt(CrInd + i - 1)) ;
    Cycles.IRZN(i + window + 1,:) = yn(Posterioryn(CrInd + i + 1)) ;
    Cycles.IRG(i + window + 1,:) = gt(Posteriorg(CrInd + i - 1)) + g;       
    Cycles.IRDtoY(i + window + 1,:) = DtoY(CrInd + i + 1); 
    Cycles.IRLambda(i + window + 1,:) = LambdaSim(CrInd + i);
    Cycles.IRgrYTSim(i + window + 1,:) = grYTSim(CrInd + i);
    Cycles.IRgrYNSim(i + window + 1,:) = grYNSim(CrInd + i);
    Cycles.IRCtoY(i + window + 1,:) = CtoY(CrInd + i + 1);
    Cycles.IRCAtoY(i + window + 1,:) = CAtoY(CrInd + i);
    Cycles.IRCTtoY(i + window + 1,:) = CTtoY(CrInd + i);
    Cycles.IRCNtoY(i + window + 1,:) = CNtoY(CrInd + i);
    Cycles.IRYtot(i + window + 1,:) = Ytot(CrInd + i);
end 
 
Cycles.IRBMean = mean(Cycles.IRB, 2);
Cycles.IRCAMean = mean(Cycles.IRCA,2);
Cycles.IRBCMean = mean(Cycles.IRBC,2);
Cycles.IRCMean = mean(Cycles.IRC,2);
Cycles.IRCTMean = mean(Cycles.IRCT,2);
Cycles.IRCNMean = mean(Cycles.IRCN,2);
Cycles.IRPMean = mean(Cycles.IRP, 2);
Cycles.IRZMean = mean(Cycles.IRZ, 2);
Cycles.IRZNMean = mean(Cycles.IRZN, 2);
Cycles.IRGMean = mean(Cycles.IRG, 2);
Cycles.IRDtoYMean = mean(Cycles.IRDtoY,2);
Cycles.IRCtoYMean = mean(Cycles.IRCtoY,2);
Cycles.IRCAtoYMean = mean(Cycles.IRCAtoY,2);
Cycles.IRCTtoYMean = mean(Cycles.IRCTtoY,2);
Cycles.IRCNtoYMean = mean(Cycles.IRCNtoY,2);
Cycles.IRLambdaMean = mean(Cycles.IRLambda,2);
Cycles.IRgrYTSimMean = mean(Cycles.IRgrYTSim,2);
Cycles.IRgrYNSimMean = mean(Cycles.IRgrYNSim,2);
Cycles.IRYtotMean = mean(Cycles.IRYtot,2);

fprintf("Done \n")

%% Output

fprintf("Saving ... \n")

FICE.V = V;
FICE.b = b;
FICE.S = S; 
FICE.yt = yt;
FICE.yn = yn;
FICE.gt = gt;
FICE.Pi = Pi;
FICE.mCT = mCT;
FICE.mCN = mCN;
FICE.mP = mP;
FICE.CA = CA;
FICE.Sim = Sim;
FICE.Pol = Pol;
FICE.AAA = AAA;
FICE.CCC = CCC;
FICE.UOpt = UOpt;
FICE.COpt = COpt;
FICE.SimB = SimB;
FICE.PSim = PSim;
FICE.CSim = CSim;
FICE.Freq  = Crisis2
FICE.Ytot = Ytot;
FICE.DtoY = DtoY;
FICE.CtoY = CtoY;
FICE.CAtoY = CAtoY;
FICE.CTtoY = CTtoY;
FICE.CNtoY = CNtoY;
FICE.mCtoY = mean(CtoY);
FICE.mDtoY = mean(DtoY)
FICE.mCTtoY = mean(CTtoY);
FICE.mCNtoY = mean(CNtoY);
FICE.CTOpt = CTOpt;
FICE.CTSim = CTSim;
FICE.BCSim = BCSim;
FICE.Cycles = Cycles;
FICE.CrInd = CrInd;
FICE.grYTSim = grYTSim;
FICE.grYNSim = grYNSim;
FICE.Index = Index;
FICE.stdZBB = ZBBT;
FICE.stdGBB = GBBT;
FICE.SimBhat = SimBhat;
FICE.LambdaSim = LambdaSim;

save('Data/FICE.mat','FICE', '-v7.3') ;
save("Data/FICEsim.mat",'-struct','FICE', '-v7.3') ;

