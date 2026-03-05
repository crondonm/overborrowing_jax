%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: Use value function iteration to solve the social planner's
% problem under perfect information
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping

clearvars
clear global
close all

maxNumCompThreads(16);

fprintf("Starting to solve the Planner´s problem under perfect information ... \n")

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

yT = repmat(S(:,1),[bn bn]);  % Grid for the transitory tradable component
yN = repmat(S(:,2),[bn bn]);  % Grid for the transitory Non-tradable component
gT = repmat(S(:,3),[bn bn]);  % Grid for the common permanent component
bbb = repmat(reshape(repmat(b,[Yn 1]),[Yn*bn 1]),[1 bn]);
bbs = repmat(b,[Yn*bn 1]);

% Parameters for policy iteration
Delta  = 1;
Conv   = 1e-08;
window = Param.window;
Tol = 1e-08;
maxiter = 1000;


%% Value Function Iteration
%  *********************************************************************************

ct = (1+r)*bbb + exp(gT + yT) - bbs.*exp(gT);
c  = ((omega)*(max(ct,1e-10).^((eta-1)/eta)) + ((1-omega))*(exp(yN+gT).^((eta-1)/eta))).^((eta/(eta-1)));
p  = (((1-omega))/(omega)).*(exp(yN+gT)./max(ct,1e-10)).^(-1/eta); 
BC = - kappa*(exp(yT+gT)+p.*exp(yN+gT));
lambda = (omega).*(c.^(-rho+1/eta)).*(ct.^(-(1/eta)));

for ii=1:(Yn*bn)
    Temp1 = find(ct(ii,:)>0,1,'last');
    Temp2 = find((bbs(ii,:) - BC(ii,:))>=0,1,'first');
    if Temp2>Temp1
        BC(ii,Temp1) = bbs(ii,Temp1);
    end
end

% If consumption is negative, the constraint is binding: incentive so that the optimization never chooses this point
Ind = (bbs < BC) | (ct<0);
u = (c.^(1-rho))/(1-rho);
u(Ind==1) = -1e12; % Never choose a point where the collateral is violated

% Initialize Value Function

V= ones(Yn,bn);

% Start iteration

if min(c(:,end))<0
    fprintf('Unfeasible Grid')
end

fprintf('Start value function iteration ...  \n')

factor = beta.*exp(gT).^(1-rho);
for iter = 1:maxiter
    v = V;
    Vtemp = repmat(Pi*V,[bn 1]);
    [V,Pol] = max(u + factor.*Vtemp,[],2);
    V = reshape(V,[Yn bn]);
    Crit = max(abs(V(:) - v(:)));
    if mod(iter,10)==0
        disp([iter Crit])
    end
    if Crit<Tol
        disp([iter Crit])
        break
    end
    finiter=iter;
end

fprintf('Done \n')

% Reshaping policy functions

Pol  = reshape(Pol,[Yn bn]);
P3d  = reshape(p,[Yn bn bn]);
C3d  = reshape(c,[Yn bn bn]);
BC3d = reshape(BC,[Yn bn bn]);
CT3d = reshape(ct,[Yn bn bn]);

BCOptInd = sub2ind(size(BC3d),repmat((1:Yn)',[bn 1]),reshape(repmat((1:bn),[Yn 1]),[bn*Yn 1]),Pol(:));
UOpt     = reshape(lambda(BCOptInd),[Yn bn]);
COpt  = reshape(c(BCOptInd),[Yn bn]);
POpt = reshape(p(BCOptInd),[Yn bn]);
CTOpt  = reshape(ct(BCOptInd),[Yn bn]);
BCOpt = reshape(BC(BCOptInd),[Yn bn]);

%% Computing Optimal Tax

fprintf('Starting computation of optimal tax ... \n')

Bopt = zeros(Yn,bn);
Utom = zeros(Yn,bn);
EESP = zeros(Yn,bn);
ExpU = Pi*UOpt;
factorx = beta.*exp(gT).^(-rho) ;

for iz = 1:Yn
  for ib = 1:bn
    Bopt(iz,ib) = b(Pol(iz,ib));
    Utom(iz,ib) = ExpU(iz, Pol(iz,ib));
    EESP(iz,ib) = (1/((1+r)*factorx(iz))).*(UOpt(iz,ib));
  end
end

TAO =  (EESP./Utom) - 1;
AAA = Bopt <= BCOpt + (b(2)-b(1))/2;  
TAO(AAA==1) = 0;
TAO(TAO<0) = 0;

fprintf('Done \n')

%% Simulation

SimB      = zeros(1,Tsim);
SimB(1)   = floor(bn/2) ;

% Assign the value of sym but in the grid (using the closest one)

Temp1 = [findClosest2(S(:,1),yt) findClosest2(S(:,2),yn) findClosest2(S(:,3),gt+g)] ;
Sim = Xvec' ; % Simulation coming from the transition matrix
Simyt = Sim(:,1) ; Simyn = Sim(:,2) ; Simg = Sim(:,3) ; 

Index    = zeros(1,Tsim);
Index(1) = find((Temp1(:,1)==findClosest2(Sim(1,1),yt)).*(Temp1(:,2)==findClosest2(Sim(1,2),yn)).*(Temp1(:,3)==findClosest2(Sim(1,3),gt))) ;

ZBBT = nstd*std(Simyt) ; 
GBBT = nstd*std(Simg) ;

% Initizalizing variables 
 
PSim = zeros(1,Tsim);
CSim = zeros(1,Tsim);
BCSim = zeros(1,Tsim);
BCSim2 = zeros(1,Tsim);
CTSim = zeros(1,Tsim);
CNSim = zeros(1,Tsim);
TAOSim = zeros(1,Tsim);
SimBhat = zeros(1,Tsim);
grYTSim = zeros(1,Tsim);
grYNSim = zeros(1,Tsim);
Posterioryt = zeros(1,Tsim);
Posterioryn = zeros(1,Tsim);
Posteriorg  = zeros(1,Tsim);

fprintf("Starting Simulation ... \n")

for i=2:Tsim
    
    Posterioryt(i) = findClosest2(Simyt(i),yt);
    Posterioryn(i) = findClosest2(Simyn(i),yn);
    Posteriorg(i) = findClosest2(Simg(i),gt);
    Index(i) = find((Temp1(:,1)==Posterioryt(i)).*(Temp1(:,2)==Posterioryn(i)).*(Temp1(:,3)==Posteriorg(i)));
    SimB(i) = Pol(Index(i),SimB(i - 1));
    SimBhat(i) = b(SimB(i))/exp(Simyt(i - 1));     % ReNormalization of b
    BCSim(i) = BC3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1)))/exp(Simyt(i - 1));  % Collateral at period i, this is BC=Borrowing Constraint
    BCSim2(i) = BC3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1)));  % Collateral at period i, this is BC=Borrowing Constraint
    CTSim(i) = CT3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1)))/exp(Simyt(i - 1));  % Consumption of tradables 
    CSim(i) = C3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1)))/exp(Simyt(i - 1));  % Consumption of tradables 
    PSim(i) = P3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1)));                    % Price
    grYTSim(i) = exp(Simg(i) + Simyt(i) - Simyt(i - 1) + g);
    grYNSim(i) = exp(Simg(i) + Simyn(i) - Simyt(i - 1) + g);
    CNSim(i) = grYNSim(i).*exp(Simyt(i - 1));  % Consumption of Non-tradables 
    TAOSim(i) = TAO(Index(i), SimB(i-1));

   
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

fprintf("Done \n")

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
BCSim2 = BCSim2(1, burn + 1:end);
CTSim = CTSim(1, burn  + 1:end);
CNSim = CNSim(1, burn + 1:end);
TAOSim = TAOSim(1, burn + 1:end);
SimBhat = SimBhat(1, burn + 1:end);
grYTSim = grYTSim(1, burn + 1:end);
grYNSim = grYNSim(1, burn + 1:end);

Posteriorg  = Posteriorg(1, burn + 1:end);
Posterioryt = Posterioryt(1, burn + 1:end);
Posterioryn = Posterioryn(1, burn + 1:end);


%% Definition of crises 

fprintf("Starting crises analysis ... \n")

AAA = b(SimB) > (BCSim2 + (b(2)-b(1))/2) ; 
CCC = CA ; 
CCCT = nstd*std(CCC);
Crisis = (CCC>CCCT).*(1-AAA) ;
Crisis2 = sum(Crisis)/(length(CCC)) ;
Crisis = [0 Crisis] ; 
CrInd = find(Crisis==1) ;
CrInd = CrInd(CrInd>window + 1) ; 
CrInd = CrInd(CrInd< Tsim - burn - window) ;  


    for i=-5:5
       Cycles.IRB(i + window + 1,:) = SimBhat(CrInd + i - 1);
       Cycles.IRCA(i + window + 1,:) = CA(CrInd + i - 1) ;
       Cycles.IRBC(i + window + 1,:) = BCSim(CrInd + i - 1) ;
       Cycles.IRCT(i + window + 1,:) = CTSim(CrInd + i - 1) ;
       Cycles.IRCN(i + window + 1,:) = CNSim(CrInd + i - 2) ;
       Cycles.IRC(i + window + 1,:) = CSim(CrInd + i - 1);
       Cycles.IRP(i + window + 1,:) = PSim(CrInd + i - 1) ;
       Cycles.IRZ(i + window + 1,:) = yt(Posterioryt(CrInd + i - 2));
       Cycles.IRZN(i + window + 1,:) = yn(Posterioryn(CrInd + i));
       Cycles.IRG(i + window + 1,:) = gt(Posteriorg(CrInd + i - 2)) + g;       
       Cycles.IRDtoY(i + window + 1,:) = DtoY(CrInd + i); 
       Cycles.IRgrYTSim(i + window + 1,:) = grYTSim(CrInd + i - 2);
       Cycles.IRgrYNSim(i + window + 1,:) = grYNSim(CrInd + i - 2);
       Cycles.IRTAO(i + window + 1,:) = TAOSim(CrInd + i - 2);
       Cycles.IRCtoY(i + window + 1,:) = CtoY(CrInd + i);
       Cycles.IRCTtoY(i + window + 1,:) = CTtoY(CrInd + i - 1);
       Cycles.IRCNtoY(i + window + 1,:) = CNtoY(CrInd + i - 1);
       Cycles.IRCAtoY(i + window + 1,:) = CAtoY(CrInd + i - 1);
       Cycles.IRYtot(i + window + 1,:) = Ytot(CrInd + i - 1);
    end 
 
    Cycles.IRBMean = mean(Cycles.IRB, 2);
    Cycles.IRCAMean = mean(Cycles.IRCA, 2);
    Cycles.IRBCMean = mean(Cycles.IRBC, 2);
    Cycles.IRCTMean = mean(Cycles.IRCT, 2);
    Cycles.IRCMean = mean(Cycles.IRC, 2);
    Cycles.IRCNMean = mean(Cycles.IRCN, 2);
    Cycles.IRPMean = mean(Cycles.IRP, 2);
    Cycles.IRZMean = mean(Cycles.IRZ, 2);
    Cycles.IRZNMean = mean(Cycles.IRZN, 2);
    Cycles.IRGMean = mean(Cycles.IRG, 2);
    Cycles.IRDtoYMean = mean(Cycles.IRDtoY, 2);
    Cycles.IRgrYTSimMean = mean(Cycles.IRgrYTSim, 2);
    Cycles.IRgrYNSimMean = mean(Cycles.IRgrYNSim, 2);
    Cycles.IRTAOMean = mean(Cycles.IRTAO, 2);
    Cycles.IRCtoYMean = mean(Cycles.IRCtoY, 2);
    Cycles.IRCTtoYMean = mean(Cycles.IRCTtoY, 2);
    Cycles.IRCNtoYMean = mean(Cycles.IRCNtoY, 2);
    Cycles.IRCAtoYMean = mean(Cycles.IRCAtoY, 2);
    Cycles.IRYtotMean = mean(Cycles.IRYtot, 2);

fprintf("Done \n")

%% Save Structure

fprintf("Saving ... \n")

FIP.S = S ;
FIP.V = V;
FIP.b = b ;
FIP.yt = yt ;
FIP.yn = yn ;
FIP.gt = gt ;
FIP.Pi  = Pi ;
FIP.mCT = mCT;
FIP.mCN = mCN;
FIP.mP = mP;
FIP.CA = CA;
FIP.Sim = Sim ;
FIP.Pol = Pol ;
FIP.AAA = AAA ;
FIP.CCC = CCC ;
FIP.TAO = TAO;
FIP.TAOSim = TAOSim;


FIP.grYTSim = grYTSim;
FIP.grYNSim = grYNSim;

FIP.POpt = POpt;
FIP.COpt = COpt;
FIP.BCOpt = BCOpt;
FIP.CTOpt = CTOpt;
FIP.SimB = SimB;
FIP.SimBhat = SimBhat;
FIP.Index = Index;
FIP.PSim = PSim;
FIP.CSim = CSim;
FIP.CTSim = CTSim;
FIP.CNSim = CNSim;
FIP.BCSim = BCSim;
FIP.BCSim2 = BCSim2;
FIP.Cycles = Cycles;
FIP.Posterioryt =  Posterioryt;
FIP.Posteriorg = Posteriorg;
FIP.Posterioryn = Posterioryn;
FIP.CrInd = CrInd;
FIP.Freq = Crisis2
FIP.DtoY = DtoY;
FIP.mDtoY = mean(DtoY)
FIP.CtoY = CtoY;
FIP.mCtoY = mean(CtoY);
FIP.CAtoY = CAtoY;
FIP.CTtoY = CTtoY;
FIP.CNtoY = CNtoY;
FIP.mCTtoY = mean(CTtoY);
FIP.mCNtoY = mean(CNtoY);
FIP.stdZBB = ZBBT;
FIP.stdGBB = GBBT;

save('Data/FIP.mat', 'FIP', '-v7.3')
save("Data/FIPsim.mat",'-struct','FIP', '-v7.3') ;


