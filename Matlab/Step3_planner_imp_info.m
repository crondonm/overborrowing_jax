%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: Use value function iteration to solve the social planner's
% problem under imperfect information
%               
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
%
% Date: March 2025
%
% Standard Convention:
% Variable named using Upper case: Variable only in terms of the original state 
% Variable named using Lower case: Variable depends on the decision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping

clearvars
clear global
close all

%maxNumCompThreads(16);

% Load parameters

fprintf("Starting to solve the Planner´s problem under imperfect information ... \n")

load('Data/Param.mat')
fprintf("Parameters loaded... \n")

% Load transition matrix and sthocastic states

load('Data/TranMatII.mat')
fprintf("Transition Matrix loaded... \n")

% Grid and stochastic parameters

A = Param.A;   % Autocorrelation matrix
Yn = length(S);
ytn = Param.ytn;       % Grid points for tradable component
ynn = Param.ynn;      % Grid points for non-tradable component
gn = Param.gn;         % Grid points for permanent component
gTn = Param.gTn;    % Growth Rate Tradable Output
gNn = Param.gNn;    % Growth Rate Nontradable Output

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

gTT = repmat(S(:,4),[bn bn]);
gNN = repmat(S(:,5),[bn bn]);
bbb = repmat(reshape(repmat(b,[Yn 1]),[Yn*bn 1]),[1 bn]);
bbs = repmat(b,[Yn*bn 1]);

% Parameters for value function iteration

Delta  = 1;
Conv   = 1e-08;
window = Param.window;
Tol = 1e-08;
finiter = 0;
maxiter = 1000;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Value Function Iteration
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First Order Conditions.

ct     = (1+r)*bbb + exp(gTT) - bbs.*exp(gTT);
c      = ((omega)*(max(ct,1e-10).^((eta-1)/eta)) + ((1-omega))*(exp(gNN).^((eta-1)/eta))).^((eta/(eta-1)));
p      = (((1-omega))/(omega)).*(exp(gNN)./max(ct,1e-10)).^(-1/eta); 
BC     = - kappa*(1+p.*exp(gNN)./exp(gTT));
lambda = (omega).*(c.^(-rho+1/eta)).*(ct.^(-(1/eta)));

for ii=1:(Yn*bn)
    Temp1 = find(ct(ii,:)>0,1,'last');
    Temp2 = find((bbs(ii,:) -BC(ii,:))>=0,1,'first');
    if Temp2>Temp1
        BC(ii,Temp1) = bbs(ii,Temp1);
    end
end

% If consumption is negative, the constraint is binding: incentive so that the optimization never chooses this point
Ind = (bbs < BC) | (ct<0);
u = (c.^(1-rho))/(1-rho);
u(Ind==1) = -1e12;

if min(c(:,end))<0
    fprintf('Unfeasible Grid')
end


% Initialize Value Function
V= ones(Yn,bn)./sqrt(15);


fprintf('Start value function iteration ...  \n')

factor = beta.*exp(gTT).^(1-rho);
for iter = 1:maxiter
    v = V;
    Vtemp = repmat(Pi*V,[bn 1]);
    [V,Pol] = max(u + factor.*Vtemp,[],2);
    V = reshape(V,[Yn bn]);
    Crit = max(abs(V(:) - v(:)));
    if mod(iter,10)==0
        disp([iter Crit])
    end
    finiter=iter;
    if Crit<Tol
        disp([iter Crit])
        break
    end
end

fprintf('Done \n')

% Reshaping policy functions

Pol  = reshape(Pol,[Yn bn]);
BC3d = reshape(BC,[Yn bn bn]);
CT3d = reshape(ct,[Yn bn bn]);
C3d  = reshape(c,[Yn bn bn]);
P3d  = reshape(p,[Yn bn bn]);

BCOptInd = sub2ind(size(BC3d),repmat((1:Yn)',[bn 1]),reshape(repmat((1:bn),[Yn 1]),[bn*Yn 1]),Pol(:));
CTOpt    = reshape(ct(BCOptInd),[Yn bn]);
COpt     = reshape(c(BCOptInd),[Yn bn]);
POpt     = reshape(p(BCOptInd),[Yn bn]);
BCOpt   = reshape(BC(BCOptInd),[Yn bn]);
UOpt     = reshape(lambda(BCOptInd),[Yn bn]);

%% Computing Optimal Tax

fprintf('Starting computation of optimal tax ... \n')

Bopt = zeros(Yn,bn);
Utom = zeros(Yn,bn);
EESP = zeros(Yn,bn);
ExpU = Pi*UOpt;
factorx = (1 + r) * beta .* exp(gTT).^( -rho);

for iz = 1:Yn
  for ib = 1:bn
            Bopt(iz,ib) = b(Pol(iz,ib));
            Utom(iz,ib) = ExpU(iz, Pol(iz,ib));
            EESP(iz,ib) = (1 / (factorx(iz))) .* (UOpt(iz,ib));
  end
end

TAO =  (EESP./Utom) - 1;
AAA = Bopt <= BCOpt  + (b(2) - b(1))/2 ; 
TAO(AAA==1) = 0;
TAO(TAO<0) = 0;


fprintf('Done \n')


%% Simulation

SimB = zeros(1,Tsim);
SimB(1) = floor(bn/2) ;
Temp1 = [findClosest2(S(:,1),yt) findClosest2(S(:,2),yn) findClosest2(S(:,3),gt+g) findClosest2(S(:,4),gT+g) findClosest2(S(:,5),gN+g)] ;
Sim = Xvec' ;
Simyt = Sim(:,1) ; Simyn = Sim(:,2) ; Simg = Sim(:,3) ; 
SimgT = Sim(:,4) ; SimgN = Sim(:,5) ;

a = [Sim(1,1); Sim(1,2); g+Sim(1,3); (Sim(1,1) + Sim(1,3) - Sim(1,4)+g)]; % Posterior for the first period is a (which is correct)
IndApp = [findClosest2(Simyt,yt) findClosest2(Simyn,yn) findClosest2(Simg,gt) findClosest2(SimgT,gT) findClosest2(SimgN,gN)];
Index = zeros(1,Tsim);
Index(1) = find((Temp1(:,1)==findClosest2(Sim(1,1),yt)).*(Temp1(:,2)==findClosest2(Sim(1,2),yn)).*(Temp1(:,3)==findClosest2(Sim(1,3),gt)).*(Temp1(:,4)==findClosest2(Sim(1,4),gT)).*(Temp1(:,5)==findClosest2(Sim(1,5),gN))) ;

% Criteria to consider a Boom-Bust Cycle in exogenous parameters

ZBBT = nstd*std(Simyt) ; % Standard Deviation to consider a Boom Bust
GBBT = nstd*std(Simg) ;

% Initizalizing variables 
 
TAOSim = zeros(1,Tsim + 1);
CSim = zeros(1,Tsim + 1);
PSim = zeros(1,Tsim + 1);
ghat = zeros(1,Tsim + 1);
ythat = zeros(1,Tsim + 1);
ynhat = zeros(1,Tsim + 1);
BCSim = zeros(1,Tsim + 1);
CTSim = zeros(1,Tsim + 1);
CNSim = zeros(1,Tsim + 1);
Posteriorg = zeros(1,Tsim + 1);
Posterioryt = zeros(1,Tsim + 1);
Posterioryn = zeros(1,Tsim + 1);
LambdaSim = zeros(1,Tsim + 1);
Simyt = zeros(1,Tsim + 1);
Simyn = zeros(1,Tsim + 1);
Simg = zeros(1, Tsim + 1);
grYTSim = zeros(1, Tsim + 1);
grYNSim = zeros(1, Tsim + 1);

fprintf("Starting Simulation ... \n")

for i=2:Tsim+1
    att1 = T*a(:,i-1) + cc ;
    a(:,i) = k1*att1 + k2*[(SimgT(i) + g);(SimgN(i) + g)] ; % Posterior is formed using the kalman weights
    ythat(i) = a(1,i) ;
    ynhat(i) = a(2,i) ;
    ghat(i)  = a(3,i) ;     
    GaT = IndApp(i,4) ;
    GaN = IndApp(i,5) ;
    Posterioryt(i) = findClosest2(ythat(i), yt) ;
    Posterioryn(i) = findClosest2(ynhat(i), yn) ;
    Posteriorg(i) = findClosest2(ghat(i) - g, gt) ;
    Simyt(i) = Sim(i, 1);
    Simyn(i) = Sim(i, 2);
    Simg(i) = Sim(i, 3);
    grYTSim(i) = SimgT(i);
    grYNSim(i) = SimgN(i);
    Index(i) = find((Temp1(:,1)==Posterioryt(i)).*(Temp1(:,2)==Posterioryn(i)).*(Temp1(:,3)==Posteriorg(i)).*(Temp1(:,4)==GaT).*(Temp1(:,5)==GaN));
    SimB(i) = Pol(Index(i), SimB(i - 1)) ;
    BCSim(i) = BC3d(Index(i), SimB(i - 1), Pol(Index(i),SimB(i - 1))) ;
    PSim(i) = P3d (Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1))) ;
    CSim(i) = C3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1))) ;
    CTSim(i) = CT3d(Index(i), SimB(i - 1), Pol(Index(i), SimB(i - 1))) ;
    CNSim(i) = exp(SimgN(i) + g).*exp(Simyt(i-1));
    TAOSim(i) = TAO(Index(i), SimB(i-1));
    
    if mod(i,Tsim/10) == 0
       disp(i);
    end
 end

CSim = CSim(2:end); PSim = PSim(2:end);
CNSim = CNSim(2:end); CTSim = CTSim(2:end);
BCSim = BCSim(2:end); SimB = SimB(2:end); 
Index = Index(2:end);
Simyt = Simyt(2:end); Simyn = Simyn(2:end); Simg = Simg(2:end);
Posteriorg = Posteriorg(2:end); Posterioryn = Posterioryn(2:end); Posterioryt = Posterioryt(2:end);
grYNSim = grYNSim(2:end); grYTSim = grYTSim(2:end); 
TAOSim = TAOSim(1, 2:end);

%% Definition of series for aggregate variables

Ytot = ((exp(Simyt(2:end)+Simg(2:end)+g) + PSim(2:end).*exp(Simyn(2:end) + Simg(2:end) + g)))./exp(Simyt(1:end- 1));
CtoY = CSim(2:end).*exp(Simyt(1:end- 1))./(exp(Simyt(2:end)+Simg(2:end) + g) + PSim(2:end).*exp(Simyn(2:end) + Simg(2:end) + g));
CTtoY = CTSim(2:end).*exp(Simyt(1:end- 1))./(exp(Simyt(2:end)+Simg(2:end)+g) + PSim(2:end).*exp(Simyn(2:end)+Simg(2:end) + g));
CNtoY = CNSim(2:end).*exp(Simyt(1:end- 1))./(exp(Simyt(2:end)+Simg(2:end)+g) + PSim(2:end).*exp(Simyn(2:end)+Simg(2:end) + g));
DtoY = b(SimB(2:end)).*exp(Simyt(1:end-1))./(exp(Simyt(2:end)+Simg(2:end)+g) + PSim(2:end).*exp(Simyn(2:end)+Simg(2:end) + g)); 
CA = (b(SimB(2:end)) -  b(SimB(1:end-1)));
CAtoY = CA.*exp(Simyt(1:end-1))./(exp(Simyt(2:end)+Simg(2:end)+g) + PSim(2:end).*exp(Simyn(2:end)+Simg(2:end) + g)) ; 
mP = mean(PSim);
mCT = mean(CTSim);
mCN = mean(CNSim);

CSim = CSim(2:end); PSim = PSim(2:end);
CNSim = CNSim(2:end); CTSim = CTSim(2:end);
BCSim = BCSim(2:end); SimB = SimB(2:end); 
Index = Index(2:end);
Simyt = Simyt(2:end); Simyn = Simyn(2:end); Simg = Simg(2:end);
Posteriorg = Posteriorg(2:end); Posterioryn = Posterioryn(2:end); Posterioryt = Posterioryt(2:end);
grYNSim = grYNSim(2:end); grYTSim = grYTSim(2:end); 
TAOSim = TAOSim(1, 2:end);


%% Making series of equal length

Simyt = Simyt(1, burn:end);
Simyn = Simyn(1, burn:end);
Simg = Simg(1, burn:end);
SimB = SimB(1, burn:end);
CSim = CSim(1, burn:end);
PSim = PSim(1, burn:end);
CTSim = CTSim(1, burn:end);
CNSim = CNSim(1, burn:end);
BCSim = BCSim(1, burn:end);
Posteriorg = Posteriorg(1, burn:end);
Posterioryt = Posterioryt(1, burn:end);
Posterioryn = Posterioryn(1, burn:end);

Ytot = Ytot(1, burn:end);
CA = CA(1, burn:end);
CAtoY = CAtoY(1, burn:end);
CtoY = CtoY(1, burn:end);
CTtoY = CTtoY(1, burn:end);
CNtoY = CNtoY(1, burn:end);
DtoY = DtoY(1, burn:end);
TAOSim = TAOSim(1, burn:end);

fprintf("Done \n")
  
%% Definition of crises 

fprintf("Starting crises analysis ... \n")

AAA = b(SimB) > BCSim + (b(2) - b(1))/2 ; 
CCC =CA;
CCCT = nstd*std(CCC);
Crisis = (CCC > CCCT).*(1 - AAA) ;
Crisis2 = sum(Crisis)/(length(CCC)) ;
CrInd = find(Crisis == 1) ;
CrInd = CrInd(CrInd > window + 1) ; 
CrInd = CrInd(CrInd < Tsim - burn - window) ;  


for i= - window : window
    Cycles.IRB(i + window + 1, :) = b(SimB(CrInd + i));
    Cycles.IRCA(i + window  + 1,:) = CA(CrInd + i);
    Cycles.IRBC(i + window + 1,:) = BCSim(CrInd + i);
    Cycles.IRC(i + window + 1,:) = CSim(CrInd + i);
    Cycles.IRCT(i + window + 1,:) = CTSim(CrInd + i);
    Cycles.IRCN(i + window + 1,:) = CNSim(CrInd + i - 1);
    Cycles.IRP(i + window + 1,:) = PSim(CrInd + i);
    Cycles.IRZ(i + window + 1,:) = yt(Posterioryt(CrInd + i));
    Cycles.IRZN(i + window + 1,:) = yn(Posterioryn(CrInd + i));
    Cycles.IRgrYTSim(i + window + 1,:) = grYTSim(CrInd + i) + g;
    Cycles.IRgrYNSim(i + window + 1,:) = grYNSim(CrInd + i) + g;
    Cycles.IRG(i + window + 1,:) = gt(Posteriorg(CrInd + i)) + g;       
    Cycles.IRZSim(i + window + 1,:) = Simyt(CrInd + i);
    Cycles.IRZNSim(i + window + 1,:) = Simyn(CrInd + i);
    Cycles.IRGSim(i + window + 1,:) = Simg(CrInd + i) + g;
    Cycles.IRDtoY(i + window + 1,:) = DtoY(CrInd + i); 
    Cycles.IRTAO(i + window + 1,:) = TAOSim(CrInd + i);
    Cycles.IRCtoY(i + window + 1,:) = CtoY(CrInd + i + 1);
    Cycles.IRCTtoY(i + window + 1,:) = CTtoY(CrInd + i);
    Cycles.IRCNtoY(i + window + 1,:) = CNtoY(CrInd + i);
    Cycles.IRCAtoY(i + window + 1,:) = CAtoY(CrInd + i);
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
Cycles.IRZSimMean = mean(Cycles.IRZSim, 2);
Cycles.IRZNSimMean = mean(Cycles.IRZNSim, 2);
Cycles.IRGSimMean = mean(Cycles.IRGSim, 2);
Cycles.IRDtoYMean = mean(Cycles.IRDtoY, 2);
Cycles.IRCtoYMean = mean(Cycles.IRCtoY, 2);
Cycles.IRCTtoYMean = mean(Cycles.IRCTtoY, 2);
Cycles.IRCTtoYMean = mean(Cycles.IRCTtoY,2);
Cycles.IRCNtoYMean = mean(Cycles.IRCNtoY,2);
Cycles.IRCAtoYMean = mean(Cycles.IRCAtoY,2);
Cycles.IRgrYTSimMean = mean(Cycles.IRgrYTSim, 2);
Cycles.IRgrYNSimMean = mean(Cycles.IRgrYNSim, 2);
Cycles.IRTAOMean = mean(Cycles.IRTAO, 2);
Cycles.IRYtotMean = mean(Cycles.IRYtot,2);


fprintf("Done \n")
  

%% Save Structure

fprintf("Saving ... \n")

IIPCC.b = b;
IIPCC.S = S;
IIPCC.V = V;
IIPCC.yt = yt;
IIPCC.yn = yn;
IIPCC.gt = gt;
IIPCC.gT = gT;
IIPCC.gN = gN;
IIPCC.Pi = Pi;
IIPCC.Ytot = Ytot;
IIPCC.Sim = Sim;
IIPCC.Pol = Pol;
IIPCC.TAO = TAO;
IIPCC.TAOSim = TAOSim;
IIPCC.grYTSim = grYTSim;
IIPCC.grYNSim = grYNSim;
IIPCC.COpt = COpt;
IIPCC.POpt = POpt;
IIPCC.BCOpt = BCOpt;
IIPCC.CTOpt = CTOpt;
IIPCC.mCT = mCT;
IIPCC.mCN = mCN;
IIPCC.mP = mP;
IIPCC.CA = CA;
IIPCC.SimB = SimB;
IIPCC.Index = Index;
IIPCC.PSim = PSim;
IIPCC.CSim = CSim;
IIPCC.CTSim = CTSim;
IIPCC.BCSim = BCSim;
IIPCC.Cycles = Cycles;
IIPCC.AAA = AAA ;
IIPCC.CCC = CCC ;
IIPCC.CrInd = CrInd;
IIPCC.Freq = Crisis2 
IIPCC.DtoY = DtoY;
IIPCC.CtoY = CtoY;
IIPCC.CTtoY = CTtoY;
IIPCC.CNtoY = CNtoY;
IIPCC.CAtoY = CAtoY;
IIPCC.mDtoY = mean(DtoY)
IIPCC.mCtoY = mean(CtoY);
IIPCC.mCAtoY = mean(CAtoY);
IIPCC.mCTtoY = mean(CTtoY);
IIPCC.mCNtoY = mean(CNtoY);
IIPCC.stdZBB = ZBBT;
IIPCC.stdGBB = GBBT;
IIPCC.Posterioryt = Posterioryt;
IIPCC.Posterioryn = Posterioryn;
IIPCC.Posteriorg = Posteriorg;

save('Data/IIPCC.mat', 'IIPCC', '-v7.3')
save("Data/IIPCCsim.mat",'-struct','IIPCC', '-v7.3') ;





