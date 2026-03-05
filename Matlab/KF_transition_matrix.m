%% Description

function [Tran, Indd]  = KF_transition_matrix(Sim,  yt, yn, gt, gT, gN, T, g, cc, k1, k2, ynn, ytn, gn, gTn, gNn) % This is the transition matrix given by the agent's beliefs

Simyt = Sim(:,1) ; Simyn = Sim(:,2) ; Simg = Sim(:,3) ; 
SimgT = Sim(:,4) ; SimgN = Sim(:,5) ; 

NumSim = length(Sim) ;

rng(5); %Seed to keep the shocks constant in different simulations.

IndApp=[findClosest2(Simyt,yt) findClosest2(Simyn,yn) findClosest2(Simg,gt) findClosest2(SimgT,gT) findClosest2(SimgN,gN)]; % Assign those observables to the closest points in the grid
States=1:1:(ytn*ynn*gn*gTn*gNn); %a count on the total number of exogenous states (Ny*Nr)
States=reshape(States,ytn,ynn,gn,gTn,gNn);
Posterioryt = IndApp(1,1) ; % Posterior is the correct realization of the transition matrix
Posterioryn = IndApp(1,2) ; % Posterior is the correct realization of the transition matrix
Posteriorg  = IndApp(1,3) ; % Posterior is the correct realization of the transition matrix
GaT         = IndApp(1,4) ; 
GaN         = IndApp(1,5) ;

ythat = zeros(1,NumSim-1);
ynhat = zeros(1,NumSim-1);
ghat  = zeros(1,NumSim-1);
I     = zeros(1,NumSim-1);
j     = zeros(1,NumSim-1);
x     = zeros(1,NumSim-1);

a = [Sim(1,1); Sim(1,2); g+Sim(1,3); (Sim(1,1) + Sim(1,3) - Sim(1,4)+g)]; % Posterior for the first period is a (which is correct)

for i=2:NumSim-1
    att1 = T*a(:,i-1) + cc ; % Prior is completely correct whis is the observed information in the previous period
    a(:,i) = k1*att1 + k2*[(SimgT(i)+g);(SimgN(i)+g)] ; % Posterior is formed using the kalman weights
    ythat(i) = a(1,i) ;
    ynhat(i) = a(2,i) ; % Assigning the values of the kalman equations
    ghat(i)  = a(3,i) ;
    State = States(Posterioryt, Posterioryn, Posteriorg, GaT, GaN) ;
    GaT = IndApp(i,4) ;
    GaN = IndApp(i,5) ;
    Posterioryt = findClosest2(ythat(i),yt) ; % Approximate the posterior to the points in the grids
    Posterioryn = findClosest2(ynhat(i),yn) ;
    Posteriorg  = findClosest2(ghat(i)-g,gt) ;
    I(i) = State;
    j(i) = States(Posterioryt,Posterioryn,Posteriorg,GaT,GaN);
    x(i) = 1; %Count each transition from state i to state j.
end

Conteo=sparse(I(2:end),j(2:end),x(2:end),ytn*ynn*gn*gTn*gNn,ytn*ynn*gn*gTn*gNn);
SumConteo=sum(Conteo,2);
Indd=find(SumConteo>=1);
Conteo2=Conteo(any(Conteo,2),:);
Conteo3=Conteo2(:,any(Conteo,1));
Conteo = full(Conteo3);
SumConteo=sum(Conteo,2);
Tran=Conteo./repmat(SumConteo,[1 length(Conteo)]); % Generate transition probability function including the bayesian update

end