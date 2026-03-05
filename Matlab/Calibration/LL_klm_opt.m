function L = LL_klm_opt(delta, varargin)


Z = varargin{1};
T = size(Z, 2);
N=4;

% 1. State-Space Form

A = [delta(1, 1)  delta(1, 2)   0  0;
       delta(1, 3)  delta(1, 4)   0  0;
       0               0            delta(1, 5)  0;
       1 0 0 0];
 
C = [delta(1,6)   delta(1, 7)       0               0;
       delta(1, 7)   delta(1, 8)        0              0;
       0                0               delta(1, 9)       0;
       0                0                      0              0];

D = [1 0 1 -1;
        0 1 1 -1];
 
E = [0; 0; 0; 0];

% Stationarity restrictions

% 2. Initialization

    n = size(A, 1);
    X = zeros(n, T);
    P_tt = 1e3*eye(n);

    % Likelihood Evaluation

    L = 0;

    for t = 1:T
        if t == 1
            X_t1 = zeros(n,1); 
        else 
            X_t1 = A*X(:,t-1) + E; 
        end

        P_tt1 = A*P_tt*A' + C*C';
        
        Omega = D*P_tt1*D';
        Omegainv = eye(size(Z,1))/Omega;
        Kt = P_tt1*D'/Omega;
        Ztilde = Z(:, t) - D*X_t1;
        X(:, t) = X_t1 + Kt*Ztilde;
        P_tt = P_tt1 - Kt*Omega*Kt';
        
        
           L = L - 0.5*(T*log(det(Omega))) - 0.5*Ztilde'*Omegainv*Ztilde;
          
    end

        L = L - 0.5*T*N*log(2*pi);
end





