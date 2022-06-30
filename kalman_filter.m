function [ x, y, K_m, P_m ] = kalman_filter( sysd, u, w, z, t, x_0, P_0, R, Q )

% discrete-time model
F = sysd.a;
G = sysd.b(:,1);
H = sysd.c;
J = sysd.b(:,2:end);

% initialization
k = 1 +1;
x = x_0;                 
P = P_0;
I = diag( ones(length(F),1) );
K_m = []; P_m = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% kalman filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:length(t)-1

    % 1. prediction
    x(:,k) = F*x(:,k-1) +G*u(k-1) +J*w(:,k-1);  % 'a priori' predicted state
    P      = F*P*F' +J*Q*J';                    % 'a priori' predicted covariance

    % 2. correction
    K      = P*H'*inv( H*P*H' +R );             % computes optimal gain   
    e(k)   = z(k) -H*x(:,k);                    % measurement residual
    x(:,k) = x(:,k) +K*e(k);                    % 'a posteriori' state estimative 
    P      = (I -K*H)*P;                        % 'a posteriori' covariance
    
    y(k)   = H*x(:,k);                          % estimated output
    k      = k +1;                              
    
    K_m    = [K_m K];                           % stores current K
    P_m    = [P_m P];                           % stores current P            
end

x = transpose(x);
end