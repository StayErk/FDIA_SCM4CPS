function [filtered_sig, residue, cov_matrix, p0, K,  x_hat,P1]= kalmanfilter(y,Q_in,R_in,t)

    time_vec=linspace(0,t,length(y));
    A=[1 0;0,1];
    w=2*pi*60;
    %initilize values
    x0=[0;0];   %assume initial state = 0;
    p0=1*eye(size(A));
    Q=Q_in;
    R=R_in;
    K =0;
    C =  [cos(omega*t), -sin(omega*t)];

    z1 = zeros(size(y));
    z2 = zeros(size(y));
    residue = zeros(size(y));
    cov_matrix = zeros(200, 2);
    for i=1:max(size(y))
        % x_hat(k|k-1) = A x_hat(k-1)
        x_hat_pred = A*x0;

        %1 prediction error covariance P(k|k-1) = A P(k-1) A'+ Q
        pre_P = A*p0*A'+Q;
        
        %2 kalman gain equation k(k) = P(k|k-1) C' (C P(k|k-1) C' + R)^-1
        K = A*pre_P*C'*inv(C*pre_P*C'+R);
        
        %3 filter equation x_hat(k) = x_hat(k|k-1) + k(k)( y(k) - C
        %x_hat(k|k-1) )
        x_hat = x_hat_pred + K*(y(i)-C*x_hat_pred);
        
        %4 error covariance P(k) = (I - k(k) C) P(k|k-1)
        P = pre_P-K*C*pre_P;
        
        z1(i)=x_hat(1);
        z2(i)=x_hat(2);

        %recreating the sinusoid   y_hat = C x_hat
        filtered_sig(i)=C*[z1(i);z2(i)];

        % residuo: r(t) = y(t)-y_hat(t|t-1)
        residue(i) = y(i)-C*x_hat_pred;
        cov_matrix(i,1) = y(i);
        cov_matrix(i,2) = C*x_hat_pred;
        x0=x_hat;
        p0=P;   
        P1(:,:,i)=P;   
        K = K;
        
    end
