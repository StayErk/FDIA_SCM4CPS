function [filtered_sig, residue, cov_matrix, p0, K, state_error, x_hat]= kalmanfilter(signal,Q_in,R_in,simulation_time,x)
    %this Kalman filter has been modified to follow an input siusoidal
    %"signal" with approximatley 60 Hz.
    %Q_in and R_in are tuning parameters.
    %simulation_time is the total simulation time(match input-signal).
    %filtered_sig is the Kalman-filtered output of the signal.
    time_vec=linspace(0,simulation_time,length(signal));
    A=[1.00 0.0;0.0 1.00];
    w=2*pi*60;
    %initilize values
    x0=[0;0];   %assume initial state = 0;
    p0=1*eye(size(A));
    Q=Q_in;
    R=R_in;
    K =0;
    state_error1 = zeros(length(signal));
    state_error2 = zeros(length(signal));
    state_error = zeros(length(signal), 2);
    z1 = zeros(size(signal));
    z2 = zeros(size(signal));
    residue = zeros(size(signal));
    cov_matrix = zeros(200, 2);
    for i=1:max(size(signal))
        % x_hat(k|k-1) = A x_hat(k-1)
        x_hat_pred = A*x0;

        C = [cos(w*time_vec(i)) -sin(w*time_vec(i))];
        
        %1 prediction error covariance P(k|k-1) = A P(k-1) A'+ Q
        pre_error_cov = A*p0*A'+Q;
        
        %2 kalman gain equation k(k) = P(k|k-1) C' (C P(k|k-1) C' + R)^-1
        K = A*pre_error_cov*C'*inv(C*pre_error_cov*C'+R);
        
        %3 filter equation x_hat(k) = x_hat(k|k-1) + k(k)( y(k) - C
        %x_hat(k|k-1) )
        x_hat = x_hat_pred + K*(signal(i)-C*x_hat_pred);
        
        %4 error covariance P(k) = (I - k(k) C) P(k|k-1)
        error_cov = pre_error_cov-K*C*pre_error_cov;
        
        z1(i)=x_hat(1);
        z2(i)=x_hat(2);

        %recreating the sinusoid   y_hat = C x_hat
        filtered_sig(i)=C*[z1(i);z2(i)];

        % residuo: r(t) = y(t)-y_hat(t|t-1)
        residue(i) = signal(i)-C*x_hat_pred;
        cov_matrix(i,1) = signal(i);
        cov_matrix(i,2) = C*x_hat_pred;
        state_error1(i) = x(1) - z1(i);
        state_error2(i) = x(2)- z2(i);
        x0=x_hat;
        p0=error_cov;   
        K = K;
        
    end
    state_error(:,1) = state_error1;
    state_error(:,2) = state_error2;
    x_hat = zeros(200, 2);
    x_hat(:,1) = z1; 
    x_hat(:,2) = z2;
end