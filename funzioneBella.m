function [filtered_sig, residue]= kalmanfilter(signal,Q_in,R_in,simulation_time)
    %this Kalman filter has been modified to follow an input siusoidal
    %"signal" with approximatley 50 Hz.
    %Q_in and R_in are tuning parameters.
    %simulation_time is the total simulation time(match input-signal).
    %filtered_sig is the Kalman-filtered output of the signal.
    time_vec=linspace(0,simulation_time,length(signal));
    A=[1.00 0.0;0.0 1.00];
    w=2*pi*60;
    %initilize values
    x0=[0;0];
    p0=1*eye(2);
    Q=Q_in;
    R=R_in;
    z1 = zeros(size(signal));
    z2 = zeros(size(signal));
    residue = zeros(size(signal));
    
    filtered_sig = zeros(size(signal));
    for i=1:max(size(signal))
        
        xn_nm1=x0;
        
        c=[cos(w*time_vec(i)) -sin(w*time_vec(i))];
        
        %1 Corresponding to flowchart chapter 3
        pn_nm1=A*p0*A'+Q;
        
        %2 Corresponding to flowchart 2
        kn=A*pn_nm1*c'*inv(c*pn_nm1*c'+R);
        
        %3 Corresponding to flowchart chapter 3
        xnn=xn_nm1+kn*(signal(i)-c*xn_nm1);
        
        %4 Corresponding to flowchart chapter 3
        pnn=pn_nm1-kn*c*pn_nm1;
        
        z1(i)=xnn(1);
        z2(i)=xnn(2);

        %recreating the sinusoid
        filtered_sig(i)=c*[z1(i);z2(i)];
        residue(i) = signal(i)-c*xn_nm1;
        
        x0=xnn;
        p0=pnn;
        
    end
end