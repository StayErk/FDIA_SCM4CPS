define_constants;

% matlab colors and preset
NMatlabRed      = [0.8500   0.3250   0.0980];
NMatlabYellow   = [0.929    0.694    0.125 ];
NMatlabBlue     = [0        0.4470   0.7410];  
NMatlabViolet   = [0.4940   0.1840   0.5560];
NMatlabGreen    = [0.4660   0.6740   0.1880];
NMatlabCyan     = [0.3010   0.7450   0.9330];
NMatlabBordeaux = [0.6350   0.0780   0.1840]; 

% ============================================= %

% ========= Table Paper ===============%
Fc = 60; %frequency 60Hz
Av =1;  % amplitude
dt = 0.0005; % 1/0.0005 = 2 kHz
StopTime = 0.1;             % seconds
t = (0:dt:StopTime-dt)';     % seconds

v = Av*cos(2*pi*Fc*t)*cos(0)-Av*sin(2*pi*Fc*t)*sin(0);
% Plot the signal versus time:
figure;
plot(t,v);
ylim([-2 2]);
xlabel('time (in seconds)');
title('Voltage equation');


% =========  Model ===============%
% x_dot=Ax+Bu+Gw
% y=Cx+Du+Hw+v

Ts = -1;  %Set the sample time to -1 to mark the plant as discrete
% Matrices
A = [1,0;0,1];
B = [0;0];
C = [cos(2*pi*Fc*t), -sin(2*pi*Fc*t)];
D = zeros(size(C,1),size(B,2));

plant = ss(A,B,C,D,Ts);  % Plant dynamics and additive input noise w
plant.InputName = 'un';
plant.OutputName = 'yt';

%sampling time
dt = 0.0005; % 1/0.0005 = 2 kHz

eig(A)
istable(sys)  %wtf?

%add disturbance
G = [0;0]; %matrix of w

Q = 1; % process noise covariance
R = 1; % measurement noise covariance
N = 0;  % w and v uncorrelated

%system with noise
Sum = sumblk('un = u + w');
sys = connect(plant,Sum,{'u','w'},'yt');
% sys_KF = ss(A,[B G],C,[D D], dt);

%Kalman Filter
[KF,KF_gain,P] = kalman(sys,Q,R,N);
%check the gain is right
p_P = (A*P*A'+Q)
K = p_P*C'*(inv(C*p_P*C' + R));
error = norm(abs(K-KF_gain));  % should be 0

%assess the stability
Acb = A-K*C;
eig(Acb);




%======== come era prima 

% Sum = sumblk('un = u + w');
% sys = connect(plant,Sum,{'u','w'},'yt');
% size(sys)
% 
% Q = 1; % process noise covariance
% R = 1; % measurement noise covariance
% N = 0;  %=uncorrelated
% [KF,KF_gain,P] = kalman(sys,Q,R,N)
% size(KF)
% KF = KF(1,:);
% 
% vIn = sumblk('y=yt+v');
% KF.InputName = {'u','y'};
% KF.OutputName = 'ye';
% 
% SimModel = connect(sys,vIn,KF,{'u','w','v'},{'yt','ye'});
% t = (0:100)';
% u = sin(t/5);
% 
% rng(10,'twister');
% w = sqrt(Q)*randn(length(t),1);
% v = sqrt(R)*randn(length(t),1);
% 
% out = lsim(SimModel,[u,w,v]);
% yt = out(:,1);   % true response
% ye = out(:,2);  % filtered response
% y = yt + v;     % measured response
% 
% %Compare the true response with the filtered response.
% clf
% subplot(211), plot(t,yt,'b',t,ye,'r--'), 
% xlabel('Number of Samples'), ylabel('Output')
% title('Kalman Filter Response')
% legend('True','Filtered')
% subplot(212), plot(t,yt-y,'g',t,yt-ye,'r--'),
% xlabel('Number of Samples'), ylabel('Error')
% legend('True - measured','True - filtered')
% 

