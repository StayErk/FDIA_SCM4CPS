% Prova Kalman Filter 

% Modello
% x[n + 1] = Ax[n] + Bu[n] + Gw[n] (in the example B = G)
%     y[n] = Cx[n] + Du[n] + Hw[n] + v[n]

A = [1.1269, -0.4940, 0.1129; 1, 0, 0; 0, 1, 0];

B = [-0.3832; 0.5919; 0.5191];

C = [1,0,0];

D = 0;

% Set up plant model:

Ts = -1; % discrete time
plant = ss(A, [B B], C, D, Ts);

process_noise_cov = 2.3;
sensor_noise_cov  = 1;

% Set up the kalman filter:

[kalman_filter, L, ~, Mx, Z] = kalman(plant, process_noise_cov, sensor_noise_cov);

% Update pant I/O and Kalman I/O

plant.InputName = {'u', 'w'};
plant.OutputName = {'yt'}; % true output of the plant

measurementNoiseAdd = sumblk('y = yt + v'); % input of kalman filter y is the true outuput of the plant plus measurement noise

kalman_filter.InputName  = {'u', 'y'};
kalman_filter.OutputName = {'y_hat'};

simulation = connect(plant, measurementNoiseAdd, kalman_filter, {'u', 'w', 'v'}, {'yt', 'y_hat'});

% To simulate the filter behavior we have to generate an input.
% use a sinusoidal input vector

t = (0:100)';
u = sin(t/5);

rng(10, 'twister');
inputNoise       = sqrt(process_noise_cov) * randn(length(t), 1);
measurementNoise = sqrt(sensor_noise_cov) * randn(length(t),1);

response = lsim(simulation, [u, inputNoise, measurementNoise]);

yt    = response(:,1);
y_hat = response(:,2);
y     = yt + measurementNoise; 

clf
subplot(211), plot(t,yt,'b',t,y_hat,'r--'), 
xlabel('Number of Samples'), ylabel('Output')
title('Kalman Filter Response')
legend('True','Filtered')
subplot(212), plot(t,yt-y,'g',t,yt-y_hat,'r--'),
xlabel('Number of Samples'), ylabel('Error')
legend('True - measured','True - filtered')