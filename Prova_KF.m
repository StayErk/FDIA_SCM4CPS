% Prova Kalman Filter 

% Modello
% x[n + 1] = Ax[n] + Bu[n] + Gw[n] (in the example B = G)
%     y[n] = Cx[n] + Du[n] + Hw[n] + v[n]

A = eye(2); % [1.1269, -0.4940, 0.1129; 1, 0, 0; 0, 1, 0];

B = [1;1];

C = [cos(2*pi*Fc*t), -sin(2*pi*Fc*t)];

D = zeros(size(C,1),size(B,2));

% Set up plant model:

Ts = -1; % discrete time
plant = ss(A, B, C, D, Ts);

process_noise_cov = 2.3;
sensor_noise_cov  = 1;

% Set up the kalman filter:

[kalman_filter, L, ~, Mx, Z] = kalman(plant, process_noise_cov, sensor_noise_cov);

% Update pant I/O and Kalman I/O



plant.InputName = {'u'};
plant.OutputName = {'yt'}; % true output of the plant

measurementNoiseAdd = sumblk('y = yt + v'); % input of kalman filter y is the true outuput of the plant plus measurement noise
inputNoiseAdd = sumblk('u = ut + w');

kalman_filter.InputName  = {'ut', 'y'};
kalman_filter.OutputName = {'y_hat', 'x1_hat', 'x2_hat', 'x3_hat'};

simulation = connect(inputNoiseAdd, plant, measurementNoiseAdd, kalman_filter, {'ut', 'w', 'v'}, {'yt', 'y_hat', 'x1_hat', 'x2_hat', 'x3_hat'});

% To simulate the filter behavior we have to generate an input.
% use a sinusoidal input vector

t = (0:100)';
ut = sin(t/5);

rng(10, 'twister');
inputNoise       = sqrt(process_noise_cov) * randn(length(t), 1);
measurementNoise = sqrt(sensor_noise_cov) * randn(length(t),1);

[response, ~, state] = lsim(simulation, [ut, inputNoise, measurementNoise]);

yt    = response(:,1);
y_hat = response(:,2);
x1_hat = response(:,3);
x2_hat = response(:,4);
x3_hat = response(:,5);
y     = yt + measurementNoise; 

clf
subplot(311), plot(t,yt,'b',t,y_hat,'r--'), 
xlabel('Number of Samples'), ylabel('Output')
title('Kalman Filter Response')
legend('True','Filtered')
subplot(312), plot(t,yt-y,'g',t,yt-y_hat,'r--'),
xlabel('Number of Samples'), ylabel('Error')
legend('True - measured','True - filtered')

subplot(313), plot(t, state,'g',t,[x1_hat, x2_hat,x3_hat],'r--'),
xlabel('Number of Samples'), ylabel('State')
legend('True - state','True - filtered')