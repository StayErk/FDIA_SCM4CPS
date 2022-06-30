define_constants;

% matlab colors and preset
NMatlabRed      = [0.8500   0.3250   0.0980];
NMatlabYellow   = [0.929    0.694    0.125 ];
NMatlabBlue     = [0        0.4470   0.7410];  
NMatlabViolet   = [0.4940   0.1840   0.5560];
NMatlabGreen    = [0.4660   0.6740   0.1880];
NMatlabCyan     = [0.3010   0.7450   0.9330];
NMatlabBordeaux = [0.6350   0.0780   0.1840]; 

s = tf('s');
om = logspace(-2,5,1000);
t = 0:0.001:10;
% ============================================= %


mpc = loadcase('case9');
results = runopf(mpc);
bus3_voltageMagnitude = results.bus(3, VM);
bus3_voltageAngles = results.bus(3, VA);
input_noise = randn(2);
value_noise = randn(2);
Ts = -1;
A = eye(2);
C = [cos(60), -sin(60)];
Q = 2.3;
R = 1;


plant = ss(A, [0;1], C, 0, Ts, 'InputName', {''}, 'OutputName', 'y');
plant.InputName = {'u', 'w'};
plant.OutputName = {'yt'};

noiseBlock = sumblk('y = yt + v');

[kalman_filter, filter_gain, prediction_err_cov, ~, err_cov] = kalman(sys, Q, R);

kalman_filter.InputName = {'u', 'y'};
kalman_filter.OutputName = {'y_hat'};

Model = connect(plant, noiseBlock, kalman_filter, {'u', 'w', 'v'}, {'yt', 'y_hat'});








