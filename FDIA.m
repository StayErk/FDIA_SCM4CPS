% Try out MATPOWER library
define_constants; % to use constants such as PD, PC, PG ecc..
% runpf('case9') % this runs a Newton power flow on the 9-bus system and print it.

% Following tutorial on Kalman Filter. Prendo la mano con matlab.

% State space matrix


A = [-0.71  0.06 -0.19 -0.17;
      0.06 -0.52 -0.03  0.30;
     -0.19 -0.03 -0.24 -0.02;
     -0.17  0.30 -0.02 -0.41];

B = [ 1.44  2.91   0;
     -1.97  0.83 -0.27;
     -0.20  1.39  1.10;
     -1.2   0    -0.28];

C = [ 0    -0.36 -1.58 0.28;
     -2.05  0     0.51 0.03];

D = zeros(2,3);


Q = 2.3;
R = 1;

[kalmf, L, P] = kalman(sys, Q, R);

sys = ss(A,B,C,D);
sys.InputName = {'u1','u2','w'};
sys.OutputName = {'y1','y2'};

noiseInsertion1 = sumblk('y=y1+y2+v');
SimModel = connect(sys, noiseInsertion1, kalmf, {'u1', 'u2', 'w', 'v'}, {'y1', 'y2', 'y1_e', 'y2_e', 'x1_e', 'x2_e'});


t = (0:100)';

u1 = sin(t/5);

u2 = sin(t/5);

rng(10,'twister');
w = sqrt(Q)*randn(length(t),1);
v = sqrt(R)*randn(length(t),1);

out = lsim(SimModel,[u1, u2, w,v]);
