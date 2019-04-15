%% Acc_parameters
% In this model, all the variables are declared

% measurement 
method = 1;     % 0 = direct position, 1 = distance 
samples = 10;    % Amount of samples before three measurements (chose 1 to have each sample)

% Sampling
Fs = 100;
Ts = 1/Fs;
Ttotal = 30;
t = 0:Ts:Ttotal;

% Constants
g = 9.81;
l = 0.046;
m = 0.027;

% UKF
n = 9;
nsigma = 2*n+1;
nout = 3;
alpha = 0.1;
kappa = 0;
beta = 2;
lambda = alpha*alpha*(n+kappa) - n;
Wm = zeros(1,nsigma);
Wc = zeros(1,nsigma);
Wm(1) = lambda/(n+lambda);
Wm(2:end) = 1/(2*(n+lambda));
Wc(1) = lambda/(n+lambda) + 1-alpha^2+beta;
Wc(2:end) = 1/(2*(n+lambda));