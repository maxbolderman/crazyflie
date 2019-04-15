%% Acc2_main_DTUKF
% This model estimates the state of the Quadcopter. The state is given as 
% x = [x,y,z, xdot,ydot,zdot, theta,phi,psi]
% The inputs used for the model are:    
%   - Acceleration measurements from IMU
%   - Gyroscopic measurements from IMU
% The output measurement are:
%   - Distance measurements with respect to given anchor positions

clear all;
close all;
clc;

%% PARAMETERS
% Standard parameters
Acc2_parameters;
resetP = 0;
% Noise parameters
muw = zeros(1,n);
muv = zeros(1,nout);
noiseAcc = 0.5;
Q = diag([0.5*Ts*Ts*noiseAcc,0.5*Ts*Ts*noiseAcc,0.5*Ts*Ts*noiseAcc, Ts*noiseAcc,Ts*noiseAcc,Ts*noiseAcc, Ts*0.1,Ts*0.1,Ts*0.1]).^2;
R = diag([0.25,0.25,0.25]).^2;
% Input parameters
freq = 1*2*pi;
angle = [pi/4*sin(freq*t);
         pi/4*sin(freq*t + 1/3*pi);
         pi/4*sin(freq*t + 2/3*pi)];
% Initialization parameters
P0 = diag([1,1,1, 0.01,0.01,0.01, 0.01,0.01,0.01]).^2;
P = P0;
x0 = [zeros(n-3,1); angle(:,1)];
initialerror = [0.1,0.1,0.1, 0.1,0.1,0.1, 0.1,0.1,0.1]';
x0hat = x0 + initialerror;
% Efficiency variables
traceP = zeros(1,size(t,2));
x = zeros(n, size(t,2));
xhat = zeros(n, size(t,2));
xhatpred = zeros(n, size(t,2));
y = zeros(nout, size(t,2));
yreal = zeros(nout, size(t,2));
yhat = zeros(nout, size(t,2));
yhatpred = zeros(nout, size(t,2));
traceP(1) = trace(P);
x(:,1) = x0;
xhat(:,1) = x0hat;

%% Keep track of maximum values
maxX = zeros(size(t));
maxY = zeros(size(t));
maxPxx = zeros(size(t));
maxPyy = zeros(size(t));
maxPxy = zeros(size(t));
maxPyyinv = zeros(size(t));
maxSpread = zeros(size(t));
maxK = zeros(size(t));
%% SIMULATION
% Running two besides each other: 1 for the simulation, 1 to estimate the
% simulation
for ii = 1:1:size(t,2)-1
    %% SIMULATE
    %[x(k+1), w, a] = Acc2_sim(simornot, time, acceleration, angle, omega,
    % state, timestep)
    [x(:,ii+1), omega, acceleration] = Acc2_sim(1, t(ii), [0;0;0], angle(:,ii+1), [0;0;0], x(:,ii), Ts);
    x(:,ii+1) = x(:,ii+1) + mvnrnd(muw,Q)';
    yreal(:,ii+1) = Acc2_calc_Output(x(:,ii+1), method);
    y(:,ii+1) = yreal(:,ii+1) + mvnrnd(muv,R)';
    
    %% PREDICT
    % define sigma points
    [status, spread] = choleskyGreiff(P);
    if (status == 0)
        P = P0;
        resetP = resetP+1;
        [status, spread] = choleskyGreiff(P);
    end
    sigmaX = [xhat(:,ii), xhat(:,ii)+sqrt(n+lambda)*transpose(spread), xhat(:,ii)-sqrt(n+lambda)*transpose(spread)]; 
    % integrate timestep
    sigmaXplus = zeros(size(sigmaX));
    sigmaYplus = zeros(nout, nsigma);
    for jj = 1:1:nsigma
        [sigmaXplus(:,jj), ~, ~] = Acc2_sim(0, 0, acceleration, [0;0;0], omega, sigmaX(:,jj), Ts);
        sigmaYplus(:,jj) = Acc2_calc_Output(sigmaXplus(:,jj), method);
    end
    xhatpred(:,ii+1) = sigmaXplus*Wm';
    yhatpred(:,ii+1) = sigmaYplus*Wm';
    
    Pxx = Q;
    Pyy = R; 
    Pxy = zeros(n, nout);
    for jj = 1:1:nsigma
        Pxx = Pxx + Wc(jj) * (sigmaXplus(:,jj)-xhatpred(:,ii+1))*(sigmaXplus(:,jj)-xhatpred(:,ii+1))';
        Pyy = Pyy + Wc(jj) * (sigmaYplus(:,jj)-yhatpred(:,ii+1))*(sigmaYplus(:,jj)-yhatpred(:,ii+1))';
        Pxy = Pxy + Wc(jj) * (sigmaXplus(:,jj)-xhatpred(:,ii+1))*(sigmaYplus(:,jj)-yhatpred(:,ii+1))';
    end   
    
    %% UPDATE
    K = Pxy / Pyy;
    
    if (mod(ii,samples) == 0)
        % Perform update with position measurement
        xhat(:,ii+1) = xhatpred(:,ii+1) + K*(y(:,ii+1)-yhatpred(:,ii+1));
        P = Pxx - K*Pyy*K';
    else 
        % Set prediction to estimation -> no position measurement available
        xhat(:,ii+1) = xhatpred(:,ii+1);
        P = Pxx;
    end
    if (xhat(7,ii+1) > pi) 
        xhat(7,ii+1) = x(7,ii+1)-2*pi;
    elseif (xhat(7,ii+1) <-pi)
        xhat(7,ii+1) = x(7,ii+1)+2*pi;
    end
    if (xhat(8,ii+1) > pi) 
        xhat(8,ii+1) = x(8,ii+1)-2*pi;
    elseif (xhat(8,ii+1) <-pi)
        xhat(8,ii+1) = x(8,ii+1)+2*pi;
    end
    if (xhat(9,ii+1) > pi) 
        xhat(9,ii+1) = x(9,ii+1)-2*pi;
    elseif (xhat(9,ii+1) <-pi)
        xhat(9,ii+1) = x(9,ii+1)+2*pi;
    end
    
    yhat(:,ii+1) = Acc2_calc_Output(xhat(:,ii+1),method);
    traceP(ii+1) = trace(P);
    if (traceP(ii+1) > 1000)
        P = P0;
        traceP(ii+1) = trace(P);
        resetP = resetP+1;
    end
    % Set the maximum values
    maxX(ii) = max(abs(xhat(:,ii)));
    maxY(ii) = max(max(abs(sigmaYplus)));
    maxPxx(ii) = max(max(abs(P)));
    maxPyy(ii) = max(max(abs(Pyy)));
    maxK(ii) = max(max(abs(K)));
    maxPxy(ii) = max(max(abs(Pxy)));
    Pyyinv = inv(Pyy);
    maxPyyinv(ii) = max(max(abs(Pyyinv)));
    maxSpread(ii) = max(max(abs(spread)));
    
end

%% RESULTS
% Compute RMS of the errors
RMS = rms(xhat-x,2)

% Here all necessary plots are presented
% Figure(1); measurement error and estimated output error
figure(1) 
subplot(3,1,1)
plot(t,y(1,:)-yreal(1,:))
hold on 
grid on
plot(t,yhat(1,:)-yreal(1,:))
ylabel('Anchor 1 [m]')
subplot(3,1,2)
plot(t,y(2,:)-yreal(2,:))
hold on
grid on
plot(t,yhat(2,:)-yreal(2,:))
ylabel('Anchor 2 [m]')
subplot(3,1,3)
plot(t,y(3,:)-yreal(3,:))
hold on
grid on
plot(t,yhat(3,:)-yreal(3,:))
ylabel('Anchor 3 [m]')
xlabel('Time [s]')
legend({'Measurement', 'Estimation'}, 'FontSize', 14)

% Figure(2); errors of the position
figure(2)
plot(t, xhat(1,:)-x(1,:))
hold on
grid on
plot(t, xhat(2,:)-x(2,:))
plot(t, xhat(3,:)-x(3,:))
xlabel('Time [s]')
ylabel('Error [m]')
title('Position error')
legend({'x', 'y', 'z'}, 'FontSize', 14)

% Figure(3); errors of the velocity
figure(3)
plot(t, xhat(4,:)-x(4,:))
hold on 
grid on
plot(t, xhat(5,:)-x(5,:))
plot(t, xhat(6,:)-x(6,:))
xlabel('Time [s]')
ylabel('Error [m/s]')
title('Velocity error')
legend({'x', 'y', 'z'}, 'FontSize', 14)

% Figure(4); errors of the attitude
figure(4)
plot(t, xhat(7,:)-x(7,:))
hold on
grid on
plot(t, xhat(8,:)-x(8,:))
plot(t, xhat(9,:)-x(9,:))
xlabel('Time [s]')
ylabel('Error [rad]')
title('Attitude error')
legend({'\theta', '\phi', '\psi'}, 'FontSize', 14)

% Figure(5); showing real and estimated angle
figure(5) 
subplot(3,1,1)
plot(t, xhat(7,:))
hold on 
grid on
plot(t, x(7,:))
ylabel('\theta [rad]')
subplot(3,1,2)
plot(t, xhat(8,:))
hold on 
grid on
plot(t, x(8,:))
ylabel('\phi [rad]')
subplot(3,1,3)
plot(t, xhat(9,:))
hold on 
grid on
plot(t, x(9,:))
ylabel('\psi [rad]')
xlabel('Time [s]')
legend({'Estimated', 'Real'}, 'FontSize', 14)
figure(5)
title('Angle and estimated angle')

%Figure(6); trace of the covariance matrix P
figure(6)
subplot(2,1,1)
plot(t, traceP)
grid on
xlabel('Time [s]')
ylabel('Trace P')
subplot(2,1,2)
plot(traceP)
grid on
xlabel('Iteration [-]')
ylabel('Trace P')

%Figure(7); plot all the maximum values
figure(7)
plot(t, maxX)
hold on 
grid on
plot(t, maxY)
plot(t, maxPxx)
plot(t, maxPyy)
plot(t, maxK)
plot(t, maxPxy)
plot(t, maxPyyinv)
plot(t, maxSpread)
legend({'maxX', 'maxY', 'maxPxx', 'maxPyy', 'maxK', 'maxPxy', 'maxPyyinv', 'maxSpread'}, 'FontSize', 14)









