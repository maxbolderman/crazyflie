%% Acc2_sim
% This model simulates the system. The same system is used for the
% simulation as well as the estimation.

    %[x(k+1), w, a] = Acc2_sim(simornot, time, acceleration, angle, omega,

function [dx, omega, acceleration] = Acc2_sim(sim, t, acceleration, angle, gyro, x, Ts)
g = 9.81;
dx = zeros(size(x));
if (sim == 1)
    %% Perform the simulation
    position = [sin(pi*(t+Ts));
                sin(pi*(t+Ts)+2/3*pi);
                sin(pi*(t+Ts)+2/3*pi)];
    velocity = [pi*cos(pi*(t+Ts));
                pi*cos(pi*(t+Ts)+2/3*pi);
                pi*cos(pi*(t+Ts)+2/3*pi)];
    acceleration = [-pi*pi*sin(pi*(t+Ts));
                    -pi*pi*sin(pi*(t+Ts)+2/3*pi);
                    -pi*pi*sin(pi*(t+Ts)+2/3*pi)];
%     position = [2; 5; 0.5];
%     velocity = [0;0;0];
%     acceleration = [0;0;0];
%     angle = [0;0;0];
    dx(1:3) = position;
    dx(4:6) = velocity;
    dx(7:9) = angle;
    % calculate the measurement of acceleration and gyroscope
    Rgb = calc_Rgb(x(7:9));
    W = calc_W(x(7:9),1);
    acceleration = Rgb*(acceleration + [0;0;g]);
    etadot = (dx(7:9)-x(7:9))/Ts;
    omega = W*etadot;
    
elseif (sim == 0)
    %% Perform the estimation
    Rgb = calc_Rgb(x(7:9));
    Rbg = transpose(Rgb);
    %W = calc_W(x(7:9),1);
    %Winv = inv(W);
    Winv = calc_W(x(7:9),-1);
    accglob = Rbg*acceleration - [0;0;g];
    dx(1:3) = x(1:3) + Ts*x(4:6) + 0.5*Ts*Ts*accglob;
    dx(4:6) = x(4:6) + Ts*accglob;
    dx(7:9) = x(7:9) + Ts*Winv*gyro;
    % calculate measurement -> not used after
    omega = [0;0;0];
    acceleration = [0;0;0];
    
    % Set between boundaries
    while (abs(dx(7)) > pi)
        if (dx(7) > pi)
            dx(7) = dx(7) - 2*pi;
        elseif (dx(7) < -pi)
            dx(7) = dx(7) + 2*pi;
        end
    end
    while (abs(dx(8)) > pi)
        if (dx(8) > pi)
            dx(8) = dx(8) - 2*pi;
        elseif (dx(8) < -pi)
            dx(8) = dx(8) + 2*pi;
        end
    end
    while (abs(dx(9)) > pi)
        if (dx(9) > pi)
            dx(9) = dx(9) - 2*pi;
        elseif (dx(9) < -pi)
            dx(9) = dx(9) + 2*pi;
        end
    end
end
end
