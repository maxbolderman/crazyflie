%% Acc2_calc_Output
% This model calculates the output, based on the state.
% Two different output equations are presented here;
%   - direct position measurement
%   - distance measurement

function y = Acc2_calc_Output(x, method)
y = zeros(3,1);
if (method == 0)
    %% Direct position measurement
    y = x(1:3);
elseif (method == 1)
    %% Distance measurement
    anchor = [0, 0, 0;             % Positions of the UWB sensors
              0, 4, 0;
              0, 0, 3];
    y(1) = sqrt((anchor(:,1)-x(1:3))'*(anchor(:,1)-x(1:3)));
    y(2) = sqrt((anchor(:,2)-x(1:3))'*(anchor(:,2)-x(1:3)));
    y(3) = sqrt((anchor(:,3)-x(1:3))'*(anchor(:,3)-x(1:3)));
end
    



end