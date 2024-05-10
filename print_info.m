function print_info(agents, agent, fieldmap, N_matrix)
%PRINT_INFO Summary of this function goes here
%   Detailed explanation goes here

% Assuming you have data for velocity and target speed
velocity = sqrt(agents(agent).vx^2 + agents(agent).vy^2);
target_speed = agents(agent).target_speed;

accelerating = [];
if velocity > target_speed
    accelerating = false;
else
    accelerating = true;
end


% Create a table
T = table({num2str(velocity)}, {num2str(target_speed)}, 'VariableNames', {'Velocity', 'Target_Speed'});
disp(T);

vx = agents(agent).vx;
vy = agents(agent).vy;
motion_angle = atan2(vy, vx);
motion_string = [];

if motion_angle >= 0 && motion_angle <= pi/8 || motion_angle < 0 && motion_angle >= -pi/8
    motion_string = "syd";
    x_plus = 1;
    y_plus = 0;
elseif motion_angle > pi/8 && motion_angle <= 3*pi/8
    motion_string =  "sydöst";
    x_plus = 1;
    y_plus = 1;
elseif motion_angle > 3*pi/8 && motion_angle <= 5*pi/8
    motion_string = "öst";
    x_plus = 0;
    y_plus = 1;
elseif motion_angle > 5*pi/8 && motion_angle <= 7*pi/8
    motion_string = "nordöst";
    x_plus = -1;
    y_plus = 1;
elseif motion_angle > 7*pi/8 && motion_angle <= pi || motion_angle < -7*pi/8 && motion_angle >= -pi
    motion_string = "norr";
    x_plus = -1;
    y_plus = 0;
elseif motion_angle < -pi/8 && motion_angle >= -3*pi/8
    motion_string = "sydväst";
    x_plus = 1;
    y_plus = -1;
elseif motion_angle < -3*pi/8 && motion_angle >= -5*pi/8
    motion_string = "väst";
    x_plus = 0;
    y_plus = -1;
elseif motion_angle < -5*pi/8 && motion_angle >= -7*pi/8
    motion_string = "nordväst";
    x_plus = -1;
    y_plus = -1;
end

change = fieldmap(ceil(agents(agent).x * N_matrix) + x_plus, ceil(agents(agent).y * N_matrix) + y_plus) - fieldmap(ceil(agents(agent).x*N_matrix), ceil(agents(agent).y*N_matrix));

% Create a table
T = table({motion_string}, {change}, {accelerating}, 'VariableNames', {'motion_string', 'change', 'Acceleration (1 true, 0 false)'});
disp(T);

% target point
x_index = ceil(agents(agent).x * N_matrix);
y_index = ceil(agents(agent).y * N_matrix);
tx_index = ceil(agents(agent).tx * N_matrix);
ty_index = ceil(agents(agent).ty * N_matrix);

% Create a table
T = table({ty_index - y_index}, {-(tx_index - x_index)}, 'VariableNames', {'x_corrected_index difference to target', 'y_corrected_index difference to target'});
disp(T);

right_or_left_angular_acc = "";
if sign(agents(agent).omega) > 0
    right_or_left_angular_acc = "left";
else
    right_or_left_angular_acc = "right";
end

% Create a table
T = table({right_or_left_angular_acc}, 'VariableNames', {'angular acceleration direction'});
disp(T);

% Create a table
T = table({ceil(agents(agent).py*N_matrix) - y_index}, {-(ceil(agents(agent).px*N_matrix) - x_index)}, {agents(agent).bresenham}, 'VariableNames', {'x_corrected_index difference', 'y_corrected_index difference', 'breshenham'});
disp(T);

% Create a table
T = table({x_index}, {y_index}, 'VariableNames', {'row position', 'column position'});
disp(T);

disp("----------------------------------------------------")

end
