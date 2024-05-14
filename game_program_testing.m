function [iterations_until_completion, reason_for_exit, pursuer_wall_deaths, evader_wall_deaths] = ...
    game_program(No_of_evaders, No_of_pursuers, obstacle_density, visualization_boolean, ...
    sighting_target_boolean, jesus_boolean, smart_target_boolean, iterations_ahead, increment_boolean, increment, M, game, max_iter, crystal_tp_ahead, pursuit_crystal, initial_matrix, scale_factor, decay_rate)
    
    disp("not testing")

    % initial_matrix = readmatrix('initial_matrix.csv');
    frames = {0}; % Initialize an empty array to store frames
    
    % Playing field matrix
    N_matrix = 302; % spikad
    remaining_blue_agents = No_of_evaders;

    % M = random_walk_map_generator(obstacle_density, N_matrix, 25);    % column width spikad
    % M(1, :) = 3*ones(1, 302);
    % M(302, :) = 3*ones(1, 302);
    % M(:, 1) = 3*ones(302, 1);
    % M(:, 302) = 3*ones(302, 1);
    % 
    % M(:, 151:152) = 3*ones(302, 2);
    % 
    % M(:, 100:200) = 3*ones(302, 101);
    % 
    % M(145:155, 2:301) = zeros(11, 300);
    %

    % M_v = -3*(readmatrix("binary_image_v.csv") - 1);
    % 
    % M = 3*ones(N_matrix);
    % M(2:end-1, 2:end-1) = M_v;
    % M(end+1, :) = ones(1, N_matrix-2);
    % M = [ones(1, N_matrix-2); M];
    % M(:, end) = ones(N_matrix, 1);
    % M = [M, ones(N_matrix, 1)];
    % 
    % M_a = M == 1;
    % M = ones(N_matrix);
    % M = M - M_a;

    base_M = M;

    color_M = zeros(N_matrix);
    % Number of agents and their information
    agents = struct('x', [], 'y', [], 'vx', [], 'vy', [], 'ax', [], 'ay', [], ...
        'theta', [], 'omega', [], 'alpha', [], 'color', [], 'tx', [], ...
        'ty', [], 'exploration_map', [], 'bresenham', false, ...
        'px', [], 'py', [], 'px_list', [], 'py_list', [], 'x_index_memo', [], ...
        'y_index_memo', [], 'target_speed', [], 'evader_index', [], 'previous_x', [], 'previous_y', [], 'change_point', [], 'x0_pred', [], 'y0_pred', [], ...
        'x_list', [], 'y_list', []); %, 'temp', []);
    
    M_new = M == 3;
    M_new = double(M_new);

    iterationer = 300;  % spikad
    fieldmap = getFieldmap(M_new, iterationer);
   
    if visualization_boolean
        figure('WindowState','maximized');
    end
    
    % Randomly initialize agents
    for agent = 1: No_of_evaders + No_of_pursuers
        % 
        % if agent == 3 
        %     agents(agent).x = .07;
        %     agents(agent).y = .07;
        % end
        % if agent == 2 
        %     agents(agent).x = .25;
        %     agents(agent).y = .92;
        % end
        % if agent == 1 
        %     agents(agent).x = .24;
        %     agents(agent).y = .94;
        % end
%         
        agents(agent).x = initial_matrix(1, agent);
        agents(agent).y = initial_matrix(2, agent);
        % while true
        %     counter = 0;
        %     agents(agent).x = rand(1);
        %     agents(agent).y = rand(1);
        % 
        %     for k1 = -2:2
        %         for k2 = -2:2
        %             if M(xy_index_check(ceil(agents(agent).x * N_matrix) + k1, N_matrix), xy_index_check(ceil(agents(agent).y * N_matrix) + k2, N_matrix)) == 0
        %                 counter = counter + 1;
        %             end
        %         end
        %     end
        % 
        %     if counter == 25
        %         break
        %     end
        % end

        % %prevents wall spawns
        % wall_found = false;
        % while ~wall_found
        %     wall_found = false;
        %     agents(agent).x = rand(1);
        %     agents(agent).y = rand(1);
        %     for k1 = -1:
        %         for k2 = -1:1
        %             if M_new(xy_index_check(ceil(agents(agent).x * N_matrix) + k1, N_matrix), xy_index_check(ceil(agents(agent).y * N_matrix) + k2, N_matrix)) == 1
        %                 wall_found = true;
        %                 break
        %             end
        %         end
        %         if wall_found
        %             break
        %         end
        %     end
        % end

        agents(agent).vx = 0; % Placeholder value for linear velocity
        agents(agent).vy = 0; % Placeholder value for linear velocity in y-direction
        agents(agent).ax = 0; % Placeholder value for linear acceleration
        agents(agent).ay = 0; % Placeholder value for linear acceleration in y-direction
        % agents(agent).theta = rand(1)*2*pi;

        agents(agent).theta = initial_matrix(3, agent);
        agents(agent).omega = 0; % Placeholder value for angular velocity
        agents(agent).alpha = 0; % Placeholder value for angular acceleration

        if agent <= No_of_evaders
            agents(agent).color = 2;
        else
            agents(agent).color = 1;
        end
        
        agents(agent).exploration_map = fieldmap;
    end

    for agent = 1:length(agents)
        disp("agent no. " + num2str(agent))
        disp(agents(agent).x)
        disp(agents(agent).y)
        disp(agents(agent).theta)
    end

    
    initial_matrix = [agents(1).x, agents(2).x, agents(3).x;
        agents(1).y, agents(2).y, agents(3).y;
        agents(1).theta, agents(2).theta, agents(3).theta];
    initial_M = M;

    
    % Simulation parameters
    h = 0.005; % Step length
    n = 500; % Number of iterations
    %n = 1000000;
    kill_radius = .0075; % spikad-ish

    % variables
    v_max = 3;    % spikad, vi ändrar inte 3.5 för basmodellen
    max_acceleration = 3*v_max;
    max_deacceleration = -14.8*v_max;
    factor_reach_max_angular_velocity = 3;
    max_radie = 25;

    % iter_ahead_own = floor(abs(v_max / (max_deacceleration*h)));

    % return values
    iterations_until_completion = 0;
    reason_for_exit = "NULL";
    pursuer_wall_deaths = 0;
    evader_wall_deaths = 0;
    % Game loop
    for K = 1:n
        % if K > 100
        %     pause(.5)
        % end
        % disp("------------------------")
        iterations_until_completion = iterations_until_completion + 1;
        % if K > 1
        %     break
        % end
        % pause(.25)

%         disp("init")
%         for agent_index = 1:length(agents)
%             disp("agent number:" + num2str(agent_index))
%             disp(agents(agent_index).x)
%             disp(agents(agent_index).y)
%         end
    
        agents = update_bresenham(agents, fieldmap, length(agents), N_matrix);
    
        % iterate agents again to choose pursuers tx, ty and evaders px, py
        for agent = 1:length(agents)
            if agents(agent).color == 1     % pursuer
                if isempty(agents(agent).px_list)   % out of sight
                    continue
                else
                    distance_vector = zeros(1,length(agents(agent).px_list));
                    for pursuer = 1:length(agents(agent).px_list)
                        distance_vector(pursuer) = norm([agents(agent).px_list(pursuer); agents(agent).py_list(pursuer)] ...
                        - [agents(agent).x; agents(agent).y]);
                    end
                    [~, index] = min(distance_vector);

                    % for evader sighting
                    agents(agent).px = agents(agent).px_list(index);
                    agents(agent).py = agents(agent).py_list(index);
                    for potential_evader = 1:length(agents)
                        if agents(potential_evader).color == 2
                            if norm([agents(potential_evader).x - agents(agent).px, agents(potential_evader).y - agents(agent).py]) < kill_radius*.1
                                    agents(agent).evader_index = potential_evader;
                                    break
                            end
                        end
                    end
                end
            end
        end 
    
        % calculate movement
        for agent = 1:length(agents)
    
            max_angular_velocity = 7*pi;
            % crash prevention
            x_index = xy_index_check(ceil(agents(agent).x * N_matrix), N_matrix);
            y_index = xy_index_check(ceil(agents(agent).y * N_matrix), N_matrix);

%             disp("agent number: " + num2str(agent) + "x_index1: " + num2str(x_index))
%             disp("agent number: " + num2str(agent) + "y_index1: " + num2str(y_index))
    
            %if ~agents(agent).bresenham || agents(agent).color == 2 % ~bresenham or evader
            % Update Matrix index
            min_pot = inf;
            if ~increment_boolean
                increment = 0;
            end
            max_pot = 1;   
            
% disp(agents);
% disp('----------------');
% disp(max_radie);
% disp('----------------');
% disp(x_index);
% disp('----------------');
% disp(y_index);
% disp('----------------');
% disp(M);
% disp('----------------');
% disp(agent);
% disp('----------------');
% disp(increment);
% disp('----------------');
% disp(max_pot);
% disp('----------------');
% disp(N_matrix);
% disp('----------------');
% disp(fieldmap);
% disp('----------------');
% disp(kill_radius);
% disp('----------------');
% disp(sighting_target_boolean);
% disp('----------------');
% disp(jesus_boolean);
% disp('----------------');
% disp(smart_target_boolean);
% disp('----------------');
% disp(h');
% disp('----------------');
% disp(iterations_ahead);
% disp('----------------');
% disp(increment_boolean);
% disp('----------------');
% disp(factor_reach_max_angular_velocity);
% disp('----------------');
% disp(max_angular_velocity);
% disp('----------------');
% disp(max_iter);
% disp('----------------');
% disp(crystal_tp_ahead);

            [agents, M] = set_target_points(agents, max_radie, x_index, y_index, M, agent, increment, max_pot, N_matrix, fieldmap, kill_radius, sighting_target_boolean, jesus_boolean, smart_target_boolean, h, iterations_ahead, increment_boolean, factor_reach_max_angular_velocity, max_angular_velocity, max_iter, crystal_tp_ahead, pursuit_crystal, scale_factor, decay_rate);
    
            % Calculate the angle between the agent's current heading and the direction to the target
            targetAngle = atan2(agents(agent).ty - agents(agent).y, agents(agent).tx - agents(agent).x);
            error = targetAngle - agents(agent).theta; %skillnaden i vinkel
            error = atan2(sin(error), cos(error)); %se till att det skall vara mellan -pi och pi
            
            if any(error>0.1*pi) || any(error<-0.1*pi)
                agents(agent).alpha = sign(error)*max_angular_velocity/(factor_reach_max_angular_velocity*h);
            else
                agents(agent).alpha = 10*error*max_angular_velocity/(factor_reach_max_angular_velocity*h);
                if any(abs((-sign(error)*max_angular_velocity/factor_reach_max_angular_velocity + agents(agent).omega))*h > error)
    
                    agents(agent).omega = sign(error)*abs(agents(agent).theta - targetAngle)/h;
                    % wanted_omega*h + sign(error)*abs(theta_0 -targetAngle) = 0 
                    % theta_0 = targetAngle -wanted_omega*h;
                    agents(agent).alpha = 0;
    
    
                    % if any(max_angular_velocity/factor_reach_max_angular_velocity > omega) && any(abs(error)<.1)
                    %     omega = 0;
                    %     alpha = 0;
                    % end
                end
            end
    
            % Update angular velocity
            agents(agent).omega = agents(agent).omega + agents(agent).alpha * h;

%             disp("agent number: " + num2str(agent) + " alpha: " + num2str(agents(agent).alpha))
            if agents(agent).omega > max_angular_velocity
                agents(agent).omega = max_angular_velocity;
            end
            if agents(agent).omega < -max_angular_velocity
                agents(agent).omega = -max_angular_velocity;
            end
             
            % Update angular position
            agents(agent).theta = agents(agent).theta + agents(agent).omega * h;
            agents(agent).theta = agents(agent).theta + pi;
            agents(agent).theta = mod(agents(agent).theta, 2*pi);
            agents(agent).theta = agents(agent).theta - pi;

%             disp("agent number: " + num2str(agent) + " omega: " + num2str(agents(agent).omega))
        
            % Update linear velocity components using the new acceleration
            v = sqrt(agents(agent).vx^2 + agents(agent).vy^2);
    
            % if gradient_linear_acceleration
            % agents = gradient_linear_acceleration_func(x_index, y_index, fieldmap, v_max, N_matrix, agents, agent);
            % target_speed = agents(agent).target_speed;
            
            % if gradient_linear_acceleration
            open_area = true;
            for_loop_flag = true;
            for k1 = [-3, 0, 3]
                for k2 = [-3, 0, 3]
                    if M_new(xy_index_check(ceil(agents(agent).x * N_matrix) + k1, N_matrix), xy_index_check(ceil(agents(agent).y * N_matrix) + k2, N_matrix)) == 1
                        open_area = false;
                        for_loop_flag = false;
                        break
                    end
                end 
                if ~for_loop_flag
                    break
                end
            end

            % disp(agents(agent).px)
            if open_area
                [crash_flag, agents] = iter_ahead_function(3, v + h*max_acceleration, agents(agent).x + h*v*cos(agents(agent).theta), agents(agent).y + h*v*sin(agents(agent).theta), max_iter, crystal_tp_ahead, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap);
                if sighting_target_boolean && ~isempty(agents(agent).px) && norm([agents(agent).tx - agents(agent).x, agents(agent).ty - agents(agent).y]) < norm([agents(agent).x_list(end-crystal_tp_ahead) - agents(agent).x, agents(agent).y_list(end - crystal_tp_ahead) - agents(agent).y]) && agents(agent).color == 1
                    v = v + h*max_deacceleration;
                else 
                    if ~crash_flag
                        v = v + h*max_acceleration;
                    else
                        [crash_flag, ~] = iter_ahead_function(1, v, agents(agent).x + h*v*cos(agents(agent).theta), agents(agent).y + h*v*sin(agents(agent).theta), max_iter, crystal_tp_ahead, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap);
                        if ~crash_flag
                            % v = v;
                        else
                            % disp(rand(1))
                            % disp(1)
                            v = v + h*max_deacceleration;
                        end
                    end
                end
            else
                if ~iter_ahead_function(1, v + h*max_acceleration, agents(agent).x + h*v*cos(agents(agent).theta), agents(agent).y + h*v*sin(agents(agent).theta), max_iter, crystal_tp_ahead, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap)
                    v = v + h*max_acceleration;
                else
                    if ~iter_ahead_function(0, v, agents(agent).x + h*v*cos(agents(agent).theta), agents(agent).y + h*v*sin(agents(agent).theta), max_iter, crystal_tp_ahead, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap)
                        % v = v;
                    else
                        v = v + h*max_deacceleration;
                    end
                end
            end

            % if v < target_speed
            % 
            %     crash_flag = true;
            %     for_loop_flag = true;
            %     for k1 = -3:3
            %         for k2 = -3:3
            %             if M_new(xy_index_check(ceil(agents(agent).x * N_matrix) + k1, N_matrix), xy_index_check(ceil(agents(agent).y * N_matrix) + k2, N_matrix)) == 1
            %                 crash_flag = iter_ahead_function(0, v + h*max_acceleration, agents(agent).x + h*v*cos(agents(agent).theta), agents(agent).y + h*v*sin(agents(agent).theta), max_iter, crystal_tp_ahead, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap);
            % 
            %                 if ~crash_flag
            %                     v = v + h*max_acceleration;
            %                 end
            % 
            %                 for_loop_flag = false;
            %                 break
            %             end
            %         end
            %         if ~for_loop_flag
            %             break
            %         end
            %     end
            % 
            %     if for_loop_flag
            %         [crash_flag, agents] = iter_ahead_function(2, v + h*max_acceleration, agents(agent).x + h*v*cos(agents(agent).theta), agents(agent).y + h*v*sin(agents(agent).theta), max_iter, crystal_tp_ahead, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap);
            %         if ~crash_flag
            %             v = v + h*max_acceleration;
            %         else
            %             if iter_ahead_function(1, v, agents(agent).x + h*v*cos(agents(agent).theta), agents(agent).y + h*v*sin(agents(agent).theta), max_iter, crystal_tp_ahead, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap)
            %                 v = v + h*max_deacceleration;
            %             else
            %                 % v = v
            %             end
            %         end
            %     end
            % 
            %     % ~krashar med accelererade radie 1
            %     %   accelerera
            %     % krashar med accelererade radie 1
            %     %   kolla krash utan acceleration radie 0
            %     %       krash utan acceleration radie 0
            %     %           deaccacellerera
            %     %       ~krash utan acceleration radie 0
            %     %           accelerera inte
            %     % 
            %     % 
            %     % [extrapolated_wall_crash_acceleration, agents, crash_flag] = own_crystal(agents, agent, fieldmap, v, h, potential_acceleration, iter_ahead_own, max_angular_velocity, factor_reach_max_angular_velocity, v_max, kill_radius);
            %     % 
            %     % if crash_flag
            %     %     v = v + h*max_deacceleration;
            %     % else
            %     %     v = v + h*extrapolated_wall_crash_acceleration;
            %     % end
            %     % % if v > target_speed
            %     % %     v = target_speed;
            %     % % end
            % 
            % else
            % 
            %     crash_flag = true;
            %     for_loop_flag = true;
            %     for k1 = -3:3
            %         for k2 = -3:3
            %             if M_new(xy_index_check(ceil(agents(agent).x * N_matrix) + k1, N_matrix), xy_index_check(ceil(agents(agent).y * N_matrix) + k2, N_matrix)) == 1
            %                 crash_flag = iter_ahead_function(0, v, agents(agent).x + h*v*cos(agents(agent).theta), agents(agent).y + h*v*sin(agents(agent).theta), max_iter, crystal_tp_ahead, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap);
            %                 for_loop_flag = false;
            %                 v = v + h*max_acceleration;
            %                 break
            %             end
            %         end
            %         if ~for_loop_flag
            %             break
            %         end
            %     end
            % 
            %     if crash_flag
            %         if iter_ahead_function(1, v, agents(agent).x + h*v*cos(agents(agent).theta), agents(agent).y + h*v*sin(agents(agent).theta), 20, 6, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap)
            %             v = v + h*max_deacceleration;
            %         end
            %     end
            % 
            % end
    
            if v > v_max
                v = v_max;
            end

            % <=
            if v <= 0
                v = v_max * 0.005;
            end

%             disp("agent number: " + num2str(agent) + " theta: " + num2str(agents(agent).theta))
            
    
            agents(agent).vx = v * cos(agents(agent).theta);
            agents(agent).vy = v * sin(agents(agent).theta);
        
            % Update linear position

            agents(agent).previous_x = agents(agent).x;
            agents(agent).previous_y = agents(agent).y;





% 
%             disp("agent number: " + num2str(agent) + " x: " + num2str(agents(agent).x))
%             disp("agent number: " + num2str(agent) + " y: " + num2str(agents(agent).y))
%             disp("agent number: " + num2str(agent) + " vx: " + num2str(agents(agent).vx))
%             disp("agent number: " + num2str(agent) + " vy: " + num2str(agents(agent).vy)) 

            agents(agent).x = agents(agent).x + agents(agent).vx * h;
            agents(agent).y = agents(agent).y + agents(agent).vy * h;

%             disp("agent number: " + num2str(agent) + " x: " + num2str(agents(agent).x))
%             disp("agent number: " + num2str(agent) + " y: " + num2str(agents(agent).y))



%             disp("-----")
%             disp("agent number: " + num2str(agent) + " x: " + num2str(agents(agent).x))
%             disp("agent number: " + num2str(agent) + " y: " + num2str(agents(agent).y))
%             disp("")
%             disp(agents(agent).x * N_matrix)
%             disp(agents(agent).y * N_matrix)
%             disp("")
%             disp(ceil(agents(agent).x * N_matrix))
%             disp(ceil(agents(agent).y * N_matrix))
%             disp("")
%             disp(xy_index_check(ceil(agents(agent).x * N_matrix), N_matrix))
%             disp(xy_index_check(ceil(agents(agent).y * N_matrix), N_matrix))
%             disp("-----")
        
            % Update Matrix index
            x_index = xy_index_check(ceil(agents(agent).x * N_matrix), N_matrix);
            y_index = xy_index_check(ceil(agents(agent).y * N_matrix), N_matrix);
% 
%             disp("agent number: " + num2str(agent) + "x_index2: " + num2str(x_index))
%             disp("agent number: " + num2str(agent) + "y_index2: " + num2str(y_index))
        
            color_M(x_index, y_index) = agents(agent).color; % Use color to represent agents

            % disp(agents(agent).x)
            % disp(agents(agent).x_index_memo)
            % disp(agents(agent).y)
            % disp(agents(agent).y_index_memo)
            % memory
            if ~isempty(agents(agent).x_index_memo)
                agents(agent).x_index_memo(end+1) = x_index;
                agents(agent).y_index_memo(end+1) = y_index;
            else
                agents(agent).x_index_memo = x_index;
                agents(agent).y_index_memo = y_index;
            end

        
            if length(agents(agent).x_index_memo)>50
                color_M(agents(agent).x_index_memo(1),agents(agent).y_index_memo(1)) = 0;
                agents(agent).x_index_memo(1) = [];
                agents(agent).y_index_memo(1) = [];
            end
    
        end
        
        [agents, wall_pursuer_counter, wall_evader_counter] = wall_killing(agents, N_matrix, base_M);
        pursuer_wall_deaths = pursuer_wall_deaths + wall_pursuer_counter;
        evader_wall_deaths = evader_wall_deaths + wall_evader_counter;
        agents = pursuer_killing(agents, kill_radius);
        
        if visualization_boolean
            % place target points
            for agent = 1:length(agents)
                if (agents(agent).bresenham && agents(agent).color == 1) && ~smart_target_boolean  % pursuer chasing evader
                    % place predicted
                    % M(agents(agent).x0_pred, agents(agent).y0_pred) = 3;

                    % % place predicted list
                    % for cell = 1:length(agents(agent).x_list)
                    %      M(xy_index_check(ceil(agents(agent).x_list(cell)*N_matrix), N_matrix), xy_index_check(ceil(agents(agent).y_list(cell)*N_matrix), N_matrix)) = 3;
                    % end
                    continue
                end

                % Update Matrix index
                x_index = ceil(agents(agent).tx * N_matrix);
                y_index = ceil(agents(agent).ty * N_matrix);
                % 
                % % Check if indices are within bounds
                % if all(x_index-1 > 0) && all(x_index+1 <= N_matrix) && all(y_index-1 > 0) && all(y_index+1 <= N_matrix)
                %     color_M(x_index, y_index) = 4; % Use color to represent agents
                % end
                % 
                % if agents(agent).color == 1 && ~isempty(agents(agent).px)
                %     % disp([ceil(N_matrix*agents(agent).px)', ceil(N_matrix*agents(agent).py)'])
                %     for p = 1:length(agents(agent).px)
                %         color_M(ceil(agents(agent).px(p) * N_matrix), ceil(agents(agent).py(p) * N_matrix)) = 1;
                %     end
                % end
                % place predicted
                % M(agents(agent).x0_pred, agents(agent).y0_pred) = 3;
                % % place predicted list
                % for cell = 1:length(agents(agent).x_list)
                %      M(xy_index_check(ceil(agents(agent).x_list(cell)*N_matrix), N_matrix), xy_index_check(ceil(agents(agent).y_list(cell)*N_matrix), N_matrix)) = 3;
                % end

            end
            % 
            % for agent = 1:length(agents)own
            %     if ~isempty(agents(agent).temp)
            %         for cell = 1:length(agents(agent).temp(:,1))
            %             M(agents(agent).temp(cell, 1), agents(agent).temp(cell, 2)) = 5;
            %         end
            %     end
            % end


            % Display the playing field matrix
            exp_walls = agents(end).exploration_map >= 1;
            exp_walls = 3*exp_walls;
            imagesc(M + color_M);
            % M = base_M;
            colormap([1 1 1; 1 0 0; 0 0 1; 0 0 0]); % White for empty, Red for red, Blue for blue, Black for wall, pink for target point

            % Counter for remaining blue agents
            remaining_blue_agents = sum([agents.color] == 2);
            % Construct the title string
            title_str = [""];
            % title_str = ["Current crystal tp iter ahead is " + num2str(crystal_tp_ahead) + ", Current density is " + num2str(obstacle_density) + "." + " Iteration no. " + num2str(K)];
            
            % Set the figure title with the constructed string
            title(title_str, 'HorizontalAlignment', 'center');
            axis square;
            pause(0.0001); % Pause for visualization

            % remove target points and predicted
            for agent = 1:length(agents)
                % 
                % if agents(agent).color == 1 && ~isempty(agents(agent).px)
                %     for p = 1:length(agents(agent).px)
                %         color_M(ceil(agents(agent).px(p)*N_matrix), ceil(agents(agent).py(p)*N_matrix)) = 0;
                %     end
                % end
                % M(agents(agent).x0_pred, agents(agent).y0_pred) = 0;

                % remove predicted list
                % for cell = 1:length(agents(agent).x_list)
                %      M(xy_index_check(ceil(agents(agent).x_list(cell)*N_matrix), N_matrix), xy_index_check(ceil(agents(agent).y_list(cell)*N_matrix), N_matrix)) = 0;
                % end


                if agents(agent).bresenham && agents(agent).color == 1 && ~smart_target_boolean  % pursuer chasing evader
                    continue
                end
                % Update Matrix index
                x_index = ceil(agents(agent).tx * N_matrix);
                y_index = ceil(agents(agent).ty * N_matrix);
    
                % % Check if indices are within bounds
                % if any(x_index > 0) && any(x_index <= N_matrix) && any(y_index > 0) && any(y_index <= N_matrix)
                %     color_M(x_index, y_index) = 0; % Use color to represent agents
                % end
            end
        end
        
        if ~visualization_boolean
            % Counter for remaining blue agents
            remaining_blue_agents = sum([agents.color] == 2);
        end

%         End the game if there are no remaining blue agents
        if remaining_blue_agents == 0
            reason_for_exit = "All evaders are dead";
            writematrix(initial_matrix, "initial_matrix.csv")
            writematrix(initial_M, "initial_m.csv")
            break;
        end
        % 
        % % break if all pursuers are dead
        if sum([agents.color] == 1) == 0
            reason_for_exit = "All pursuers are dead";
            writematrix(initial_matrix, "initial_matrix.csv")
            writematrix(initial_M, "initial_m.csv")
            break;
        end

        % if agents(agent).color == 1
        %     print_info(agents, agent, fieldmap, N_matrix);
        % end
        
        if K == n
            %writematrix(initial_matrix, "initial_matrix.csv")
            %writematrix(initial_M, "initial_m.csv")
            reason_for_exit = "Maximum iterations reached";
        end

        for agent = 1:length(agents)
            agents(agent).previous_x = agents(agent).x;
            agents(agent).previous_y = agents(agent).y;
        
            % if agents(agent).color == 1
            %     disp(K)
            %     disp("px_list")
            %     disp(agents(agent).px_list)
            % end
        end

    	frame = getframe(gcf);
        frames{K} = frame;

    end

    % writematrix(initial_matrix, "initial_matrix.csv")
    % writematrix(initial_M, "initial_m.csv")

    % Create a VideoWriter object
    writerObj = VideoWriter('testing_increment_true.avi');
    writerObj.FrameRate = 30; % Adjust frame rate as needed
    open(writerObj);

    % % Write frames to the video file
    for i = 1:length(frames)
        % Convert frame data to uint8
        frame_data = im2uint8(frames{i}.cdata);
        writeVideo(writerObj, frame_data);
    end

    % % Close the VideoWriter object
    close(writerObj);

    writematrix(initial_matrix, "initial_matrix.csv")
    writematrix(initial_M, "initial_m.csv")

end
