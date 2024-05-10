function [agents] = crystal_tp(agents, agent, h, fieldmap, N_matrix, tp_list, factor_reach_max_angular_velocity, max_angular_velocity, kill_radius, max_iter, crystal_tp_ahead)
agents(agent).x_list = [];
agents(agent).y_list = [];

for k = size(tp_list, 1):-1:1

        % if agent == 1
    %     disp(k)
    % end
    % 
    v0 = sqrt(agents(agent).vx^2 + agents(agent).vy^2);
    x0 = agents(agent).x;
    y0 = agents(agent).y;
    theta_0 = agents(agent).theta;
    omega = agents(agent).omega;
    tx = tp_list(k, 1)/N_matrix;
    ty = tp_list(k, 2)/N_matrix;
%      if agents(agent).color == 1
%         disp("tx som testas: "+num2str(tx*N_matrix))
%         disp("ty som testas: "+num2str(ty*N_matrix))
%      end
    if ~bresenham(xy_index_check(ceil(agents(agent).x*N_matrix), N_matrix), xy_index_check(ceil(agents(agent).y*N_matrix), N_matrix), xy_index_check(ceil(tx*N_matrix), N_matrix), xy_index_check(ceil(ty*N_matrix), N_matrix), fieldmap, 0)
        continue
    end
    
    if agents(agent).exploration_map(xy_index_check(ceil(agents(agent).x*N_matrix), N_matrix), xy_index_check(ceil(agents(agent).y*N_matrix), N_matrix)) == 1
        [crash_flag, agents] = iter_ahead_function(0, v0, x0, y0, max_iter, crystal_tp_ahead, tx, ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap);
    else
        [crash_flag, agents] = iter_ahead_function(0, v0, x0, y0, max_iter, crystal_tp_ahead, tx, ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, agents(agent).exploration_map);
    end
    if ~crash_flag
        if agents(agent).color == 1
%             disp(crash_flag)
%             disp("tx som lyckades: "+num2str(tx*N_matrix))
%             disp("ty som lyckades: "+num2str(ty*N_matrix))
        end
    
        agents(agent).tx = tx;
        agents(agent).ty = ty;
        return

    end

    % % direction of angular acceleration
    % % Calculate the angle between the agent's current heading and the direction to the target
    % targetAngle = atan2(ty - y0, tx- x0);
    % error = targetAngle - theta_0; %skillnaden i vinkel
    % error = atan2(sin(error), cos(error)); %se till att det skall vara mellan -pi och pi
    % 
    % if any(error>0.1*pi) || any(error<-0.1*pi)
    %     alpha = sign(error)*max_angular_velocity/(factor_reach_max_angular_velocity*h);
    % else
    %     alpha = 10*error*max_angular_velocity/(factor_reach_max_angular_velocity*h);
    %     if any(max_angular_velocity/factor_reach_max_angular_velocity > omega) && any(abs(error)<.03)
    %         omega = 0;
    %         alpha = 0;
    %     end
    % 
    % end
    % 
    % % Update angular velocity
    % omega = omega + alpha * h;
    % target_point_found = false;
    % for j = 1:tp_crystal_iter
    %     crash_flag = false;
    % 
    %     % direction of angular acceleration
    %     % Calculate the angle between the agent's current heading and the direction to the target
    %     targetAngle = atan2(ty - y0, tx - x0);
    %     error = targetAngle - theta_0; %skillnaden i vinkel
    %     error = atan2(sin(error), cos(error)); %se till att det skall vara mellan -pi och pi
    % 
    %     if norm([ty- y0, tx - x0]) < 2*.0075
    %         target_point_found = true;
    %     end
    % 
    %     if any(error>0.1*pi) || any(error<-0.1*pi)
    %         alpha = sign(error)*max_angular_velocity/(factor_reach_max_angular_velocity*h);
    %     else
    %         alpha = 10*error*max_angular_velocity/(factor_reach_max_angular_velocity*h);
    %         if any(max_angular_velocity/factor_reach_max_angular_velocity > omega) && any(abs(error)<.1)
    %             omega = 0;
    %             alpha = 0;
    %         end
    %     end
    % 
    %     if target_point_found
    %         alpha = 0;
    %     end
    % 
    %     % Update angular velocity
    %     omega = omega + alpha * h;
    % 
    %     if omega > max_angular_velocity
    %         omega = max_angular_velocity;
    %     end
    %     if omega < -max_angular_velocity
    %         omega = -max_angular_velocity;
    %     end
    % 
    %     % omega = omega + alpha * h;
    %     theta_0 = theta_0 + h * omega;
    % 
    %     theta_0 = theta_0 + pi;
    %     theta_0 = mod(theta_0, 2*pi);
    %     theta_0 = theta_0 - pi;
    % 
    %     % komposantuppdelning
    %     vx = v0*cos(theta_0);
    %     vy = v0*sin(theta_0);
    %     x0 = x0 + vx * h;
    % 
    %     agents(agent).x_list(end+1) = x0;
    %     y0 = y0 + vy * h;
    %     agents(agent).y_list(end+1) = y0;
    % 
    %     if fieldmap(xy_index_check(ceil(x0*N_matrix), N_matrix), xy_index_check(ceil(y0*N_matrix), N_matrix)) == 1
    %         crash_flag = true;
    %         agents(agent).x_list = [];
    %         agents(agent).y_list = [];
    % 
    %     end
    % 
    % 
    %     if crash_flag
    %         break
    %     end
    % end
    % 
    % if crash_flag
    %     continue
    % end
    % 
    % agents(agent).tx = tp_list(k, 1)/N_matrix;
    % agents(agent).ty = tp_list(k, 2)/N_matrix;
    % % if agents(agent).color == 1
    % %     disp("LYCKAS")
    % % end
    % return
end

% disp("Out of points in tp_list")

end

