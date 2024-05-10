function [crash_flag, agents] = iter_ahead_function(radius, v, x0, y0, max_iter, tp_crystal_iter_ahead, tx, ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap)
    theta_0 = agents(agent).theta;
    omega = agents(agent).omega;
    
    % agents(agent).x_list = [];
    % agents(agent).y_list = [];

    crash_flag = false;
    for iter = 1:max_iter
        % direction of angular acceleration
        % Calculate the angle between the agent's current heading and the direction to the target
        targetAngle = atan2(ty - y0, tx - x0);
        error = targetAngle - theta_0; %skillnaden i vinkel
        error = atan2(sin(error), cos(error)); %se till att det skall vara mellan -pi och pi
        
        if any(error>0.1*pi) || any(error<-0.1*pi)
            alpha = sign(error)*max_angular_velocity/(factor_reach_max_angular_velocity*h);
        else
            alpha = 10*error*max_angular_velocity/(factor_reach_max_angular_velocity*h);
            if any(abs((-sign(error)*max_angular_velocity/factor_reach_max_angular_velocity + omega))*h > error)

                omega = sign(error)*abs(theta_0 - targetAngle)/h;
                % wanted_omega*h + sign(error)*abs(theta_0 -targetAngle) = 0 
                % theta_0 = targetAngle -wanted_omega*h;
                alpha = 0;


                % if any(max_angular_velocity/factor_reach_max_angular_velocity > omega) && any(abs(error)<.1)
                %     omega = 0;
                %     alpha = 0;
                % end
            end
        end
    
        % Update angular velocity
        omega = omega + alpha * h;

        if omega > max_angular_velocity
            omega = max_angular_velocity;
        end
        if omega < -max_angular_velocity
            omega = -max_angular_velocity;
        end


        % omega = omega + alpha * h;
        theta_0 = theta_0 + h * omega;

        theta_0 = theta_0 + pi;
        theta_0 = mod(theta_0, 2*pi);
        theta_0 = theta_0 - pi;
    
        % komposantuppdelning
        vx = v*cos(theta_0);
        vy = v*sin(theta_0);
        x0 = x0 + vx * h;
        agents(agent).x_list(end+1) = x0;
        y0 = y0 + vy * h;
        agents(agent).y_list(end+1) = y0;
    
        for k1 = -radius:radius
            for k2 = -radius:radius
                if fieldmap(xy_index_check(ceil(x0*length(fieldmap)) + k1, length(fieldmap)), xy_index_check(ceil(y0*length(fieldmap)) + k2, length(fieldmap))) == 1
                    crash_flag = true;
                    % agents(agent).x_list = [];
                    % agents(agent).y_list = [];
                    return
                end
            end
        end
    
        if norm([x0 - tx, y0 - ty]) < 2*kill_radius
            for ii = 1:tp_crystal_iter_ahead
                % omega = omega + alpha * h;
                theta_0 = theta_0 + h * omega;
        
                theta_0 = theta_0 + pi;
                theta_0 = mod(theta_0, 2*pi);
                theta_0 = theta_0 - pi;
        
                % komposantuppdelning
                vx = v*cos(theta_0);
                vy = v*sin(theta_0);
                x0 = x0 + vx * h;
                agents(agent).x_list(end+1) = x0;
                y0 = y0 + vy * h;
                agents(agent).y_list(end+1) = y0;
    
                for k1 = -radius:radius
                    for k2 = -radius:radius
                        if fieldmap(xy_index_check(ceil(x0*length(fieldmap)) + k1, length(fieldmap)), xy_index_check(ceil(y0*length(fieldmap)) + k2, length(fieldmap))) == 1
                            crash_flag = true;
                            % agents(agent).x_list = [];
                            % agents(agent).y_list = [];
                            return
                        end
                    end
                end
            end
            return
        end
    end
end
