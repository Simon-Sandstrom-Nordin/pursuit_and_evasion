function [agents, M] = set_target_points(agents, max_radie, x_index, y_index, M, agent, increment, max_pot, N_matrix, fieldmap, kill_radius, sighting_target_boolean, jesus_boolean, smart_target_boolean, h, iterations_ahead, increment_boolean, factor_reach_max_angular_velocity, max_angular_velocity, max_iter, crystal_tp_ahead, pursuit_crystal, increment_scale_factor, increment_decay_rate, jesus_scale_factor)  
     
    tp_list = [];

    if agents(agent).bresenham  % pursuer that sights an evader
        if agents(agent).color == 1

            % choice of algorithm
            if smart_target_boolean
                agents = smart_target_function(agents, agent, h, fieldmap, N_matrix, iterations_ahead);
                
            else

                agents(agent).tx = agents(agent).px;
                agents(agent).ty = agents(agent).py;
               
            end

            tx = agents(agent).tx;
            ty = agents(agent).ty;
            v = sqrt(agents(agent).vx^2 + agents(agent).vy^2);
            x0 = agents(agent).x;
            y0 = agents(agent).y;

            [crash_flag, agents] = iter_ahead_function(0, v, x0, y0, max_iter, pursuit_crystal, tx, ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap);
            if ~crash_flag
                return
            end

        end

    end

    if ~agents(agent).bresenham && agents(agent).color == 1 && ~isempty(agents(agent).px) % sightless pursuer
        if sighting_target_boolean
            agents = sighting_target_function(agents, agent, kill_radius, fieldmap, max_iter, crystal_tp_ahead, max_angular_velocity, factor_reach_max_angular_velocity, h);
 
            if ~isempty(agents(agent).tx)
                return
            end
        end
    end

    % only if not returned at this point / evaders only use this part
    min_pot = inf;
    radie = 2;

    % 1 - 0 från 0 till max_radie
    max_distance = sqrt(2*(max_radie/N_matrix)^2);

    while radie < max_radie % is_searching && 
        
        k = -radie;
        for j = -radie: radie
                            
            %undviker krash
            if (x_index+k <= 1) || (x_index+k >= length(M)) || (y_index+j <= 1) || (y_index+j >= length(M))
                continue
            end
            
            distance = 0;

            if agents(agent).bresenham
                for pursuer = 1:length(agents(agent).px_list)
                    distance = distance + norm([agents(agent).px_list(pursuer), agents(agent).py_list(pursuer)] - [(x_index+k)/N_matrix, (y_index+j)/N_matrix])/sqrt(2);
                end
                distance = distance/length(agents(agent).px_list); % / (2*);
            end
            
            if agents(agent).exploration_map(x_index+k, y_index+j) + 1 - 2*distance < min_pot
                min_pot = agents(agent).exploration_map(x_index+k, y_index+j) + 1 - 2*distance;
                min_pot_xi = x_index+k;
                min_pot_yi = y_index+j;
                tp_list(end+1, :) = [min_pot_xi, min_pot_yi];
            end

            if agents(agent).exploration_map(x_index+k, y_index+j) <= max_pot
                 if agents(agent).color == 1
                     agents(agent).exploration_map(x_index+k, y_index+j) = agents(agent).exploration_map(x_index+k, y_index+j) + increment;
                      if agents(agent).exploration_map(x_index+k, y_index+j) > 1
                          agents(agent).exploration_map(x_index+k, y_index+j) = 1;
                      end
                 end
            end
           
            if size(tp_list, 1) < 1
                tp_list(end+1, :) = [x_index+k, y_index+j];
            end
        end
      
        j = radie;
        for k = -radie + 1: radie
                                 
            %undviker krash
            if (x_index+k <= 1) || (x_index+k >= length(M)) || (y_index+j <= 1) || (y_index+j >= length(M))
                 continue
            end

            distance = 0;

            if agents(agent).bresenham
                for pursuer = 1:length(agents(agent).px_list)
                    distance = distance + norm([agents(agent).px_list(pursuer), agents(agent).py_list(pursuer)] - [(x_index+k)/N_matrix, (y_index+j)/N_matrix])/sqrt(2);
                end
                distance = distance/length(agents(agent).px_list); % / (2*);
            end
               
            if agents(agent).exploration_map(x_index+k, y_index+j) + 1 - 2*distance < min_pot
                min_pot = agents(agent).exploration_map(x_index+k, y_index+j) + 1 - 2*distance;
                min_pot_xi = x_index+k;
                min_pot_yi = y_index+j;
                tp_list(end+1, :) = [min_pot_xi, min_pot_yi];
            end

            if agents(agent).exploration_map(x_index+k, y_index+j) <= max_pot
                 if agents(agent).color == 1
                     agents(agent).exploration_map(x_index+k, y_index+j) = agents(agent).exploration_map(x_index+k, y_index+j) + increment;
                      if agents(agent).exploration_map(x_index+k, y_index+j) > 1
                          agents(agent).exploration_map(x_index+k, y_index+j) = 1;
                      end
                 end
            end
        
            if size(tp_list, 1) < 2
                tp_list(end+1, :) = [x_index+k, y_index+j];
            end

        end
        k = radie;
        for j = radie - 1: -1 :-radie
        
            %undviker krash
            if (x_index+k <= 1) || (x_index+k >= length(M)) || (y_index+j <= 1) || (y_index+j >= length(M))
                continue
            end
            
            distance = 0;

            if agents(agent).bresenham
                for pursuer = 1:length(agents(agent).px_list)
                    distance = distance + norm([agents(agent).px_list(pursuer), agents(agent).py_list(pursuer)] - [(x_index+k)/N_matrix, (y_index+j)/N_matrix])/sqrt(2);
                end
                distance = distance/length(agents(agent).px_list); % / (2*);
            end
            
            if agents(agent).exploration_map(x_index+k, y_index+j) + 1 - 2*distance < min_pot
                min_pot = agents(agent).exploration_map(x_index+k, y_index+j) + 1 - 2*distance;
                min_pot_xi = x_index+k;
                min_pot_yi = y_index+j;
                tp_list(end+1, :) = [min_pot_xi, min_pot_yi];
            end

            if agents(agent).exploration_map(x_index+k, y_index+j) <= max_pot
                 if agents(agent).color == 1
                     agents(agent).exploration_map(x_index+k, y_index+j) = agents(agent).exploration_map(x_index+k, y_index+j) + increment;
                     if agents(agent).exploration_map(x_index+k, y_index+j) > 1
                          agents(agent).exploration_map(x_index+k, y_index+j) = 1;
                     end
                 end
            end

            if size(tp_list, 1) < 3
                tp_list(end+1, :) = [x_index+k, y_index+j];
            end

        end

        j = -radie;
        for k = radie - 1: -1 : -radie + 1
                            
            %undviker krash
            if (x_index+k <= 1) || (x_index+k >= length(M)) || (y_index+j <= 1) || (y_index+j >= length(M))
                continue
            end
            
            distance = 0;

            if agents(agent).bresenham
                for pursuer = 1:length(agents(agent).px_list)
                    distance = distance + norm([agents(agent).px_list(pursuer), agents(agent).py_list(pursuer)] - [(x_index+k)/N_matrix, (y_index+j)/N_matrix])/sqrt(2);
                end
                distance = distance/length(agents(agent).px_list); % / (2*);
            end
            
            if agents(agent).exploration_map(x_index+k, y_index+j) + 1 - 2*distance < min_pot
                min_pot = agents(agent).exploration_map(x_index+k, y_index+j) + 1 - 2*distance;
                min_pot_xi = x_index+k;
                min_pot_yi = y_index+j;
                tp_list(end+1, :) = [min_pot_xi, min_pot_yi];
            end

            if agents(agent).exploration_map(x_index+k, y_index+j) <= max_pot
                 if agents(agent).color == 1
                     agents(agent).exploration_map(x_index+k, y_index+j) = agents(agent).exploration_map(x_index+k, y_index+j) + increment;
                      if agents(agent).exploration_map(x_index+k, y_index+j) > 1
                          agents(agent).exploration_map(x_index+k, y_index+j) = 1;
                      end
                 end
            end

        end

        if size(tp_list, 1) < 4
            tp_list(end+1, :) = [x_index+k, y_index+j];
        end

        radie = radie + 1;
    end

    if jesus_boolean

        x0 = agents(agent).x;
        y0 = agents(agent).y;
        motion_angle = agents(agent).theta;

        m = tan(motion_angle);

        k = y0 - m * x0;
        str_line_y = @(x) m*x + k;
        str_line_x = @(y) (y - k) / m;

        % finn vilka två sidor på randen vi letar efter
        if motion_angle >= 0 && motion_angle < pi/2
            phi = atan((1-y0) / (1- x0));
            if abs(motion_angle) < abs(phi)
                x = 1;
                y = str_line_y(1);
            else
                y = 1;
                x = str_line_x(1);
            end
        
        elseif motion_angle >= pi/2 && motion_angle <= pi
            phi = pi - atan((1-y0) / x0);
            if abs(motion_angle) < abs(phi)
                y = 1;
                x = str_line_x(1);
            else
                x = 0;
                y = str_line_y(0);
            end
        
        elseif motion_angle < 0 && motion_angle > -pi/2
            phi = atan(y0/ (1-x0));
            if abs(motion_angle) < abs(phi)
                x = 1;
                y = str_line_y(1);
            else
                y = 0;
                x = str_line_x(0);
            end
        else
            phi = pi - atan(y0 / x0);
            if abs(motion_angle) < abs(phi)
                y = 0;
                x = str_line_x(0);
            else
                x = 0;
                y = str_line_y(0);
            end
        end

        if agents(agent).color == 2
            % px_list has all pursuers
            % for each pursuer, make a set of invalid angles
            invalid_angles = [];
            for pursuer = 1:length(agents(agent).px_list)
                % Calculate the angle between the agent's current heading and the direction to the target
                angle_to_pursuer = atan2(agents(agent).py_list(pursuer) - agents(agent).y, agents(agent).px_list(pursuer)- agents(agent).x);
                invalid_angles(end+1) = angle_to_pursuer;
            end
            for angle_index = 1:length(invalid_angles)
                if agents(agent).theta <= invalid_angles(angle_index) + pi/8 && agents(agent).theta >= invalid_angles(angle_index) - pi/8
                    % agents(agent).temp = [];
                    agents(agent).tx = min_pot_xi / N_matrix;
                    agents(agent).ty = min_pot_yi / N_matrix;
                    return
                end
            end
        end

        cells = bresenham_cells(ceil(agents(agent).x*N_matrix), ceil(agents(agent).y*N_matrix), xy_index_check(ceil(x*N_matrix), N_matrix), xy_index_check(ceil(y*N_matrix), N_matrix), fieldmap);
        temp_min_pot = inf;

        if ~isempty(cells)
            for cell = 1:length(cells(:, 1))
                distance = 0;

                if agents(agent).color == 1 && ~agents(agent).bresenham && fieldmap(cells(cell, 1), cells(cell, 2)) ~= 1
                    distance = jesus_scale_factor*norm([agents(agent).x - cells(cell, 1)/N_matrix; agents(agent).y - cells(cell, 2)/N_matrix]);
                end

                if agents(agent).bresenham
                    for pursuer = 1:length(agents(agent).px_list)
                        distance = distance + norm([agents(agent).px_list(pursuer), agents(agent).py_list(pursuer)] - [cells(cell, 1)/N_matrix, cells(cell, 2)/N_matrix])/sqrt(2);
                    end
                    distance = distance/length(agents(agent).px_list); % / (2*);
                end

                if agents(agent).exploration_map(cells(cell,1), cells(cell, 2)) + 1 - 2*distance < min_pot
                    temp_min_pot = agents(agent).exploration_map(cells(cell,1), cells(cell, 2)) + 1 - 2*distance;
                    temp_min_pot_xi = cells(cell,1);
                    temp_min_pot_yi = cells(cell, 2);
                end
            end

            if temp_min_pot < min_pot
                if bresenham(ceil(agents(agent).x*N_matrix), ceil(agents(agent).y*N_matrix), temp_min_pot_xi, temp_min_pot_yi, fieldmap, 1)
                    min_pot = temp_min_pot;
                    min_pot_xi = temp_min_pot_xi;
                    min_pot_yi = temp_min_pot_yi;
                    tp_list(end+1, :) = [min_pot_xi, min_pot_yi];
                    if fieldmap(temp_min_pot_xi, temp_min_pot_yi)<0.01
                        agents(agent).change_point = 1;
                    end
                end
            end
        end
      
    end

    if min_pot >= 2
        x_coord = ceil(agents(agent).x*N_matrix);
        y_coord = ceil(agents(agent).y*N_matrix);
        start_x = xy_index_check(x_coord - 25, N_matrix);
        end_x = xy_index_check(x_coord + 25, N_matrix);
        start_y = xy_index_check(y_coord - 25, N_matrix);
        end_y = xy_index_check(y_coord + 25, N_matrix);
        agents(agent).exploration_map(start_x:end_x, start_y:end_y) = fieldmap(start_x:end_x, start_y:end_y);
        
        % disp("partially erase explotation map")
        [agents, M] = set_target_points(agents, max_radie, x_index, y_index, M, agent, increment, max_pot, N_matrix, fieldmap, kill_radius, sighting_target_boolean, jesus_boolean, smart_target_boolean, h, iterations_ahead, increment_boolean, factor_reach_max_angular_velocity, max_angular_velocity, max_iter, crystal_tp_ahead, pursuit_crystal, increment_scale_factor, increment_decay_rate, jesus_scale_factor);
        return
    end 

    agents(agent).tx = min_pot_xi / N_matrix;
    agents(agent).ty = min_pot_yi / N_matrix;

    tp_list = [(x_index+4), (y_index); 
    (x_index), (y_index+4); 
    (x_index-4), (y_index); 
    (x_index), (y_index-4); 
    (x_index+4), (y_index+4); 
    (x_index+4), (y_index-4);
    (x_index-4), (y_index+4); 
    (x_index-4), (y_index-4);
    tp_list];

    % check list
    agents = crystal_tp(agents, agent, h, fieldmap, N_matrix, tp_list, factor_reach_max_angular_velocity, max_angular_velocity, kill_radius, max_iter, crystal_tp_ahead);

end
