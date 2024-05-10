function agents = update_bresenham(agents, fieldmap, N, N_matrix)
    % empty line-of-sight vectors
    for agent = 1:N
        agents(agent).px_list = [];
        agents(agent).py_list = [];
    end
        
    % check line of sight
    for agent = 1:N
        % Persuer
        if agents(agent).color == 1
        
            for agent2 = 1:N
                % Evader
                if agents(agent2).color == 2

%                     disp("update_bresenham")
                    % crash prevention
                    x0 = xy_index_check(ceil(agents(agent).x * N_matrix), N_matrix);
%                     disp("x0:" + num2str(x0))
                    x1 = xy_index_check(ceil(agents(agent2).x * N_matrix), N_matrix);
%                     disp("x1:" + num2str(x1))
                    y0 = xy_index_check(ceil(agents(agent).y * N_matrix), N_matrix);
%                     disp("y0:" + num2str(y0))
                    y1 = xy_index_check(ceil(agents(agent2).y * N_matrix), N_matrix);
%                     disp("y1:" + num2str(y1))
%                     disp("end of update_bresenham")
                    
                    if bresenham(x0, y0, x1, y1, fieldmap, 1)
                        agents(agent).px = [];
                        agents(agent).py = [];
                        agents(agent).exploration_map = fieldmap;
                        % each pursuer has a list of evaders within
                        % sight's positions
                        agents(agent).px_list(end+1) = agents(agent2).x;
                        agents(agent).py_list(end+1) = agents(agent2).y;
                        agents(agent).bresenham = true;
                        agents(agent2).bresenham = true;
                        % each evader has a list of pursuers within
                        % sight's positions
                        agents(agent2).px_list(end+1) = agents(agent).x;
                        agents(agent2).py_list(end+1) = agents(agent).y;
                    end
        
                end
            end
        else
            %
        end
    end
       
    % iteratate and set breshham to false if px_list is empty
    for agent = 1:N
        if isempty(agents(agent).px_list)
            agents(agent).bresenham = false;
        end
    end
end