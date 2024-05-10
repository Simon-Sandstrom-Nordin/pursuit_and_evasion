function agents = sighting_target_function(agents, agent, kill_radius, fieldmap, max_iter, tp_crystal_iter_ahead, max_angular_velocity, factor_reach_max_angular_velocity, h)
len_field = length(fieldmap);
% if length(agents(agent).px)>1
%      disp([302.*agents(agent).px(end-4:end), 302.*agents(agent).py(end-4:end)])
% end
% disp("using sighting")
prevpx = agents(agent).px;
prevpy = agents(agent).py;
v = sqrt(agents(agent).vx^2 + agents(agent).vy^2);

% if ~iter_ahead_function(1, v, agents(agent).x, agents(agent).y, max_iter, tp_crystal_iter_ahead, agents(agent).px(1), agents(agent).py(1), agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap)
%     agents(agent).tx = agents(agent).px(1);
%     agents(agent).ty = agents(agent).py(1);
%     agents(agent).px = agents(agent).px(1);
%     agents(agent).py = agents(agent).py(1);
%     return
% end

if length(agents(agent).px) == 1
    
    if ~iter_ahead_function(1, v, agents(agent).x, agents(agent).y, max_iter, tp_crystal_iter_ahead, agents(agent).px, agents(agent).py, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap) && bresenham(ceil(len_field*agents(agent).x),ceil(len_field*agents(agent).y), ceil(len_field*agents(agent).px),ceil(len_field*agents(agent).py), fieldmap,0)
        %disp("using one")
        agents(agent).tx = agents(agent).px;
        agents(agent).ty = agents(agent).py;
        agents(agent).px(end+1) = agents(agent).x;
        agents(agent).py(end+1) = agents(agent).y;
        return

    else
        %disp("using two")
        agents(agent).px(end+1) = agents(agent).x;
        agents(agent).py(end+1) = agents(agent).y;
        agents(agent).tx = [];
        agents(agent).ty = [];
        return

    end

else

    v = sqrt(agents(agent).vx^2 + agents(agent).vy^2);
    if ~iter_ahead_function(1, v, agents(agent).x, agents(agent).y, max_iter, tp_crystal_iter_ahead, agents(agent).px(end-1), agents(agent).py(end-1), agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap) && bresenham(ceil(len_field*agents(agent).x),ceil(len_field*agents(agent).y), ceil(len_field*agents(agent).px(end-1)),ceil(len_field*agents(agent).py(end-1)), fieldmap,0)
        %disp("using 1")
        agents(agent).tx = agents(agent).px(end-1);
        agents(agent).ty = agents(agent).py(end-1);
        return
    else
        %disp("using 2")
        if ~iter_ahead_function(1, v, agents(agent).x, agents(agent).y, max_iter, tp_crystal_iter_ahead, agents(agent).px(end), agents(agent).py(end), agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap)
            agents(agent).tx = agents(agent).px(end);
            agents(agent).ty = agents(agent).py(end);
            return
        else
            %disp("using three")
            agents(agent).px(end+1) = agents(agent).x;
            agents(agent).py(end+1) = agents(agent).y;
            agents(agent).tx = [];
            agents(agent).ty = [];
            return
        end
    end
end

% agents(agent).px(end) = agents(agent).x;
% agents(agent).py(end) = agents(agent).y;

x_difference = abs(agents(agent).x - agents(agent).px);
y_difference = abs(agents(agent).y - agents(agent).py);

% disp("if statement 3")
if ~isempty(agents(agent).px)
    
    
    while norm([x_difference(end-1), y_difference(end-1)]) < 10*kill_radius && bresenham(ceil(len_field*agents(agent).x),ceil(len_field*agents(agent).y), ceil(len_field*agents(agent).px(end-1)),ceil(len_field*agents(agent).py(end-1)), fieldmap,0)
        agents(agent).px(end) = [];
        x_difference(end) = [];
        agents(agent).py(end) = [];
        y_difference(end) = [];
        if length(agents(agent).px) == 1
            %disp("using 7")
            agents(agent).tx = agents(agent).px;
            agents(agent).ty = agents(agent).py;
            agents(agent).px = [];
            agents(agent).py = [];   
            return
        end
    end 
    %disp("using 8")
    agents(agent).tx = agents(agent).px(end-1);
    agents(agent).ty = agents(agent).py(end-1);
   
    v = sqrt(agents(agent).vx^2 + agents(agent).vy^2);
    if iter_ahead_function(1, v, agents(agent).x, agents(agent).y, max_iter, tp_crystal_iter_ahead, agents(agent).tx, agents(agent).ty, agents, agent, max_angular_velocity, factor_reach_max_angular_velocity, h, kill_radius, fieldmap)
        %disp("using 9")
        agents(agent).tx = [];
        agents(agent).ty = [];
        return
    end

    
    % if isempty(agents(agent).tx) && sempty(agents(agent).tx)
    %     % disp("ingen tx eller ty")
    % else 
    %     % disp("vi har tx eller ty")
    % end
end
end
