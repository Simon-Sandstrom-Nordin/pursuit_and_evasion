function agents = pursuer_killing(agents, kill_radius)
    % Game logic: pursuer kill radius
    blue_agents_to_kill = [0];
    blue_agents_within_distance = [];
    for agent = 1:length(agents)
        if agents(agent).color == 1 && agents(agent).bresenham % Check if agent is red
            blue_agents = find([agents.color] == 2); % Find indices of blue agents
            distances = sqrt((agents(agent).x - [agents(blue_agents).x]).^2 + ...
                (agents(agent).y - [agents(blue_agents).y]).^2);
            blue_agents_within_distance = blue_agents(distances < kill_radius);
            if ~isempty(blue_agents_within_distance)
                
                agents(agent).px = [];
                agents(agent).py = [];
                agents(agent).bresenham = false;
                for blue_to_kill = 1: length(blue_agents_within_distance)
                    blue_agents_to_kill(end+1) = blue_agents_within_distance(blue_to_kill);
                end
            end
        end
    end

    blue_agents_to_kill = blue_agents_to_kill(2:end);
    blue_agents_to_kill = unique(blue_agents_to_kill(:).');
    blue_agents_to_kill = sort(blue_agents_to_kill, 'descend');
    % Check if there are blue agents to kill
    if ~isempty(blue_agents_to_kill)
        % Remove killed blue agents from the simulation
        for killed_agent = blue_agents_to_kill
            % disp("persuer_killing")
            agents(killed_agent) = [];
        end
    end
end
