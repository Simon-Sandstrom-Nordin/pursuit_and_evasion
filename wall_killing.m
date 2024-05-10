function [agents, wall_pursuer_counter, wall_evader_counter] = wall_killing(agents, N_matrix, base_M)

    wall_pursuer_counter = 0;
    wall_evader_counter = 0;
    % Game logic: agents within walls die
    to_be_killed = [];
    for agent = 1:length(agents)
        x_index = xy_index_check(ceil(agents(agent).x * N_matrix), N_matrix);
        y_index = xy_index_check(ceil(agents(agent).y * N_matrix), N_matrix);
        if base_M(x_index, y_index) == 3
            to_be_killed(end+1) = agent;
        end
    end

    for agent = length(to_be_killed):-1:1
        if agents(to_be_killed(agent)).color == 1
            wall_pursuer_counter = wall_pursuer_counter + 1;
        else
            wall_evader_counter = wall_evader_counter + 1;
        end
       % disp("wall_killing ----------------------------------------------")
        agents(to_be_killed(agent)) = [];
        for k = 1:length(agents)    % fula sättet som inte generaliserar till fler pursuers, syftet är för sighting target så att den inte följer in i väggen efter död evader
            if agents(k).color == 1
                agents(k).px = [];
            end
        end
    end
end
