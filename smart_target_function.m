function agents = smart_target_function(agents,agent, h, fieldmap, N_matrix, iterations_ahead)
    evader_index = agents(agent).evader_index;
    x_init = agents(evader_index).x;
    y_init = agents(evader_index).y;
    vx_init = agents(evader_index).vx;
    vy_init = agents(evader_index).vy;
    theta_init = agents(evader_index).theta;
    % euler forwards
    for iter = 1:iterations_ahead
        x_next = x_init + vx_init * h;
        y_next = y_init + vy_init * h;

        % wall checker
        % out of bounds
        if x_next <= 0 || x_next >= 1 || y_next <= 0 || y_next >= 1
            break
        elseif fieldmap(ceil(x_next * N_matrix), ceil(y_next * N_matrix)) == 1
            break
        elseif ~bresenham(ceil(agents(agent).x*N_matrix), ceil(agents(agent).y*N_matrix), ceil(x_next * N_matrix), ceil(y_next * N_matrix), fieldmap, 1)
            break
        end

        x_init = x_next;
        y_init = y_next;
        theta_init = theta_init + agents(evader_index).omega * h;
        
        v = sqrt(agents(evader_index).vx^2 + agents(evader_index).vy^2);
        % composant projection
        vx_init = v * cos(theta_init);
        vy_init = v * sin(theta_init);
    end

    agents(agent).tx = x_init;
    agents(agent).ty = y_init;

end

