function M = random_walk_map_generator(obstacle_density, N_matrix, column_width)
% restrict column_width to divisors to N_matris

obstacle_density = 1 - obstacle_density;

% scale down by factor column_width
M = ones(ceil((N_matrix-2) / column_width));  % inner portion only
number_of_blanc_pixles = obstacle_density * length(M)^2 - 1;

x = ceil(rand(1)*length(M));
y = ceil(rand(1)*length(M));
M(x, y) = 0;
% random walk
previous_direction = 0;
while number_of_blanc_pixles > 0
    % inital random walk
    direction = randi([1, 4]);
    
    if direction == 1 && previous_direction ~= 3   % up
        if y + 1 > length(M)
            continue
        end
        y = y + 1;
    elseif direction == 2 && previous_direction ~= 4   % right
        if x + 1 > length(M)
            continue
        end
        x = x + 1;
    elseif direction == 3 && previous_direction ~= 1   % down
        if y - 1 < 1
            continue
        end
        y = y - 1;
    elseif direction == 4 && previous_direction ~= 2   % left
        if x - 1 < 1
            continue
        end
        x = x - 1;
    end
    
    if M(x, y) == 1
        M(x, y) = 0;
        number_of_blanc_pixles = number_of_blanc_pixles - 1;
    end

    previous_direction = direction;
end

% scale up
scaled_up_M = kron(M, ones(column_width));

% add walls
ones_only = ones(N_matrix);
row_indices = 2:N_matrix-1;
col_indices = 2:N_matrix-1;
ones_only(row_indices, col_indices) = scaled_up_M;
M = 3*ones_only;
end
