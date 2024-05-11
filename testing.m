clc; clear all; close all; format long;

N_matrix = 302;

M = zeros(N_matrix-2);

% scale up
% add walls
ones_only = ones(N_matrix);
row_indices = 2:N_matrix-1;
col_indices = 2:N_matrix-1;
ones_only(row_indices, col_indices) = M;
%M = 3*ones_only;

density_vec = 0:.1:0.9;
max_iter = 22;
crystaltp = 6;
pursuit_crystal = 0;
iterations_ahead = 3;

No_of_games = 1000;  % ger ungefär 1 % fluct
results_matrix = zeros(No_of_games+4, 1); % 1p vs 3e

% constant parameters
No_of_evaders = 2;  % välj kontrollvärde                      
No_of_pursuers = 1; % välj kontrollvärde

% booleans
visualization_boolean = false;
sighting_target_boolean = false;
increment_boolean = false;
smart_target_boolean = false;

% nya tester för increment of iter_ahead
% n_vec = [0];
percentage_vec = [];
% base model + _ + no.purs + _ + no.ev + _ + altalg
% columns are different target densities
% iterations until completion
% successful games procent
% median and average no. of iterations of these successful games
% [10, 12, 13, 12, 11, 10, procent, average, percentage_of_games_stopped_due_to_pursuer_death]';


% successful_games_vector = [0];
pursuer_death_vector = [];

% data matrices
pursuer_wall_deaths_matrix = zeros(1, No_of_games);
evader_wall_deaths_matrix = zeros(1, No_of_games);

jesus_scale_factor = 0;

increment = .3;
increment_scale_factor = 1;
increment_decay_rate = 2;

column_to_be_adjusted = 1;
win_counter = 0;
pursuer_death_counter = 0;
evader_wall_deaths_counter = 0;

for density = density_vec
    iterations_until_completion_matrix = [];
    for game = 1:No_of_games
        disp("Current game is no. " + num2str(game) + " so far " + num2str(win_counter) + " completed, or fractionally " + num2str(win_counter/game))
        
        [iterations_until_completion, reason_for_exit, pursuer_wall_deaths, evader_wall_deaths] = game_program(No_of_evaders, No_of_pursuers, density, visualization_boolean, ...
    sighting_target_boolean, true, smart_target_boolean, iterations_ahead, increment_boolean, increment, M, game, max_iter, crystaltp, pursuit_crystal, increment_scale_factor, increment_decay_rate);
        evader_wall_deaths_counter = evader_wall_deaths_counter + evader_wall_deaths;
        if (reason_for_exit == "All pursuers are dead")
            %break
            pursuer_death_counter = pursuer_death_counter + 1;
            iterations_until_completion = -1;
            iterations_until_completion_matrix(end+1) = iterations_until_completion;
            continue
        end
        if (reason_for_exit == "Maximum iterations reached")
            iterations_until_completion = -1;
            iterations_until_completion_matrix(end+1) = iterations_until_completion;

            continue
        end
        iterations_until_completion_matrix(end+1) = iterations_until_completion;
        disp("end of game" + num2str(game))
        win_counter = win_counter + 1;
    end
    
    no_of_failed_games = 0;
    no_of_games_ended_due_to_pursuer_death = 0;
    no_of_games_ended_due_to_max_iterations_reached = 0;
    
    mean_var = 0;
    for k = 1:length(iterations_until_completion_matrix)
        if iterations_until_completion_matrix(k) == -1
            no_of_failed_games = no_of_failed_games + 1;
        else
            mean_var = mean_var + iterations_until_completion_matrix(k);
        end
    end
    
    mean_var = mean_var / (No_of_games - no_of_failed_games);
    
    percentage = (No_of_games- no_of_failed_games) / No_of_games;

    results_matrix(:, column_to_be_adjusted) = [iterations_until_completion_matrix, percentage, mean_var, pursuer_death_counter, evader_wall_deaths_counter]';
    column_to_be_adjusted = column_to_be_adjusted + 1;
    name = "increment_10000_"+num2str(density)+".csv";
    writematrix(results_matrix, name)
end
