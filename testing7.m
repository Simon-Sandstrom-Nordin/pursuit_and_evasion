clc; clear all; close all; format long;

N_matrix = 1;

%M = readmatrix('binary_image_alt_1.csv');
M = zeros(N_matrix-2);

% rand density, 0.2-0.7, 200 per spel, crystaltp 3-15

% scale up
% add walls
ones_only = ones(N_matrix);
row_indices = 2:N_matrix-1;
col_indices = 2:N_matrix-1;
ones_only(row_indices, col_indices) = M;
%M = 3*ones_only;

%density_vec = [.65]; % 
density_vec = .2:.1:0.8;
% crystal_tp_ahead = 5;
max_iter = 22;
%max_iter = 16;
% crystaltp = 7;
% max_iter_list = 2:2:6;

%crystaltp_list = 3:6;
crystaltp = 6;
%pursuit_crystal_list = 0:1;
pursuit_crystal = 0;
iterations_ahead = 3;

No_of_games = 1000;  % ger ungefär 1 % fluct
results_matrix = zeros(No_of_games + 4, 10); % 1p vs 3e

% evaders kör rakt in i väggar
% evaders target point skiftar fram och tillbaka med jesus
% evaders väljer inte punkt efter jesus
% 

% temp vars behövs inte annat är för testning
% behöver ej skriva ut varje gång längre

% iterations_ahead = 5;
iterations_ahead_list = 1:2;

% constant parameters
No_of_evaders = 2;  % välj kontrollvärde                      
No_of_pursuers = 1; % välj kontrollvärde
visualization_boolean = false;
sighting_target_boolean = true;
% sighting + increment
increment_boolean = true;

% increment = 0.05;
increment_list = .3;

%kolla korrekta incrementations värdet 
%incrementation_values = 0.1:0.1:1;

smart_target_boolean = true;

% results_matrix = zeros(No_of_games + 4, length(crystal_tp_ahead)); % 1p vs 3e

% incrementation_values = 0.6;
column_to_be_adjusted = 1;
iterations_until_completion_matrix = [];
for increment = increment_list
    % sighting target

    % nya tester för increment of iter_ahead
    % n_vec = [0];
    percentage_vec = [];
    % base model + _ + no.purs + _ + no.ev + _ + altalg
    % columns are different target densities
    % iterations until completion
    % successful games procent
    % median and average no. of iterations of these successful games
    % [10, 12, 13, 12, 11, 10, procent, average, percentage_of_games_stopped_due_to_pursuer_death]';

    % procent = ,10
   
    pursuer_death_counter = 0;
    evader_wall_deaths_counter = 0;
    
    % successful_games_vector = [0];
    pursuer_death_vector = [];
    
    % obstacle_density = .5;
    
    % data matrices
    iterations_until_completion_matrix = [];
    pursuer_wall_deaths_matrix = zeros(1, No_of_games);
    evader_wall_deaths_matrix = zeros(1, No_of_games);
    
    for game = 1:No_of_games

        obstacle = randi([2, 7])/10;
        disp("Max_iter is " + num2str(max_iter))
        disp("start of game" + num2str(game))
        % n_vec(end+1) = game;
        disp("Current game is no. " + num2str(game))
        [iterations_until_completion, reason_for_exit, pursuer_wall_deaths, evader_wall_deaths] = game_program(No_of_evaders, No_of_pursuers, rand(1), visualization_boolean, ...
sighting_target_boolean, true, smart_target_boolean, iterations_ahead, increment_boolean, increment, M, game, max_iter, crystaltp, pursuit_crystal);
        evader_wall_deaths_counter = evader_wall_deaths_counter + evader_wall_deaths;
        if (reason_for_exit == "All pursuers are dead")
            %break
            pursuer_death_counter = pursuer_death_counter + 1;
            iterations_until_completion = -1;
            iterations_until_completion_matrix(end+1) = iterations_until_completion;
            
            % successful_games_vector(end+1) = successful_games_vector(end);
            continue
        end
        if (reason_for_exit == "Maximum iterations reached")
            iterations_until_completion = -1;
            iterations_until_completion_matrix(end+1) = iterations_until_completion;
            % successful_games_vector(end+1) = successful_games_vector(end);
            continue
        end
        % successful_games_vector(end+1) = successful_games_vector(end) + 1;
        iterations_until_completion_matrix(end+1) = iterations_until_completion;
        % no_of_failed_games = 0;
        % for k = 1:length(iterations_until_completion_matrix)
        %     if iterations_until_completion_matrix(k) == -1
        %         no_of_failed_games = no_of_failed_games + 1;
        %     end
        % end
        % percentage_vec(end+1) = (game - no_of_failed_games) / game;
        disp("end of game" + num2str(game))
    end
    
    % figure(1)
    % pause(.00001)
    % plot(n_vec, successful_games_vector./n_vec, 'r-')
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
    %name = "basemodel" + num2str(incrementation_value) + "_200.csv";
    %writematrix(results_matrix,name)
    name = "SIC_1000.csv";
    writematrix(results_matrix, name)

end
% name = "tp_crystal_" + num2str(crystaltp) + "_200.csv";
% writematrix(results_matrix, name)
% name = "tp_crystal_" + num2str(incrementation_value) + "_density_" + obstacle + "_800.csv";
% writematrix(results_matrix, name)

tot_val = mean(results_matrix(end-2,:));
% disp("total result for increment value " + num2str(incrementation_value)+ " is " + num2str(tot_val))

%writematrix(results_matrix,"basemodel"+ "_200.csv")