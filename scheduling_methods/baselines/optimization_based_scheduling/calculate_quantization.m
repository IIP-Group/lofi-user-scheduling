% =========================================================================
% -- Quantization function
% -------------------------------------------------------------------------
% Reference:
%           V. Palhares and C. Studer, "An Optimization-Based User 
%           Scheduling Framework for mmWave Massive MU-MIMO Systems," 2022 
%           IEEE 23rd International Workshop on Signal Processing Advances 
%           in Wireless Communication (SPAWC), Oulu, Finland, 2022, 
%           pp. 1-5, doi: 10.1109/SPAWC51304.2022.9833954.
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function [status,C_quant] = calculate_quantization(par,C)

% Simple quantization
C_2 = (C>= 0.5);

[status,constraints] = satisfy_constraints(par,C_2);

% Second attempt
while(status ~= 1)
    if(constraints(3) ~= 1) % Columns are not correct
        users_per_t = sum(C_2,1);
        ind_2 = find(users_per_t < par.min_sum_columns); % how many columns need more users
        ind_3 = find(users_per_t > par.max_sum_columns); % how many columns need less users
        
        % Columns that need more users
        closer_to_threshold_min = [];
        for u=1:length(ind_2)
            changed_users = par.min_sum_columns - users_per_t(ind_2(u)); % how many need to be scheduled
            ind_4 = find(C_2(:,ind_2(u)) == 0); % find the ones set to 0
            [~,closer_to_threshold_min(:,u)] = maxk(C(ind_4,ind_2(u)),changed_users); % get the closer ones from the threshold
            C_2(ind_4(closer_to_threshold_min),ind_2(u)) = 1; % Set 1's to reach min sum of columns
        end
        
        % Columns that need less users
        closer_to_threshold_max = [];
        for u=1:length(ind_3)
            changed_users = users_per_t(ind_3(u)) - par.max_sum_columns; % how many need to be removed
            ind_4 = find(C_2(:,ind_3(u)) == 1); % find the ones set to 1
            [~,closer_to_threshold_max(:,u)] = mink(C(ind_4,ind_3(u)),changed_users); % get the closer ones from the threshold
            C_2(ind_4(closer_to_threshold_max),ind_3(u)) = 0; % Set 1's to reach max sum of columns
        end
    end
    
    % Check if row constraint is satisfied now
    [status,constraints] = satisfy_constraints(par,C_2);
    
    if(constraints(2) ~= 1) % Rows are not correct
        users_per_t = sum(C_2,1);
        t_per_user = sum(C_2,2);
        if(~isempty(ind_2)) % Means that we scheduled users
            timeslots = closer_to_threshold_min(:); % timeslots that have changed
            for u=1:length(timeslots)
                changed_timeslots = t_per_user(timeslots(u)) - par.max_sum_rows; % how many need to be removed
                for uu = 1:changed_timeslots
                    ind_4 = find(C_2(timeslots(u),:) == 1); % See what other timeslots this user is in
                    [~,ind_5] = max(users_per_t(ind_4)); % get the one with the largest number of users already selected
                    C_2(timeslots(u),ind_5) = 0; % Set 0's to reach max sum of rows
                    users_per_t = sum(C_2,1);
                end
            end
        end
        if(~isempty(ind_3)) % Means that we removed users
            timeslots = closer_to_threshold_max(:); % timeslots that have changed
            for u=1:length(timeslots)
                changed_timeslots = t_per_user(timeslots(u)) - par.max_sum_rows; % how many need to be removed
                for uu = 1:changed_timeslots
                    ind_4 = find(C_2(timeslots(u),:) == 1); % See what other timeslots this user is in
                    [~,ind_5] = min(users_per_t(ind_4)); % get the one with the largest number of users already selected
                    C_2(timeslots(u),ind_5) = 1; % Set 0's to reach max sum of rows
                    users_per_t = sum(C_2,1);
                end
            end
        end
    end
    
    % Check if column constraint is satisfied now
    [status,constraints] = satisfy_constraints(par,C_2);
end

C_quant = C_2;


end





