% =========================================================================
% -- Orthogonal projection onto simplex with inequality constraints
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
function [X] = Proximal_Algorithm(Z,a,b)

[size_Z] = size(Z);

X = zeros(size_Z);

for t=1:size_Z(2)

    %%%% Simplex Projection
    z = real(Z(:,t));
    n = length(z);
    z_sorted = sort(z,'descend');
    z_shift_lambda = z_sorted-z_sorted(1); % shift to 0

    %% Start linear search

    [cumulative_shifts,idx] = sort([-z_shift_lambda;(1-z_shift_lambda)]); % list all the shifts
    shifts = (idx > n); % know if you shifted to 0 or to 1 in that iteration

    difference = zeros(length(cumulative_shifts),1);
    difference(1) = -z_sorted(1)+cumulative_shifts(1);
    intervals = zeros(length(cumulative_shifts),2);
    intervals(1,:) = [-inf,-z_sorted(1)];

    x = z_shift_lambda;
    store_x = zeros(n,length(cumulative_shifts));
    store_x(:,1) = x;
    sum_x = zeros(length(cumulative_shifts),1);
    sum_x(1) = sum(x);
    final_MSE = inf;

    valid_lambda = [];

    m = sum(max(min(x,1),0)==1);
    p = sum(max(min(x,1),0)==0);

    % Find intervals
    for i=2:length(cumulative_shifts)
        current_valid_lambda = [];

        if(sum_x(i-1,1) > b)
            break;
        end
        % Finding the intervals where the # of 0's and 1's remain constant
        difference(i) = -z_sorted(1)+cumulative_shifts(i); % absolute lambda
        intervals(i,:) = [difference(i-1),difference(i)]; % intervals where the # of 0's and 1's remain constant

        x = z+difference(i); % calculate x using lambda = upper limit of the interval
        store_x(:,i) = x;
        sum_x(i,1) = sum(max(min(x,1),0)); % calculate the sum of x

        if(shifts(i-1) == 1)
            m = m+1; % get the # of 1's
        end
        if(shifts(i-1) == 0)
            p = p-1; % get the # of 0's
        end

        if(sum_x(i) >= a)
            endpoints = [intervals(i,1),intervals(i,2)]; % endpoints of the current interval

            z_j = sum(z_sorted(m+1:n-p)); % sum of the z's between 0 and 1
            current_lambdas = [(1/(n-p-m))*(a-m-z_j),(1/(n-p-m))*(b-m-z_j),0]; % possible lambdas
            if(sum_x(i-1) < a) % if the sum in the previous iteration was lower than the lower bound, have to update the interval of search
                endpoints(1) = current_lambdas(1); % new lower bound will be the lambda where sum(x) = a
            end
            if(sum_x(i) > b) % if the sum in the previous iteration was higher than the upper bound, have to update the interval of search
                endpoints(2) = current_lambdas(2); % new upper bound will be the one lambda where sum(x) = b
            end

            check_bound_1 = current_lambdas(1) >= endpoints(1) && current_lambdas(1) <= endpoints(2); % checking if first lambda is inside the interval
            check_bound_2 = current_lambdas(2) >= endpoints(1) && current_lambdas(2) <= endpoints(2); % checking if second lambda is inside the interval
            check_bound_3 = current_lambdas(3) >= endpoints(1) && current_lambdas(3) <= endpoints(2); % checking if third lambda is inside the interval

            if(check_bound_1 == 1 || check_bound_2 == 1 || check_bound_3 == 1) % if at least one of the two is inside the interval
                if(check_bound_1 == 1) % if the first lambda is inside the interval
                    current_valid_lambda = [current_valid_lambda;current_lambdas(1)]; % add lambda to the current list of valid possibilities
                end
                if(check_bound_2 == 1) % if the second lambda is inside the interval
                    current_valid_lambda = [current_valid_lambda;current_lambdas(2)]; % add lambda to the current list of valid possibilities
                end
                if(check_bound_3 == 1) % if the third lambda is inside the interval
                    current_valid_lambda = [current_valid_lambda;current_lambdas(3)]; % add lambda to the current list of valid possibilities
                end
            else % if both lambda are outside the interval, test the lower endpoint (upper endpoint will be tested in the next iteration)
                current_valid_lambda = [current_valid_lambda;intervals(i,1)];
            end
            valid_lambda = [valid_lambda;current_valid_lambda]; % store all valid lambdas from all iterations
            for j=1:length(current_valid_lambda) % for current valid lambdas
                possible_x = max(min(z+current_valid_lambda(j),1),0); % calculate x
                MSE = 1/2*(norm(possible_x-z))^2; % calculate the objective function for the given x
                if(MSE < final_MSE) % if the current MSE is lower than the smallest one so far
                    final_MSE = MSE; % update the MSE
                    final_lambda = current_valid_lambda(j); % update the lambda
                    X(:,t) = possible_x; % get the current best solution
                    final_m = m;
                    final_p = p;
                end
            end
        end
    end
end


end


