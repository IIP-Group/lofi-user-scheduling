% =========================================================================
% -- Gradient descent-based scheduler
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

function [par,var,chosen_C_final,chosen_C_final_quantized,best_result,best_initialization,status_C,f_vec] = optimization_based_scheduling(par,var)

var = param_optimization(var,par);

% Initializing matrices for all random initializations
C = zeros(par.U,par.timeslots,var.grad_iterations,var.num_initializations); % initializing scheduling matrix
C_quantized = zeros(par.U,par.timeslots,var.grad_iterations,var.num_initializations); % initializing quantized scheduling matrix
C_final = zeros(par.U,par.timeslots,var.num_initializations); % C after var.grad_iterations of gradient descent - all random initilizations
C_final_quantized = zeros(par.U,par.timeslots,var.num_initializations); % quantized final C - all random initilizations
f_with_regularizer = zeros(var.grad_iterations,var.num_initializations); % objective function with regularization
f_without_regularizer = zeros(var.grad_iterations,var.num_initializations); % objective function without regularization
f_with_regularizer_quant = zeros(var.grad_iterations,var.num_initializations); % objective function with regularization
f_without_regularizer_quant = zeros(var.grad_iterations,var.num_initializations); % objective function without regularization
f_with_regularizer_quant_final = zeros(var.num_initializations,1); % final value of objective function - all random initilizations
f_without_regularizer_quant_final = zeros(var.num_initializations,1); % final value of objective function - all random initilizations

for l=1:var.num_initializations % For different random initilizations
    C(:,:,1,l) = rand(par.U,par.timeslots); % random initialization of scheduling matrix
    [C(:,:,1,l)] = ProjMatrix(par,C(:,:,1,l)); % projection of C on simplexes
    [~,C_quantized(:,:,1,l)] = calculate_quantization(par,C(:,:,1,l));

    % Calculate objective function for current C
    var.calc_grad = 'yes';
    [f_with_regularizer(1,l),f_without_regularizer(1,l),regularization,...
        gradient_f_with_regularizer,~,~] = calculate_objective_and_gradient(par,var,C(:,:,1,l));
    var.calc_grad = 'no';
    [f_with_regularizer_quant(1,l),f_without_regularizer_quant(1,l),~,~,~,~] = calculate_objective_and_gradient(par,var,C_quantized(:,:,1,l));

    for i=1:var.grad_iterations
        fprintf('l : %d, i : %d\n', l, i);
        %% Computing the gradient of the objective function

        % Gradient descent
        Z = C(:,:,i,l) - var.tau.*gradient_f_with_regularizer; % Minimization

        % Projection on Simplexes
        [C(:,:,i+1,l)] = ProjMatrix(par,Z);
        %C_quantized(:,:,i+1,l) = (C(:,:,i+1,l)>= 0.5);
        [~,C_quantized(:,:,i+1,l)] = calculate_quantization(par,C(:,:,i+1,l));

        % Calculate objective function for current C
        var.calc_grad = 'yes';
        [f_with_regularizer(i+1,l),f_without_regularizer(i+1,l),regularization,...
            gradient_f_with_regularizer,~,~] = calculate_objective_and_gradient(par,var,C(:,:,i+1,l));
        var.calc_grad = 'no';
        [f_with_regularizer_quant(i+1,l),f_without_regularizer_quant(i+1,l),~,~,~,~] = calculate_objective_and_gradient(par,var,C_quantized(:,:,i+1,l));

    end
    C_final(:,:,l) = C(:,:,i+1,l); % Final matrix C for different random initializations
    [status_C,C_final_quantized(:,:,l)] = calculate_quantization(par,C_final(:,:,l));% Adjusting final matrix C to 0's and 1's
    % Calculate objective function for final C's
    var.calc_grad = 'no';
    [f_with_regularizer_quant_final(l),f_without_regularizer_quant_final(l),~,~,~,~] = calculate_objective_and_gradient(par,var,C_final_quantized(:,:,l));
end

best_result = inf;

for l=1:var.num_initializations
    if(f_without_regularizer_quant_final(l) < best_result)
        best_result = f_without_regularizer_quant_final(l);
        best_initialization = l;
    end
end

% Choose the random initilization that generates the smallest objective function
chosen_C_final = C_final(:,:,best_initialization);
chosen_C_final_quantized = C_final_quantized(:,:,best_initialization);

%% Plotting MSE vs. iterations

if (strcmp(par.visualize_opt_graphs,'on'))
    figure(1)
    plot(1:length(f_without_regularizer(:,best_initialization)),real(f_without_regularizer(:,best_initialization)),'-b','LineWidth',2);
    hold on;
    plot(1:length(f_with_regularizer(:,best_initialization)),real(f_with_regularizer(:,best_initialization)),'-r','LineWidth',2);
    legend('Real Cost Function','Total Cost Function');
    xlabel('Iterations');
    ylabel('MSE');
    title('Without quantization');

    figure(2)
    plot(1:length(f_without_regularizer_quant(:,best_initialization)),real(f_without_regularizer_quant(:,best_initialization)),'-b','LineWidth',2);
    hold on;
    plot(1:length(f_with_regularizer_quant(:,best_initialization)),real(f_with_regularizer_quant(:,best_initialization)),'-r','LineWidth',2);
    legend('Real Cost Function','Total Cost Function');
    xlabel('Iterations');
    ylabel('MSE');
    title('With quantization');
end

f_vec = real(f_without_regularizer(:,best_initialization)); % keep to take the mean out of many iterations
%toc;
end