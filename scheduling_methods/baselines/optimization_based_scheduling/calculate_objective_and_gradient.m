% =========================================================================
% -- Optimization-based Scheduling
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

function [f_with_regularizer,f_without_regularizer,regularization,gradient_f_with_regularizer,...
    grad_f_without_regularizer,grad_regularization] = calculate_objective_and_gradient(par,var,C)

% Regularizer
[regularization,grad_regularization] = L2_norm(par,var,C);

% Objective function
[f_without_regularizer,grad_f_without_regularizer] = MMSE(par,var,C);

f_with_regularizer = f_without_regularizer - regularization;
gradient_f_with_regularizer = grad_f_without_regularizer - grad_regularization;


end
