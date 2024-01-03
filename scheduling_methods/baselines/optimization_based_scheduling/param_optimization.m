% =========================================================================
% -- Defining optimization parameters
% -------------------------------------------------------------------------
%
% Description:
% Defines regularizer and gradient descent stepsize for "Opt.-based" 
% scheduler
%
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

function var = param_optimization(var,par)

var.num_initializations = 80; % # of different random initializations
var.grad_iterations = 100; % # of iterations of the gradient descent

switch(par.sim_scenario)
    case 'a' % B = 16, U = 16
        var.lambda = 0.01; % regularizer
        var.tau = 1*10^(-1); % gradient descent stepsize
    otherwise
        error('Variable par.sim_scenario is not defined.');
end

end