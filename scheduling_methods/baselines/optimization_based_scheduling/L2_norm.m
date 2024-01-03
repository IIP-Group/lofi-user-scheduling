% =========================================================================
% -- L2 regularizer
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

function [regularization,grad_regularization] = L2_norm(par,var,C)

regularization = 0;
grad_regularization = zeros(par.U,par.timeslots);
for t = 1:par.timeslots
    regularization = regularization+(var.lambda*norm(C(:,t)-(0.5.*ones(par.U,1)),'fro')^2);

    if (strcmp(var.calc_grad,'yes'))
        grad_regularization(:,t) = (var.lambda*(2*C(:,t)-ones(par.U,1)));
    end
end

end