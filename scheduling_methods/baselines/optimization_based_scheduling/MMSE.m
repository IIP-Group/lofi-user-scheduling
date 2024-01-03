% =========================================================================
% -- MMSE objective function and gradient
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

function [f_without_regularizer,grad_f_without_regularizer] = MMSE(par,var,C)

f_without_regularizer = 0;
grad_f_without_regularizer = zeros(par.U,par.timeslots);
for t = 1:par.timeslots
    Hc = var.H*diag(C(:,t));
    L = pinv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));
    M = Hc*diag(C(:,t).^2)*Hc';

    f_without_regularizer = f_without_regularizer+var.Es*trace(diag(C(:,t).^2))-var.Es*trace(L*M);

    % Gradient calculation
    if (strcmp(var.calc_grad,'yes'))
        for u=1:par.U
            L_H = L*var.H(:,u);
            H_L = var.H(:,u)'*L;
            H_L_M = H_L*M;
            % gradient of the objective function without regularizer
            grad_f_without_regularizer(u,t) = var.Es*(2*C(u,t)+2*C(u,t)*H_L_M*L_H-4*C(u,t)^3*var.H(:,u)'*L_H);
        end
    end

end

end

