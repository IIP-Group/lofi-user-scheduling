% =========================================================================
% -- Random scheduling
% -------------------------------------------------------------------------
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
function [C_matrix,result] = random_scheduling(par,var)

C_matrix = zeros(par.U,par.timeslots);

% Assigning first Us UEs to first timeslot and the remaining to the second
for i=1:par.timeslots
    for j=1:par.Us
        C_matrix((i-1)*par.Us+j,i) = 1;
    end
end

% Computing the MSE of this scheduling
MSE = zeros(par.timeslots,1);
for t=1:par.timeslots
    C = C_matrix(:,t);
    Hc = var.H*diag(C);
    L = pinv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));
    M = Hc*diag(C.^2)*Hc';
    MSE(t) = var.Es*trace(diag(C.^2))-var.Es*trace(L*M);
end
result = real(sum(MSE));

end






