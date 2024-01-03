% =========================================================================
% -- Exhaustive search with MMSE criterion
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

function [best_C,best_measure] = exhaustive_search(par,var)

combos = nchoosek(1:par.U,par.Us);
vec_k = 1:par.U;
[var.num_poss,~] = size(combos);
C_matrix = zeros(par.U,par.timeslots,var.num_poss);

for i=1:par.U
    for j=1:var.num_poss
        C_matrix(i,1,j) = ismember(vec_k(i),combos(j,:));
        C_matrix(i,2,j) = (C_matrix(i,1,j) == 0);
    end
end

measure = zeros(par.timeslots,var.num_poss);
sum_measure = zeros(var.num_poss,1);

for n=1:var.num_poss
    for t=1:par.timeslots
        C = C_matrix(:,t,n);
        Hc = var.H*diag(C);
        % MSE Calculation
        L = pinv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));
        M = Hc*diag(C.^2)*Hc';
        measure(t,n) = var.Es*trace(diag(C.^2))-var.Es*trace(L*M);
    end
    sum_measure(n) = sum(measure(:,n));
end


best_measure = min(sum_measure);
best_poss = find(sum_measure == best_measure);
best_C = C_matrix(:,:,best_poss);

end












