% =========================================================================
% -- Greedy User Scheduling
% -------------------------------------------------------------------------
% Reference:
%         J. Choi, G. Lee and B. L. Evans, "User Scheduling for
%         Millimeter Wave Hybrid Beamforming Systems With Low-Resolution
%         ADCs," in IEEE Transactions on Wireless Communications, vol. 18,
%         no. 4, pp. 2401-2414, April 2019, doi: 10.1109/TWC.2019.2904030.
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function [C_matrix] = greedy(par,var)

% Obtaining H_b
A_0 = dftmtx(par.B); % DFT matrix
A = zeros(par.B); % normalizing each column of F to 1
for i=1:par.B
    A(:,i) = A_0(:,i)/norm(A_0(:,i));
end
H_b = A*var.H; % doing the DFT of the channel matrix H

% Step 1
U_group = 1:par.U;
H_gre = [];
order = [];
i = 1;

while (i<=par.Us)
    % Step 2
    sum_rate = zeros(1,length(U_group));
    for j=1:length(U_group)
        H_test = [H_gre H_b(:,U_group(j))];
        [sum_rate(j)] = calculate_sum_rate(par,var,H_test);
    end
    [~,sel] = max(sum_rate);
    order(i) = U_group(sel);

    % Step 3
    H_gre = [H_gre H_b(:,U_group(sel))];

    % Update U_group
    U_group(sel) = [];

    i=i+1;
end

% Getting the correct order
order_asc = sort(order,'ascend');

C_matrix = zeros(par.U,par.timeslots);

C_matrix(order_asc,1) = 1;

C_matrix(:,2) = (C_matrix(:,1) == 0);

end
