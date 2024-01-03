% =========================================================================
% -- Channel Structure-based Scheduling (CSS)
% -------------------------------------------------------------------------
% Reference:
%         J. Choi, G. Lee and B. L. Evans, "User Scheduling for
%         Millimeter Wave Hybrid Beamforming Systems With Low-Resolution
%         ADCs," in IEEE Transactions on Wireless Communications, vol. 18,
%         no. 4, pp. 2401-2414, April 2019, doi: 10.1109/TWC.2019.2904030.
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria M. T. Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function [C_matrix] = CSS(par,var)

% Obtaining H_b
A_0 = dftmtx(par.B); % DFT matrix
A = zeros(par.B); % normalizing each column of F to 1
for i=1:par.B
    A(:,i) = A_0(:,i)/norm(A_0(:,i));
end
H_b = A*var.H; % doing the DFT of the channel matrix H

L = 4;
Nb = 2*L;
epsilon = 0.3;
Nol = 2;
Ns = par.Us;

% Step 2
[~, temp] = sort(H_b,'descend');
Bk = Nb+1 > temp; % Choosing the highest Nb beam indexes
idx = Bk == 1; % copying
Bk = rand(size(H_b));
Bk(idx) = 0;
Hbeam = [];

order_asc = 0;
while (length(order_asc) < Ns)
    Nol = Nol+1;
    epsilon = epsilon+0.1;

    % Step 1
    U_group = 1:par.U;
    H_css = [];
    f_0 = [];
    order = [];

    i = 1;

    while (i<=Ns && (sum(U_group) ~= 0))
        % Step 3
        R_u = zeros(1,length(U_group));
        SINR_u = zeros(1,length(U_group));
        f = zeros(par.B,length(U_group));
        for j=1:length(U_group)
            H_test = [H_css H_b(:,U_group(j))];
            [SINR_u(j)] = calculate_SINR(par,var,H_test);
        end
        [~,sel] = max(SINR_u);
        order(i) = U_group(sel);
        H_css = [H_css H_b(:,U_group(sel))];

        % Step 4
        % Compute norm orthogonal to span{g}
        temp = zeros(par.B);
        for j=1:i-1
            temp = temp + (f_0(:,j)*f_0(:,j)')/norm(f_0(:,j))^2;
        end
        f(:,sel) = (eye(par.B)-temp) * H_b(:,U_group(sel));
        f_0 = [f_0 f(:,sel)];

        % Step 5
        Hbeam = [Hbeam, Bk(:,U_group(sel))]; % save beam index for selected user

        % Update U_group
        U_group(sel) = [];
        temp = [];
        for j = 1:length(U_group)
            if abs(f_0(:,i)'*H_b(:,U_group(j)))/(norm(f_0(:,i))*norm(H_b(:,U_group(j)))) < epsilon
                if sum(Hbeam(:,i) == Bk(:,U_group(j)) ) <= Nol % represents the intersection
                    temp = [temp U_group(j)];
                end
            end
        end
        U_group = temp;

        % Step 6
        i=i+1;
    end

    % Getting the correct order
    order_asc = sort(order,'ascend');
end

C_matrix = zeros(par.U,par.timeslots);

C_matrix(order_asc,1) = 1;

C_matrix(:,2) = (C_matrix(:,1) == 0);

end











