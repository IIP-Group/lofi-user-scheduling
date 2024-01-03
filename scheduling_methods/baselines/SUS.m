% =========================================================================
% -- Semiorthogonal User Scheduling (SUS)
% -------------------------------------------------------------------------
% Reference:
%          Taesang Yoo and A. Goldsmith, "On the optimality of
%          multiantenna broadcast scheduling using zero-forcing
%          beamforming," in IEEE Journal on Selected Areas in
%          Communications, vol. 24, no. 3, pp. 528-541, March 2006,
%          doi: 10.1109/JSAC.2005.862421.
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function  [C_matrix] = SUS(par,var)

epsilon = 0.5;
Ns = par.Us;

order_asc = 0;
while (length(order_asc) < Ns)
    epsilon = epsilon+0.1;

    % Step 1
    U_group = 1:par.U;
    H_sus = [];
    g_0 = [];
    order = [];

    i = 1;

    Tlen = [];
    while(i<=Ns && isempty(U_group)~=1)
        Tlen(i) = length(U_group);
        g = zeros(par.B,length(U_group));
        g_norm = zeros(length(U_group),1);
        for k = 1:length(U_group)

            % Step 2
            temp = zeros(par.B);
            for j=1:i-1
                temp = temp + (g_0(:,j)*g_0(:,j)')/norm(g_0(:,j))^2;
            end

            g(:,k) = (eye(par.B)-temp) * var.H(:,U_group(k));
            g_norm(k) = norm(g(:,k));

        end
        % Step 3
        [~,sel] = max(g_norm);
        order(i) = U_group(sel);

        H_sus = [H_sus var.H(:,U_group(sel))];
        g_0 = [g_0 g(:,sel)];

        % Step 4
        U_group(sel) = [];
        temp = [];
        for j = 1:length(U_group)
            if abs(var.H(:,U_group(j))'*g_0(:,i))/(norm(var.H(:,U_group(j)))*norm(g_0(:,i))) < epsilon
                temp = [temp U_group(j)];
            end

        end
        U_group = temp;

        i=i+1;
    end

    % Getting the correct order
    order_asc = sort(order,'ascend');

end

C_matrix = zeros(par.U,par.timeslots);
C_matrix(order_asc,1) = 1;

C_matrix(:,2) = (C_matrix(:,1) == 0);

end