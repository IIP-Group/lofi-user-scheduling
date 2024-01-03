% =========================================================================
% -- LoFi++ (K=4)
% -------------------------------------------------------------------------
%
% Description :
% For all K random restarts:
% (i) Draw one random schedule
% (ii) Determine the UE with the worst SINR in each timeslot
% (iii) Swap these two worst UEs to create a second schedule
% (iv) Among the two schedules, pick the schedule with the largest minimal 
% per-UE SINR, i.e.,max(min(SINR_u)), u=1,...,U
% end loop
% (v) Among all the K random restarts, pick the schedule with the largest
% minimal per-UE SINR, i.e., max(min(SINR_u)), u=1,...,U
% The dimensions of C_matrix are # of UEs x # of timeslots x # of candidate
% configurations, which in this case is 2K for the 2K schedules. We
% output C_chosen with size # of UEs x # of timeslots, where C_(u,t) = 0
% means that the UE u doesn't transmit a signal in timeslot t and
% C_(u,t) = 1 means that the UE u transmits a signal in timeslot t.
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Reinhard Wiesmayr, Christoph Studer
% -- e-mails: <wiesmayr@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================
function [C_chosen] = lofi_2_K_N(par,var,num_restarts)

snr_min_ar = zeros(num_restarts*2,1);

C_matrix = zeros(par.U,par.timeslots,2*num_restarts);

for i=1:num_restarts
    % only one candidate
    chosen_ues_1 = randsample(par.U, par.U/par.timeslots).';
    complement_1 = 1:par.U;
    complement_1(chosen_ues_1) = [];

    C_matrix(chosen_ues_1,1,2*i-1) = 1;
    C_matrix(complement_1,2,2*i-1) = 1;


    t=1;
    Hc = var.H*diag(C_matrix(:,t,2*i-1));
    L = inv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));

    L_H = zeros(par.B,par.U);
    H_L = zeros(par.U,par.B);
    for n = 1:par.U
        L_H(:,n) = L*var.H(:,n);
        H_L(n,:) = var.H(:,n)'*L;
    end


    SINR_u = zeros(par.U,2);
    for u=chosen_ues_1
        % equalizer vector of the u-th user with scheduling
        w_c_u = C_matrix(u,t,2*i-1)*(H_L(u,:)); 
        % channel vector of the u-th user with scheduling
        h_c_u = C_matrix(u,t,2*i-1)*(var.H(:,u)); 

        d_c = (var.Es*abs(w_c_u*h_c_u)^2); % desired signal
        % interference
        i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- ...
            (abs(w_c_u*h_c_u)^2))); 
        n_c = (var.average_N0*norm(w_c_u)^2); % noise

        SINR_u(u,1) = d_c/(i_c+ n_c);
    end

    t=2;
    Hc = var.H*diag(C_matrix(:,t,2*i-1));
    L = inv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));

    L_H = zeros(par.B,par.U);
    H_L = zeros(par.U,par.B);
    for n = 1:par.U
        L_H(:,n) = L*var.H(:,n);
        H_L(n,:) = var.H(:,n)'*L;
    end

    for u=complement_1
        % equalizer vector of the u-th user with scheduling
        w_c_u = C_matrix(u,t,2*i-1)*(H_L(u,:)); 
        % channel vector of the u-th user with scheduling
        h_c_u = C_matrix(u,t,2*i-1)*(var.H(:,u)); 

        d_c = (var.Es*abs(w_c_u*h_c_u)^2); % desired signal
        i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- ...
            (abs(w_c_u*h_c_u)^2))); % interference
        n_c = (var.average_N0*norm(w_c_u)^2); % noise

        SINR_u(u,1) = d_c/(i_c+ n_c);
    end


    % find swap index
    [min_candidate_1_1, idx_1] = min(SINR_u(chosen_ues_1,1)); % timeslot 1
    [min_candidate_1_2, idx_2] = min(SINR_u(complement_1,1)); % timeslot 2
    idx_1 = chosen_ues_1(idx_1); % user index
    idx_2 = complement_1(idx_2);

    % swap
    chosen_ues_2 = chosen_ues_1;
    chosen_ues_2(chosen_ues_1==idx_1) = idx_2;
    complement_2 = complement_1;
    complement_2(complement_1==idx_2) = idx_1;

    % evaluate after swap
    C_matrix(chosen_ues_2,1,2*i) = 1;
    C_matrix(complement_2,2,2*i) = 1;

    t=1;
    Hc = var.H*diag(C_matrix(:,t,2*i));
    L = inv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));

    L_H = zeros(par.B,par.U);
    H_L = zeros(par.U,par.B);
    for n = 1:par.U
        L_H(:,n) = L*var.H(:,n);
        H_L(n,:) = var.H(:,n)'*L;
    end

    % SINR_u = zeros(par.U,par.timeslots,2);
    for u=chosen_ues_2
        % equalizer vector of the u-th user with scheduling
        w_c_u = C_matrix(u,t,2*i)*(H_L(u,:)); 
        % channel vector of the u-th user with scheduling
        h_c_u = C_matrix(u,t,2*i)*(var.H(:,u));

        d_c = (var.Es*abs(w_c_u*h_c_u)^2); % desired signal
        i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- ...
            (abs(w_c_u*h_c_u)^2))); % interference
        n_c = (var.average_N0*norm(w_c_u)^2); % noise

        SINR_u(u,2) = d_c/(i_c+ n_c);
    end

    t=2;
    Hc = var.H*diag(C_matrix(:,t,2*i));
    L = inv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));

    L_H = zeros(par.B,par.U);
    H_L = zeros(par.U,par.B);
    for n = 1:par.U
        L_H(:,n) = L*var.H(:,n);
        H_L(n,:) = var.H(:,n)'*L;
    end

    for u=complement_2
        % equalizer vector of the u-th user with scheduling
        w_c_u = C_matrix(u,t,2*i)*(H_L(u,:)); 
        % channel vector of the u-th user with scheduling
        h_c_u = C_matrix(u,t,2*i)*(var.H(:,u)); 

        d_c = (var.Es*abs(w_c_u*h_c_u)^2); % desired signal
        i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- ...
            (abs(w_c_u*h_c_u)^2))); % interference
        n_c = (var.average_N0*norm(w_c_u)^2); % noise

        SINR_u(u,2) = d_c/(i_c+ n_c);
    end

    SINR_mins = min(SINR_u, [], 1); % get the minimum of both trials
    snr_min_ar(2*i-1:2*i) = SINR_mins(:);
end

[~, idx_max] = max(snr_min_ar);

C_chosen = C_matrix(:,:,idx_max);  % idx_chosen must be 1

end









