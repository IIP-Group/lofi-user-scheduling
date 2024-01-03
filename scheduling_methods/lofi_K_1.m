% =========================================================================
% -- LoFi (K=1)
% -------------------------------------------------------------------------
%
% Description :
% (i) Draw two random schedules
% (ii) Pick the schedule that has the largest minimal per-UE SINR, i.e.,
% max(min(SINR_u)), u=1,...,U
% The dimensions of C_matrix are # of UEs x # of timeslots x # of candidate
% configurations, which in this case is 2 for the two random schedules. We
% output C_chosen with size # of UEs x # of timeslots, where C_(u,t) = 0
% means that the UE u doesn't transmit a signal in timeslot t and
% C_(u,t) = 1 means that the UE u transmits a signal in timeslot t.
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Reinhard Wiesmayr, Victoria Palhares, Christoph Studer
% -- e-mails: <wiesmayr@iis.ee.ethz.ch, palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================
function [C_chosen] = lofi_K_1(par,var)

C_matrix = zeros(par.U,par.timeslots,2);

% 1st candidate
chosen_ues_1 = randsample(par.U, par.U/par.timeslots).';
complement_1 = 1:par.U;
complement_1(chosen_ues_1) = [];

C_matrix(chosen_ues_1,1,1) = 1;
C_matrix(complement_1,2,1) = 1;


t=1;
Hc = var.H*diag(C_matrix(:,t,1));
L = inv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));

L_H = zeros(par.B,par.U);
H_L = zeros(par.U,par.B);
for n = 1:par.U
    L_H(:,n) = L*var.H(:,n);
    H_L(n,:) = var.H(:,n)'*L;
end


SINR_u = zeros(par.U,2);
for u=chosen_ues_1
    w_c_u = C_matrix(u,t,1)*(H_L(u,:)); % equalizer vector of the u-th user with scheduling
    h_c_u = C_matrix(u,t,1)*(var.H(:,u)); % channel vector of the u-th user with scheduling

    d_c = (var.Es*abs(w_c_u*h_c_u)^2); % desired signal
    i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- (abs(w_c_u*h_c_u)^2))); % interference
    n_c = (var.average_N0*norm(w_c_u)^2); % noise

    SINR_u(u,1) = d_c/(i_c+ n_c);
end

t=2;
Hc = var.H*diag(C_matrix(:,t,1));
L = inv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));

L_H = zeros(par.B,par.U);
H_L = zeros(par.U,par.B);
for n = 1:par.U
    L_H(:,n) = L*var.H(:,n);
    H_L(n,:) = var.H(:,n)'*L;
end

for u=complement_1
    w_c_u = C_matrix(u,t,1)*(H_L(u,:)); % equalizer vector of the u-th user with scheduling
    h_c_u = C_matrix(u,t,1)*(var.H(:,u)); % channel vector of the u-th user with scheduling

    d_c = (var.Es*abs(w_c_u*h_c_u)^2); % desired signal
    i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- (abs(w_c_u*h_c_u)^2))); % interference
    n_c = (var.average_N0*norm(w_c_u)^2); % noise

    SINR_u(u,1) = d_c/(i_c+ n_c);
end

min_candidate_1 = min(SINR_u(:,1));

% 2nd candidate
chosen_ues_2 = randsample(par.U, par.U/par.timeslots).';
complement_2 = 1:par.U;
complement_2(chosen_ues_2) = [];

C_matrix(chosen_ues_2,1,2) = 1;
C_matrix(complement_2,2,2) = 1;

t=1;
Hc = var.H*diag(C_matrix(:,t,2));
L = inv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));

L_H = zeros(par.B,par.U);
H_L = zeros(par.U,par.B);
for n = 1:par.U
    L_H(:,n) = L*var.H(:,n);
    H_L(n,:) = var.H(:,n)'*L;
end

% SINR_u = zeros(par.U,par.timeslots,2);
for u=chosen_ues_2
    w_c_u = C_matrix(u,t,2)*(H_L(u,:)); % equalizer vector of the u-th user with scheduling
    h_c_u = C_matrix(u,t,2)*(var.H(:,u)); % channel vector of the u-th user with scheduling

    d_c = (var.Es*abs(w_c_u*h_c_u)^2); % desired signal
    i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- (abs(w_c_u*h_c_u)^2))); % interference
    n_c = (var.average_N0*norm(w_c_u)^2); % noise

    SINR_u(u,2) = d_c/(i_c+ n_c);
end

t=2;
Hc = var.H*diag(C_matrix(:,t,2));
L = inv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));

L_H = zeros(par.B,par.U);
H_L = zeros(par.U,par.B);
for n = 1:par.U
    L_H(:,n) = L*var.H(:,n);
    H_L(n,:) = var.H(:,n)'*L;
end

for u=complement_2
    w_c_u = C_matrix(u,t,2)*(H_L(u,:)); % equalizer vector of the u-th user with scheduling
    h_c_u = C_matrix(u,t,2)*(var.H(:,u)); % channel vector of the u-th user with scheduling

    d_c = (var.Es*abs(w_c_u*h_c_u)^2); % desired signal
    i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- (abs(w_c_u*h_c_u)^2))); % interference
    n_c = (var.average_N0*norm(w_c_u)^2); % noise

    SINR_u(u,2) = d_c/(i_c+ n_c);
end

min_candidate_2 = min(SINR_u(:,2));

[~,idx_chosen] = max([min_candidate_1,min_candidate_2]);

C_chosen = C_matrix(:,:,idx_chosen);



end










