% =========================================================================
% -- LoFi (K=4)
% -------------------------------------------------------------------------
%
% Description :
% For all 2*K random restarts:
% (i) Draw one random schedule
% end loop
% (ii) Pick the schedule with the largest minimal per-UE SINR, i.e.,
% max(min(SINR_u)), u=1,...,U
% The dimensions of C_matrix are # of UEs x # of timeslots x # of candidate
% configurations, which is by definition of the paper 2K random schedules. 
% We output C_chosen with size # of UEs x # of timeslots, where
% C_(u,t) = 0 means that the UE u doesn't transmit a signal in timeslot t 
% and C_(u,t) = 1 means that the UE u transmits a signal in timeslot t.
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Reinhard Wiesmayr, Christoph Studer
% -- e-mails: <wiesmayr@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================
function [C_chosen] = lofi_K_N(par,var,num_restarts)

C_matrix = zeros(par.U,par.timeslots,num_restarts);
min_candidate = zeros(num_restarts,1);
SINR_u = zeros(par.U,num_restarts);

for i=1:num_restarts % random initializations
    % first time slot
    chosen_ues_i(1,:) = randsample(par.U, par.U/par.timeslots).'; 
    complement_i = 1:par.U;
    complement_i(chosen_ues_i(1,:)) = [];
    % second time slot
    chosen_ues_i(2,:) = complement_i; 

    C_matrix(chosen_ues_i(1,:),1,i) = 1;
    C_matrix(chosen_ues_i(2,:),2,i) = 1;

    for t=1:par.timeslots
        Hc = var.H*diag(C_matrix(:,t,i));
        L = inv(Hc*Hc'+(var.average_N0/var.Es)*eye(par.B));

        L_H = zeros(par.B,par.U);
        H_L = zeros(par.U,par.B);
        for n = 1:par.U
            L_H(:,n) = L*var.H(:,n);
            H_L(n,:) = var.H(:,n)'*L;
        end

        for u=chosen_ues_i(t,:)
            % equalizer vector of the u-th user with scheduling
            w_c_u = C_matrix(u,t,i)*(H_L(u,:)); 
            % channel vector of the u-th user with scheduling
            h_c_u = C_matrix(u,t,i)*(var.H(:,u)); 

            d_c = (var.Es*abs(w_c_u*h_c_u)^2); % desired signal
            % interference
            i_c = (var.Es*((abs(w_c_u*Hc).^2* ...
                ones(par.U,1))- (abs(w_c_u*h_c_u)^2))); 
            n_c = (var.average_N0*norm(w_c_u)^2); % noise

            SINR_u(u,i) = d_c/(i_c+ n_c);
        end
    end

    min_candidate(i) = min(SINR_u(:,i));
end


[~,idx_chosen] = max(min_candidate);

C_chosen = C_matrix(:,:,idx_chosen);



end










