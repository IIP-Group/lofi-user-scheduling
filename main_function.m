% =========================================================================
% -- LoFi User Scheduling for Multiuser MIMO Wireless Systems
% -------------------------------------------------------------------------
%
% Description :
% User scheduling simulator for each channel realization, SNR point and
% scheduling method
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function main_function(par,proc_ID)

% Finding channel ID, SNR and scheduling method for this process
var.id = proc_ID+1;
parameters = par.parameters(var.id,:); % getting columns of process line
var.ch = str2num(parameters(1)); % channel id
var.SNR_dB = str2num(parameters(2)); % SNR in dB
var.scheduling_method = parameters(3); % chosen scheduling method

snr_idx = find(par.SNRdB_list == var.SNR_dB); % getting SNR index
% getting scheduler index
scheduler_idx = find(par.scheduling_options == var.scheduling_method);

index = (var.ch-1)*length(par.SNRdB_list)+snr_idx;
rng(1000*index + 0, 'twister');

% Setting modulation scheme
% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.modulation)
    case 'BPSK'
        par.symbols = [ -1 1 ];
    case 'QPSK'
        par.symbols = [ -1-1i,-1+1i, ...
            +1-1i,+1+1i ];
    case '16QAM'
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM'
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    otherwise
        error('The variable par.modulation is not defined.')
end
% extract average symbol energy
var.Es = mean(abs(par.symbols).^2);

% precompute bit labels
var.Q = log2(length(par.symbols)); % number of bits per symbol
var.bits = de2bi(0:length(par.symbols)-1,var.Q,'left-msb');

% generate channels with power control
[var] = gen_channel(par,var);
var.H_full_CSI = var.H; % perfect CSI

% Average noise variance
var.H_norm = norm(var.H,'fro')^2*(par.Us/par.U);
var.average_N0 = (var.H_norm/par.B)*var.Es*10^(-var.SNR_dB/10);

% Channel estimation or perfect channel knowledge
switch(par.channel_estimation)
    case 'BEACHES'
        N0_est = var.average_N0/(par.Us*var.Es);
        channel_noise = (randn(par.B,par.U)+1i*randn(par.B,par.U));
        Hnoisy = var.H + channel_noise*sqrt(N0_est/2);
        [Hest, ~] = BEACHES(par,Hnoisy,N0_est,'antenna_domain');
        var.H = Hest; % replace ground truth channel by estimated channel
    case 'Perfect'
        var.H = var.H_full_CSI;
    otherwise
        error('Please choose a channel estimation method');
end

% Reinitalizing seed
rng(1000*index + 0, 'twister');
tic;
% Select scheduling algorithm:
% the schedulers below define whether a UE transmits in timeslot t or
% not. If C_(u,t) = 0, the UE u doesn't transmit a signal in timeslot t.
% If C_(u,t) = 1, the UE u transmits a signal in timeslot t. The
% scheduling matrix C has size par.U X par.timeslots
switch (var.scheduling_method)
    case "No Scheduling"
        % No Scheduling
        C = ones(par.U,par.timeslots);
    case "Random"
        % Random scheduling
        [C_rand,~] = random_scheduling(par,var);
        C = C_rand;
    case "SUS"
        % SUS based scheduling
        [C_sus] = SUS(par,var);
        C = C_sus;
    case "CSS"
        % CSS based scheduling
        [C_css] = CSS(par,var);
        C = C_css;
    case "Greedy"
        % Greedy scheduling
        [C_greedy] = greedy(par,var);
        C = C_greedy;
    case "Opt.-based"
        % Optimization-based scheduling
        var.regularization = "L2 Norm";
        [par,var,~,C_quant,~,~,~,~] = ...
            optimization_based_scheduling(par,var);
        C = C_quant;
    case "Exhaustive"
        % Exhaustive Search
        [C_es,~] = exhaustive_search(par,var);
        C = C_es(:,:,1);
    case "LoFi (K=1)"
        % LoFi (K=1)
        [C_lofi_K_1] = lofi_K_1(par,var);
        C = C_lofi_K_1;
    case "LoFi (K=4)"
        % LoFi (K=4)
        % definition in the paper assumes K=2*num_restarts
        num_restarts=8;
        [C_lofi_K_N] = lofi_K_N(par,var,num_restarts);
        C = C_lofi_K_N;
    case "LoFi++ (K=1)"
        % LoFi++ (K=1)
        [C_lofi_2_K_1] = lofi_2_K_1(par,var);
        C = C_lofi_2_K_1;
    case "LoFi++ (K=4)"
        % LoFi++ (K=4)
        num_restarts = 4;
        [C_lofi_2_K_N] = lofi_2_K_N(par,var,num_restarts);
        C = C_lofi_2_K_N;
    otherwise
        error('The variable var.scheduling_method is not defined');
end
elapsedTime = toc; % runtime of scheduling methods

% Reinitalizing seed to generate exactly the same symbols and noise,
% considering a certain C.
rng(1000*index + 0, 'twister');

% Generating random transmitted bits
var.fixed_bits = randi([0 1],par.U,var.Q,par.timeslots,par.trials); % one symbol vector for every trial and timeslot

% initialize results
ber = zeros(par.U,par.timeslots,par.trials);

for t=1:par.timeslots % Transmission timeslots. In this work, timeslots = 2
    Dc = diag(C(:,t));
    Hc = var.H*Dc;

    % Real noise variance
    var.real_H_norm(t) = norm(Hc,'fro')^2;
    var.real_N0(t) = (var.real_H_norm(t)/par.B)*var.Es*10^(-var.SNR_dB/10);

    %MMSE equalizer with scheduling
    Wc = pinv(Hc*Hc'+(var.real_N0(t)/var.Es)*eye(par.B))*Hc;

    % Transmission trials. We only vary symbol and noise vector
    for tt=1:par.trials
        % generate transmit symbol
        idx = bi2de(var.fixed_bits(:,:,t,tt),'left-msb')+1;
        s = par.symbols(idx).';

        % Generating Noise
        n = sqrt(0.5*var.real_N0(t))*(randn(par.B,1)+1i*randn(par.B,1));

        yc = Hc*s+n;
        s_schedule = Dc*s;
        s_hat_schedule = Dc*Wc'*yc;
        var.vector = real(diag(Dc*Wc'*Hc));

        % Calculate BER
        [ber(:,t,tt)] = calculate_ber(par,var,C(:,t),s_hat_schedule,t,tt);
    end
end

% Average N0
results.avg_real_N0 = mean(var.real_N0);

% BER
timeslots_per_user = sum(C,2);
sum_all_C = sum(sum(C));
% sum over timeslots and then sum over all trials
results.ber_per_user = sum(sum(ber,2),3)/(var.Q*par.trials.*timeslots_per_user);
% average rate per user
results.ber_all_users = sum(ber(:))/(var.Q*par.trials*sum_all_C);


%% Save file
index_process = ((scheduler_idx-1)*(length(par.SNRdB_list)*par.channels)+index)-1;

FileName=[par.results_path,'USER_SCHEDULING_',num2str(index_process),'.mat'];
save (FileName,'par','var','results','elapsedTime');


end



