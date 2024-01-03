% =========================================================================
% -- Defining the simulation parameters
% -------------------------------------------------------------------------
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function par = param_config(par)

switch(par.sim_scenario)
    case 'a' % Simulation setup presented in the paper
        par.channels = 100; % # of channel realizations. We compute the user scheduler once per channel realization
        par.SNRdB_list = 5:5:25; % SNR list in dB (5 points)
        par.scheduling_options = ["No Scheduling","Random","SUS","CSS","Greedy","Opt.-based","Exhaustive",...
            "LoFi (K=1)","LoFi (K=4)","LoFi++ (K=1)","LoFi++ (K=4)"];
        par.scheduling_chosen = ["No Scheduling","Random","SUS","CSS","Greedy","Opt.-based","Exhaustive",...
            "LoFi (K=1)","LoFi (K=4)","LoFi++ (K=1)","LoFi++ (K=4)"];
        par.total_jobs = par.channels*length(par.SNRdB_list)*length(par.scheduling_options);
        par.trials = 100000; % # transmission trials: # of times that we vary the symbol and noise vector per channel realization
        % System parameters
        par.B = 16; % # of basestation antennas
        par.U = 16; % # of users equipments
        par.timeslots = 2; % # of timeslots
        % Power control
        par.pwr_cntrl = 'on';
        par.P_db = 6; % dynamic range of the receive power from the weakest to the strongest UE in dB
        % Modulation scheme
        par.modulation = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
        par.channel_estimation = 'BEACHES'; % "Perfect" for Perfect-CSI or "BEACHES" for channel estimation with the BEACHES algorithm
    otherwise
        error('Variable par.sim_scenario is not defined.');
end

% Users per timeslot
par.min_sum_columns = (par.U/par.timeslots); % Minimum # of users per timeslot
par.max_sum_columns = par.min_sum_columns; % Maximum # of users per timeslot
par.min_sum_rows = 1; % Minimum # of transmissions per user
par.max_sum_rows = par.min_sum_rows; % Maximum # of transmissions per user

par.Us = par.min_sum_columns; % For the case of 16 UEs, we have 8 UEs per timeslot

par.visualize_opt_graphs = 'off'; % 'on' or 'off'. If you want to see gradient descent graphs from opt.-based scheduler

end