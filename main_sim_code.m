% =========================================================================
% -- LoFi User Scheduling for Multiuser MIMO Wireless Systems
% -------------------------------------------------------------------------
% editor: Victoria Palhares
%
% Description :
% Top level file of this user scheduling simulator that sets the parameters
% and launches the simulation (locally on this machine)
% -------------------------------------------------------------------------
% Last Updated: 23/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

% Set up scenario described in param_config:
% Scenario 'a' is the one described in the paper with par.B=16, par.U=16, 
% and par.timeslots=2.
% To define scenarios other than the default one, add a new scenario 
% 'b','c',... and customize the parameters in the switch(par.sim_scenario) 
% block in the param_config file. 
% From the provided channel dataset, we can set par.B<=16, par.U<=16, and 
% par.channels <=100. Our simulator is fixed to par.timeslots = 2, but our 
% framework could be extended in future work to more timeslots
par.sim_scenario = 'a';

% Define parameters of par.sim_scenario
par = param_config(par);

% Add necessary folders
addpath('channel/');
addpath('results/');
addpath('results/figures/');
addpath('scheduling_methods/');
addpath('scheduling_methods/baselines/');
addpath('scheduling_methods/baselines/optimization_based_scheduling/');

%% Calculating total number of processes
par.parameters = []; % matrix describes channel realization, SNR point in dB and scheduling method
% Test different combinations:
for type=1:length(par.scheduling_chosen)
    for c = 1:par.channels
        for int_snr=1:length(par.SNRdB_list)
            par.parameters = [par.parameters;c,par.SNRdB_list(int_snr),par.scheduling_chosen(type)];
        end
    end
end

%%
% Modify this path to the folder where you save your processes
par.results_path = '/usr/scratch/pedrabonita/vmenescal/ICASSP/jobs_github/';

%% Launch the simulator
for kk=0:(par.total_jobs-1)
    fprintf('Process ID : %d\n', kk);
    % Run user scheduling simulator
    main_function(par,kk);
end
% Aggregate results from processes, plot Fig. 3 and 4, and create runtime
% table
generate_results(par);

fprintf('Simulation ended at %s. Check results in results/figures/ \n', datestr(now,'HH:MM:SS.FFF'));


