# Simulator for LoFi User Scheduling for Multiuser MIMO Wireless Systems 
(c) 2023 Victoria Palhares, Reinhard Wiesmayr
e-mail: palhares@iis.ee.ethz.ch, wiesmayr@iis.ee.ethz.ch


### Important information 

If you are using this simulator (or parts of it) for a publication, then you *must* cite our paper:

A. Gallyas-Sanhueza, G.Marti, V. Palhares, R. Wiesmayr, and C. Studer, "LoFi User Scheduling for Multiuser MIMO Wireless Systems," 2024 IEEE International Conference on Acoustics, Speech and Signal Processing, Seoul, South Korea, 2024.

and clearly mention this in your paper. In the following, by *paper*, we mean the above paper.

### How to start a simulation:

To regenerate the plots in figures 3 and 4 of the paper, simply assign the letter `'a'` to the variable `par.sim_scenario` in the `main_sim_code` 
and run the script. It will generate two plots with BER curves.
Also note that the list of schedulers to be simulated is set in the `par.scheduling_chosen` variable in the `param_config.m` file. 

In order to simulate with parameters other than our predefined values, define a new set of parameters in the `switch case` block of `param_config.m` file and assign a value to the variable `par.sim_scenario` in the `main_sim_code.m` that corresponds to the new set of parameters. Please consider the following guidelines when creating a new set of simulation parameters:
- Set the number of BS antennas and UEs and other parameters in the `param_config.m` file
- Our channel dataset, which has been generated with REMCOM Wireless InSite, can be found in `channel/Channels_B_16_U_16_Ch_100.mat`. The channel dataset is of size `16*16*100`, where the first dimension specifies 16 basestation receive antennas, the second 16 single-antenna user equipments (UEs), and the third dimension represents 100 channel realizations. In the simulated outdoor scenario, we consider 22448 possible UE positions. For every channel realization, we picked 16 UE positions uniformly and independently at random. Note that the number of channel realizations in `par.channels` has to be set to an integer less than or equal to 100, unless other channels are provided in the `channel` folder
- The variable `par.scheduling_options` denote all the available methods in the simulator while `par.scheduling_chosen` denote all the methods that will be simulated
- We perform `par.trials` transmission trials, where  we fix the channel and scheduler and vary the symbols and noise
- We perform user scheduling once per channel realization
- Our simulator is fixed to `par.timeslots = 2`, but our framework could be extended in future work to more timeslots
- The dynamic range of the receive power from the weakest to the strongest UE can be set in the variable `par.P_db` as a decibel value
- The modulation scheme can be changed to BPSK, QPSK, or 64QAM in the variable `par.modulation`
- The channel estimation can be disabled by setting the variable `par.channel_estimation` to 'Perfect'
   
	
- The overall number of processes is given by the number of `par.channels*length(par.SNRdB_list)*length(par.scheduling_chosen)`
- You need to modify your `par.results_path` variable in `main_sim_code` according to where you want to save your result files

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.


### Version history
* Version 1 (Jan. 03 2024) - palhares@iis.ee.ethz.ch  - initial publication of the simulator
