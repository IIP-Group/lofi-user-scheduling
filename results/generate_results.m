% =========================================================================
% -- Generate Figure 3,4 and Table 1
% -------------------------------------------------------------------------
%
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function generate_results(par)

addpath(par.results_path);

current_path = pwd;


index_process = [];
for type=1:length(par.scheduling_chosen)
    scheduler_idx = find(strcmp(par.scheduling_options,par.scheduling_chosen{type}));
    for c = 1:par.channels
        for int_snr=1:length(par.SNRdB_list)
            index = (c-1)*length(par.SNRdB_list)+int_snr;
            index_process = [index_process;((scheduler_idx-1)*(length(par.SNRdB_list)*par.channels)+index)-1];
        end
    end
end

% Getting all the files for all the processes
file = cell(length(index_process));
for i=1:length(index_process)
    file_name = [par.results_path,['USER_SCHEDULING_',num2str(index_process(i)),'.mat']];
    file{i} = file_name;
end

res.BER_all_users = cell(length(par.scheduling_chosen),1);
res.BER_per_SNR = cell(length(par.scheduling_chosen),1);
res.runtime = cell(length(par.scheduling_chosen),1);
res.runtime_per_SNR = cell(length(par.scheduling_chosen),1);

% Collecting the data for all setups and all processes
for n_proc=1:length(index_process)
    load(file{n_proc},'par','var','results','elapsedTime');
    par.scheduling_options = ["No Scheduling","Random","SUS","CSS","Greedy","Opt.-based","Exhaustive",...
        "LoFi (K=1)","LoFi (K=4)","LoFi++ (K=1)","LoFi++ (K=4)"];
    par.scheduling_chosen = ["No Scheduling","Random","SUS","CSS","Greedy","Opt.-based","Exhaustive",...
        "LoFi (K=1)","LoFi (K=4)","LoFi++ (K=1)","LoFi++ (K=4)"];

    ind_style = find(par.scheduling_chosen == var.scheduling_method);

    res.BER_all_users{ind_style} = [res.BER_all_users{ind_style};results.ber_all_users]; % for BER vs. SNR
    res.runtime{ind_style} = [res.runtime{ind_style};elapsedTime]; % for runtime table
end


par.parameters = [];
for type=1:length(par.scheduling_chosen)
    for c = 1:par.channels
        for int_snr=1:length(par.SNRdB_list)
            par.parameters = [par.parameters;c,par.SNRdB_list(int_snr),par.scheduling_chosen(type)];
        end
    end
end

% Average over channel realizations
parameters = cell(length(par.scheduling_chosen));
for i=1:length(par.scheduling_chosen)
    indexes = find(par.parameters(:,3) == par.scheduling_chosen(i));
    parameters{i} = par.parameters(indexes,:);

    for j=1:length(par.SNRdB_list)
        indexes_2 = find(str2double(parameters{i}(:,2)) == par.SNRdB_list(j));
        res.BER_per_SNR{i}(:,j) = mean(res.BER_all_users{i}(indexes_2));
        res.runtime_per_SNR{i}(:,j) = mean(res.runtime{i}(indexes_2));
    end
end

% Prepare Table 1 for Latex
res.table_string = '';
for i=1:length(par.scheduling_chosen)
    res.table_string = ([res.table_string,par.scheduling_chosen{i},' & $' num2str(round(res.runtime_per_SNR{i}(1,4),3)),'\\']);
end

%% Figure 3 and 4

colors = define_colors();

% Figure 3
style_color_1  = {colors.dark_grey,colors.light_grey,colors.green,colors.green,colors.green,colors.orange,colors.gold,colors.navy};
style_line_1   = {'-','-','--','--','-','-','-','-'};
style_marker_1 = {'none','s','+','o','d','h','s','s'};
legends_1 = {'No scheduling','Random','SUS [4]','CSS [3]','Greedy [3]','Opt.-based [17]','Exhaustive','LoFi++ (K=4)'};
scheduling_chosen_1 = ["No Scheduling","Random","SUS","CSS","Greedy","Opt.-based","Exhaustive","LoFi++ (K=4)"];

% Figure 4
style_color_2  = {colors.light_grey,colors.magenta,colors.magenta,colors.navy,colors.navy};
style_line_2   = {'-','--','-','--','-'};
style_marker_2 = {'s','d','d','s','s'};
legends_2 = {'Random','LoFi (K=1)','LoFi (K=4)','LoFi++ (K=1)','LoFi++ (K=4)'};
scheduling_chosen_2 = ["Random","LoFi (K=1)","LoFi (K=4)","LoFi++ (K=1)","LoFi++ (K=4)"];


% Figure 3
h = figure(1);
for i=1:length(scheduling_chosen_1)
    ind_curve = find(par.scheduling_chosen == scheduling_chosen_1(i));
    ind_style = find(scheduling_chosen_1 == scheduling_chosen_1(i));
    temp = res.BER_per_SNR{ind_curve};
    semilogy(par.SNRdB_list,temp,'Color',style_color_1{ind_style},'LineStyle',style_line_1{ind_style},...
        'Marker',style_marker_1{ind_style},'LineWidth',2,'DisplayName', legends_1{ind_style});
    hold on;
end
legend('Location','southwest')
legend('show', 'FontSize',12);
xlabel('receive SNR [dB]');
ylabel('uncoded bit error-rate (BER)');
axis([5 25 1e-3 1]);
grid;
savefig(h,[current_path,'/results/figures/Figure_3'])
clear('h')

% Figure 4
h = figure(2);
for i=1:length(scheduling_chosen_2)
    ind_curve = find(par.scheduling_chosen == scheduling_chosen_2(i));
    ind_style = find(scheduling_chosen_2 == scheduling_chosen_2(i));
    temp = res.BER_per_SNR{ind_curve};
    semilogy(par.SNRdB_list,temp,'Color',style_color_2{ind_style},'LineStyle',style_line_2{ind_style},...
        'Marker',style_marker_2{ind_style},'LineWidth',2,'DisplayName', legends_2{ind_style});
    hold on;
end
legend('Location','southwest')
legend('show', 'FontSize',12);
xlabel('receive SNR [dB]');
ylabel('uncoded bit error-rate (BER)');
axis([5 25 1e-3 1]);
grid;
savefig(h,[current_path,'/results/figures/Figure_4'])
clear('h');


FileName=[datestr(now, 'yyyy_mmm_dd'),'_USER_SCHEDULING.mat'];
save (FileName,'res');

end