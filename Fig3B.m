% time intervals (in seconds)
timesteps = (0:0.2:10);
network_time_intervals = mat2cell([timesteps(1:end-1)',timesteps(2:end)'],ones([1,50]),2);

mice=length(mouse_IDs);

% define thresholds
network_p_threshold = 0.01;
network_performance_threshold = 0.7;
pagerank_constant = 0.85;
network_frequency_threshold = 0.5;
% define as 0 for naive animals i.e. pool sessions below the performance threshold
% define as 1 for expert animals i.e. pool sessions above the performance threshold
perf_def = 0;

% cannel labels for pooling
% get selected channels
label=[];
channel_label = {'GP', 'DG', 'VL','Cpu', 'CA1_Py', 'S1_bf','VM', 'LDVL', 'Rt','M1','CA1_Mol','dCA1'};

comb=length(channel_label);
        
pool_mice=[1 2 3 4 5 6 7 8 9 10 11 12 13];
clear outdegree_vector;

Accuracy=[];
link_to_error=[];
f1_score=[];

for count_time_interval = 1:length(network_time_intervals)
    
count_i_mice=7

count_mice=pool_mice(count_i_mice);
session_names=sessions_IDs{1,count_mice};
sessions=length(session_names);

% --- import channel labels ---
channel_labels=channel_labels_all{1,count_mice};

% --- Choose time intervals (in seconds) ---
% CUE=[1,1.5], TEX=[2,3.5], LIK=[3.5,6]
% pre_texture [1.5,2.0], textue [3.0,3.5]
%network_time_interval = {[4.0,4.5]};

% --- choose NIFTy, bivariate and multivariate measures ---
% --- Generate adjacency matrices for GO and NOGO trials ---
    
clear network_adj_matrices

for trials_type = {'Hit'}

ii=1;
count=[];
% check that every folder contains results file

count_sessions = 16
            
                if numel(fieldnames(data_TE{1,count_mice}{1,count_sessions}))==3
                
                if and(numel(fieldnames(data_TE{1,count_mice}{1,count_sessions}.NIfTy))==2,...
                    numel(fieldnames(data_TE{1,count_mice}{1,count_sessions}.bivar))==2)
                
                count(ii)=count_sessions;
                clear temp_data
                
                temp_data.parameteres=data_TE{1,count_mice}{1,count(ii)}.NIfTy.Hit.results.parameters;
                % select timesteps in the desired time interval
                %temp_data.parameteres.samples_original/20;
%                 network_time_steps = find(temp_data.parameteres.samples_original/20 >= min(network_time_intervals{count_time_interval,1})...
%                     & temp_data.parameteres.samples_original/20 <= max(network_time_intervals{count_time_interval,1}));

                network_time_steps =count_time_interval;
                
                if length(network_time_steps) < 1
                    error('ERROR: no data points found in the specified time interval.')
                end
                
                % store the new adjacency matrix (assumption: all training sessions have the same number of channels)
                if strcmp(trials_type{1,1},'Hit')
                    
                    temp_data.bivar=data_TE{1,count_mice}{1,count(ii)}.bivar.Hit.results;
                    network_adj_matrices.GO.bivar(:,:,count(ii)) = network_bivar_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.multivar=data_TE{1,count_mice}{1,count(ii)}.multivar.Hit.results;
                    test_adj_matrices.GO.multivar(:,:,count(ii)) = network_bivar_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);

                elseif strcmp(trials_type{1,1},'CR')
                    
                    temp_data.bivar=data_TE{1,count_mice}{1,count(ii)}.bivar.CR.results;
                    network_adj_matrices.NOGO.bivar(:,:,count(ii)) = network_bivar_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.multivar=data_TE{1,count_mice}{1,count(ii)}.multivar.CR.results;
                    test_adj_matrices.NOGO.multivar(:,:,count(ii)) = network_bivar_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);
                                        
                    ii=ii+1;
                end
                end
                
                end
end
end
   
%% barplot
channel_labels_all{1, 7} 
subplot(6,1,1)
% Sbf to CPu
bar(squeeze(temp_data.multivar.TE_table(:,4,1)))
title('Sbf to CPu')
ylim([0 0.1])

subplot(6,1,2)
% CPu to GP
bar(squeeze(temp_data.multivar.TE_table(:,1,2)))
title('CPu to GP')
ylim([0 0.1])

subplot(6,1,3)
% GP to VL
bar(squeeze(temp_data.multivar.TE_table(:,2,6)))
title('GP to VL')
ylim([0 0.1])

subplot(6,1,4)
% VL to S1bf
bar(squeeze(temp_data.multivar.TE_table(:,6,4)))
title('VL to S1bf')
ylim([0 0.1])

subplot(6,1,5)
% VM to S1bf
bar(squeeze(temp_data.multivar.TE_table(:,5,4)))
title('VM to S1bf')
ylim([0 0.1])

subplot(6,1,6)
% VM to S1bf
bar(squeeze(temp_data.multivar.TE_table(:,4,8)))
title('S1bf to M1')
ylim([0 0.1])

%% 

imagesc(squeeze(temp_data.multivar.TE_table(:,5,:)))