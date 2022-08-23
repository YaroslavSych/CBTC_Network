% time intervals (in seconds)
timesteps = (0:0.2:10);
network_time_intervals = mat2cell([timesteps(1:end-1)',timesteps(2:end)'],ones([1,50]),2);

mice=length(mouse_IDs);
lcc_GO_size = NaN([100,mice,length(network_time_intervals)]);
lcc_NOGO_size = NaN([100,mice,length(network_time_intervals)]);

% cannel labels for pooling
% get selected channels
label=[];
channel_label = {'GP', 'DG', 'VL','Cpu', 'CA1_Py', 'S1_bf','VM', 'LDVL', 'Rt','M1','CA1_Mol','dCA1'};
comb=length(channel_label);

% define thresholds
network_p_threshold = 0.01;
network_performance_threshold = 0.7;
pagerank_constant = 0.85;
network_frequency_threshold = 0.5;
      
total_links=cell(1,length(network_time_intervals));
shared_links=cell(1,length(network_time_intervals));
network_links_shared_p_value=cell(1,length(network_time_intervals));

pool_mice=[1 2 3 4 5 6 7 8 9 10 11 12 13];

for count_time_interval = 1:length(network_time_intervals)

for count_i_mice=1:length(pool_mice)

count_mice=pool_mice(count_i_mice);
session_names=sessions_IDs{1,count_mice};
sessions=length(session_names);

% --- import channel labels ---
channel_labels=channel_labels_all{1,count_mice};

% --- Choose time intervals (in seconds) ---
%network_time_interval = {[1,1.5],[3,3.5],[6,6.5]};

% --- choose NIFTy, bivariate and multivariate measures ---
% --- Generate adjacency matrices for GO and NOGO trials ---
    
clear network_adj_matrices

for trials_type = {'Hit','CR'}

ii=1;
count=[];
% check that every folder contains results file

    for count_sessions = 1:sessions
            
                if numel(fieldnames(data_TE{1,count_mice}{1,count_sessions}))==3
                
                if and(numel(fieldnames(data_TE{1,count_mice}{1,count_sessions}.NIfTy))==2,...
                    numel(fieldnames(data_TE{1,count_mice}{1,count_sessions}.bivar))==2)
                
                count(ii)=count_sessions;
                clear temp_data
                
                temp_data.parameteres=data_TE{1,count_mice}{1,count(ii)}.NIfTy.Hit.results.parameters;
                % select timesteps in the desired time interval
                network_time_steps = find(temp_data.parameteres.samples_timesteps >= min(network_time_intervals{count_time_interval,1})...
                    & temp_data.parameteres.samples_timesteps < max(network_time_intervals{count_time_interval,1}));
                if length(network_time_steps) < 1
                    error('ERROR: no data points found in the specified time interval.')
                end
                
                % store the new adjacency matrix (assumption: all training sessions have the same number of channels)
                if strcmp(trials_type{1,1},'Hit')
                    
                    temp_data.bivar=data_TE{1,count_mice}{1,count_sessions}.bivar.Hit.results;
                    network_adj_matrices.GO.bivar(:,:,ii) = network_bivar_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.multivar=data_TE{1,count_mice}{1,count_sessions}.multivar.Hit.results;
                    network_adj_matrices.GO.multivar(:,:,ii) = network_bivar_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);

%                     temp_data.NIfTy=data_TE{1,count_mice}{1,count_sessions}.NIfTy.Hit.results;
%                     network_adj_matrices.GO.NIfTy(:,:,ii) = network_threshold(temp_data.NIfTy.p_table,network_time_steps,network_p_threshold);
                   
                    
                elseif strcmp(trials_type{1,1},'CR')
                    
                    temp_data.bivar=data_TE{1,count_mice}{1,count_sessions}.bivar.CR.results;
                    network_adj_matrices.NOGO.bivar(:,:,ii) = network_bivar_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.multivar=data_TE{1,count_mice}{1,count_sessions}.multivar.CR.results;
                    network_adj_matrices.NOGO.multivar(:,:,ii) = network_bivar_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);
                    
%                     temp_data.NIfTy=data_TE{1,count_mice}{1,count_sessions}.NIfTy.CR.results;
%                     network_adj_matrices.NOGO.NIfTy(:,:,ii) = network_threshold(temp_data.NIfTy.p_table,network_time_steps,network_p_threshold);
                    
                end
                ii=ii+1;

                end
                
                end
    end
end
   
% --- Group sessions by performance threshold ---
high_perf_sessions = cell2mat(performance{1,count_mice}) > network_performance_threshold;
high_perf_sessions_indices =find(high_perf_sessions(count) > 0);
    
low_perf_sessions = cell2mat(performance{1,count_mice}) <= network_performance_threshold;
low_perf_sessions_indices = find(low_perf_sessions(count) > 0);


% --- Give bivarite adj matrx to calculate network properties for network time intervals ---

if size(network_adj_matrices.GO.multivar,1)==12
    adj_input_GO=network_adj_matrices.GO.multivar;
    adj_input_NOGO=network_adj_matrices.NOGO.multivar;
else
    channel_sub_network=[26 31 30 25 34 46 44 47 15 48 33 35];
    adj_input_GO=network_adj_matrices.GO.multivar(channel_sub_network,channel_sub_network,:);
    adj_input_NOGO=network_adj_matrices.NOGO.multivar(channel_sub_network,channel_sub_network,:);
end

% --- calculate network properties for EXPERT sessoins ---
% --- Hit trials ---
network_links_GO=calculate_links_evolution(adj_input_GO(:,:,high_perf_sessions_indices));

total_links{1,count_time_interval}{1,count_i_mice}.expert.Hit= mean(network_links_GO{1,1},1,'omitnan');
shared_links{1,count_time_interval}{1,count_i_mice}.expert.Hit= mean(network_links_GO{1,2},1,'omitnan');
network_links_shared_p_value{1,count_time_interval}{1,count_i_mice}.expert.Hit= mean(network_links_GO{1,3},1,'omitnan');
% network_transitivity{1,count_time_interval}{1,count_i_mice}.expert.Hit= mean(network_links_GO{1,4},1,'omitnan');
% network_assortativity{1,count_time_interval}{1,count_i_mice}.expert.Hit=mean(network_links_GO{1,5},1,'omitnan');
% network_modularity{1,count_time_interval}{1,count_i_mice}.expert.Hit=mean(network_links_NOGO{1,6},1,'omitnan');

% --- CR trials ---
network_links_NOGO=calculate_links_evolution(adj_input_NOGO(:,:,high_perf_sessions_indices));

total_links{1,count_time_interval}{1,count_i_mice}.expert.CR= mean(network_links_NOGO{1,1},1,'omitnan');
shared_links{1,count_time_interval}{1,count_i_mice}.expert.CR= mean(network_links_NOGO{1,2},1,'omitnan');
network_links_shared_p_value{1,count_time_interval}{1,count_i_mice}.expert.CR= mean(network_links_NOGO{1,3},1,'omitnan');
% network_transitivity{1,count_time_interval}{1,count_i_mice}.expert.CR= mean(network_links_NOGO{1,4},1,'omitnan');
% network_assortativity{1,count_time_interval}{1,count_i_mice}.expert.CR=mean(network_links_NOGO{1,5},1,'omitnan');
% network_modularity{1,count_time_interval}{1,count_i_mice}.expert.CR=mean(network_links_NOGO{1,6},1,'omitnan');


% --- calculate network properties for NAIVE sessoins ---
% --- Hit trials ---
if isempty(low_perf_sessions_indices)
total_links{1,count_time_interval}{1,count_i_mice}.naive.Hit= NaN;
shared_links{1,count_time_interval}{1,count_i_mice}.naive.Hit= NaN;
network_links_shared_p_value{1,count_time_interval}{1,count_i_mice}.naive.Hit= NaN;
% network_transitivity{1,count_time_interval}{1,count_i_mice}.naive.Hit= NaN;
% network_assortativity{1,count_time_interval}{1,count_i_mice}.naive.Hit=NaN;
% network_modularity{1,count_time_interval}{1,count_i_mice}.naive.Hit=NaN;

% --- CR trials ---
total_links{1,count_time_interval}{1,count_i_mice}.naive.CR= NaN;
shared_links{1,count_time_interval}{1,count_i_mice}.naive.CR= NaN;
network_links_shared_p_value{1,count_time_interval}{1,count_i_mice}.naive.CR= NaN;
% network_transitivity{1,count_time_interval}{1,count_i_mice}.naive.CR= NaN;
% network_assortativity{1,count_time_interval}{1,count_i_mice}.naive.CR=NaN;
% network_modularity{1,count_time_interval}{1,count_i_mice}.naive.CR=NaN;    
else
    
% --- Hit trials ---    
network_links_GO=calculate_links_evolution(adj_input_GO(:,:,low_perf_sessions_indices));

total_links{1,count_time_interval}{1,count_i_mice}.naive.Hit= mean(network_links_GO{1,1},1,'omitnan');
shared_links{1,count_time_interval}{1,count_i_mice}.naive.Hit= mean(network_links_GO{1,2},1,'omitnan');
network_links_shared_p_value{1,count_time_interval}{1,count_i_mice}.naive.Hit= mean(network_links_GO{1,3},1,'omitnan');
% network_transitivity{1,count_time_interval}{1,count_i_mice}.naive.Hit= mean(network_links_GO{1,4},1,'omitnan');
% network_assortativity{1,count_time_interval}{1,count_i_mice}.naive.Hit=mean(network_links_GO{1,5},1,'omitnan');
% network_modularity{1,count_time_interval}{1,count_i_mice}.naive.Hit=mean(network_links_GO{1,6},1,'omitnan');

% --- CR trials ---
network_links_NOGO=calculate_links_evolution(adj_input_NOGO(:,:,low_perf_sessions_indices));

total_links{1,count_time_interval}{1,count_i_mice}.naive.CR= mean(network_links_NOGO{1,1},1,'omitnan');
shared_links{1,count_time_interval}{1,count_i_mice}.naive.CR= mean(network_links_NOGO{1,2},1,'omitnan');
network_links_shared_p_value{1,count_time_interval}{1,count_i_mice}.naive.CR= mean(network_links_NOGO{1,3},1,'omitnan');
% network_transitivity{1,count_time_interval}{1,count_i_mice}.naive.CR= mean(network_links_NOGO{1,4},1,'omitnan');
% network_assortativity{1,count_time_interval}{1,count_i_mice}.naive.CR=mean(network_links_NOGO{1,5},1,'omitnan');
% network_modularity{1,count_time_interval}{1,count_i_mice}.naive.CR=mean(network_links_NOGO{1,6},1,'omitnan');
end
    end
    
end
    

%% plot Total links shaded error bar

input=[];
for j=1:length(network_time_intervals)
    for i=1:mice
        input.expert.Hit(j,i)= total_links{1,j}{1,i}.expert.Hit;
        input.naive.Hit(j,i)= total_links{1,j}{1,i}.naive.Hit;
        input.expert.CR(j,i)= total_links{1,j}{1,i}.expert.CR;
        input.naive.CR(j,i)= total_links{1,j}{1,i}.naive.CR;
    end
end

time= cell2mat(network_time_intervals);

% --- subplot 1 ---
subplot(2,1,1)

input_to_plot=[];
input_to_plot.naive.mean=mean(input.naive.Hit,2,'omitnan');
input_to_plot.expert.mean=mean(input.expert.Hit,2,'omitnan');

input_to_plot.naive.std=std(input.naive.Hit,[],2,'omitnan')./sqrt(mice);
input_to_plot.expert.std=std(input.expert.Hit,[],2,'omitnan')./sqrt(mice);

shadedErrorBar(time(:,1),input_to_plot.naive.mean,input_to_plot.naive.std,'k',0)
hold on

shadedErrorBar(time(:,1),input_to_plot.expert.mean,input_to_plot.expert.std,'b',0);
xlim([1 8])
ylim([0 10])
title('Hit')

hold off

% --- subplot 2 ---
subplot(2,1,2)

input_to_plot=[];
input_to_plot.naive.mean=mean(input.naive.CR,2,'omitnan');
input_to_plot.expert.mean=mean(input.expert.CR,2,'omitnan');

input_to_plot.naive.std=std(input.naive.CR,[],2,'omitnan')./sqrt(mice);
input_to_plot.expert.std=std(input.expert.CR,[],2,'omitnan')./sqrt(mice);

shadedErrorBar(time(:,1),input_to_plot.naive.mean,input_to_plot.naive.std,'k',0)
hold on

shadedErrorBar(time(:,1),input_to_plot.expert.mean,input_to_plot.expert.std,'r',0);
xlim([1 8])
ylim([0 10])
title('CR')
hold off

%% plot Shared links shaded error bar

input=[];
for j=1:length(network_time_intervals)
    for i=1:mice
        input.expert.Hit(j,i)= shared_links{1,j}{1,i}.expert.Hit;
        input.naive.Hit(j,i)= shared_links{1,j}{1,i}.naive.Hit;
        input.expert.CR(j,i)= shared_links{1,j}{1,i}.expert.CR;
        input.naive.CR(j,i)= shared_links{1,j}{1,i}.naive.CR;
    end
end

time= cell2mat(network_time_intervals);

% --- subplot 1 ---
subplot(2,1,1)

input_to_plot=[];
input_to_plot.naive.mean=mean(input.naive.Hit,2,'omitnan');
input_to_plot.expert.mean=mean(input.expert.Hit,2,'omitnan');

input_to_plot.naive.std=std(input.naive.Hit,[],2,'omitnan')./sqrt(mice);
input_to_plot.expert.std=std(input.expert.Hit,[],2,'omitnan')./sqrt(mice);

shadedErrorBar(time(:,2),input_to_plot.naive.mean,input_to_plot.naive.std,'k',0)
hold on

shadedErrorBar(time(:,2),input_to_plot.expert.mean,input_to_plot.expert.std,'b',0);
xlim([0 8])
ylim([0 5])
title('Hit')

hold off

% --- subplot 2 ---
subplot(2,1,2)

input_to_plot=[];
input_to_plot.naive.mean=mean(input.naive.CR,2,'omitnan');
input_to_plot.expert.mean=mean(input.expert.CR,2,'omitnan');

input_to_plot.naive.std=std(input.naive.CR,[],2,'omitnan')./sqrt(mice);
input_to_plot.expert.std=std(input.expert.CR,[],2,'omitnan')./sqrt(mice);

shadedErrorBar(time(:,2),input_to_plot.naive.mean,input_to_plot.naive.std,'k',0)
hold on

shadedErrorBar(time(:,2),input_to_plot.expert.mean,input_to_plot.expert.std,'r',0);
xlim([0 8])
ylim([0 5])
title('CR')
hold off

% plot p-value of shared links shaded error bar

% input=[];
% for j=1:length(network_time_intervals)
%     for i=1:mice
%         input.expert.Hit(j,i)=network_links_shared_p_value{1,j}{1,i}.expert.Hit;
%         input.naive.Hit(j,i)= network_links_shared_p_value{1,j}{1,i}.naive.Hit;
%         input.expert.CR(j,i)= network_links_shared_p_value{1,j}{1,i}.expert.CR;
%         input.naive.CR(j,i)= network_links_shared_p_value{1,j}{1,i}.naive.CR;
%     end
% end
% 
% time= cell2mat(network_time_intervals);
% 
% % --- subplot 1 ---
% subplot(2,1,1)
% 
% input_to_plot=[];
% input_to_plot.naive.mean=mean(input.naive.Hit,2,'omitnan');
% input_to_plot.expert.mean=mean(input.expert.Hit,2,'omitnan');
% 
% input_to_plot.naive.std=std(input.naive.Hit,[],2,'omitnan')./sqrt(mice);
% input_to_plot.expert.std=std(input.expert.Hit,[],2,'omitnan')./sqrt(mice);
% 
% shadedErrorBar(time(:,2),input_to_plot.naive.mean,input_to_plot.naive.std,'k',0)
% hold on
% 
% shadedErrorBar(time(:,2),input_to_plot.expert.mean,input_to_plot.expert.std,'b',0);
% title('Hit')
% 
% hold off
% 
% % --- subplot 2 ---
% subplot(2,1,2)
% 
% input_to_plot=[];
% input_to_plot.naive.mean=mean(input.naive.CR,2,'omitnan');
% input_to_plot.expert.mean=mean(input.expert.CR,2,'omitnan');
% 
% input_to_plot.naive.std=std(input.naive.CR,[],2,'omitnan')./sqrt(mice);
% input_to_plot.expert.std=std(input.expert.CR,[],2,'omitnan')./sqrt(mice);
% 
% shadedErrorBar(time(:,2),input_to_plot.naive.mean,input_to_plot.naive.std,'k',0)
% hold on
% 
% shadedErrorBar(time(:,2),input_to_plot.expert.mean,input_to_plot.expert.std,'r',0);
% title('CR')
% hold off