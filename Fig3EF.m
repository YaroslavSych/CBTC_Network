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
    
for count_i_mice=1:length(pool_mice)

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
                    network_adj_matrices.GO.multivar(:,:,count(ii)) = network_bivar_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);

                elseif strcmp(trials_type{1,1},'CR')
                    
                    temp_data.bivar=data_TE{1,count_mice}{1,count(ii)}.bivar.CR.results;
                    network_adj_matrices.NOGO.bivar(:,:,count(ii)) = network_bivar_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.multivar=data_TE{1,count_mice}{1,count(ii)}.multivar.CR.results;
                    network_adj_matrices.NOGO.multivar(:,:,count(ii)) = network_bivar_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);
                                        
                    ii=ii+1;
                end
                end
                
                end
    end
end
   
% --- Choose high-performance sessions ---

    if perf_def == 1
        high_perf_sessions = [performance{1,count_mice}{1,count}] > network_performance_threshold;
        high_perf_sessions_indices = count(high_perf_sessions > 0);

    else
        high_perf_sessions = [performance{1,count_mice}{1,count}] < network_performance_threshold;
        high_perf_sessions_indices = count(high_perf_sessions > 0);

    end

% --- Calculate in-, out-degree and clustering coefficient ---
% if size(network_adj_matrices.GO.NIfTy,1)>12
%     channel_sub_network=[26 31 30 25 34 46 44 47 15 48 33 35];
% else
%     channel_sub_network=1:12;
% end

for count_channel_labels=1:comb
                    
    label = strncmp(channel_labels_all{count_mice},channel_label{1,count_channel_labels},5);
    channel_sub_network = find(label>0);
    if length(channel_sub_network)>1
        channel_sub_network=channel_sub_network(1,1);
    end
    
    if isempty(channel_sub_network)        
    else
        
    % compute outdegree of nodes
    % ---------- Hit ----------
    scale_to_network=length(channel_label);
    channel_labels_per_mouse=channel_labels_all{count_mice};

    % --- NIfTy ---
    % reorder adj matrix to a channel_label
%     temp_input= network_adj_matrices.GO.NIfTy(:,:,:);
%     output_aligned=reorder_channels_by_name(temp_input,high_perf_sessions_indices,channel_sub_network,channel_labels_per_mouse,channel_label);
%     outdegree_vector.Hit.NIfTy{1,count_channel_labels}{1,count_i_mice}=output_aligned;
%     
%     network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions_indices)>0,2));
%     
%     outdegree.mean_GO.NIfTy(count_channel_labels,count_i_mice) = mean(network_outdegree,1)/scale_to_network;
%     outdegree.var_GO.NIfTy(count_channel_labels,count_i_mice) = std(network_outdegree,0,1)/scale_to_network;

    % --- bivar ---
    temp_input= network_adj_matrices.GO.bivar(:,:,:);
    %output_aligned=reorder_channels_by_name(temp_input,high_perf_sessions_indices,channel_sub_network,channel_labels_per_mouse,channel_label);
    output_aligned=reorder_adj_matrix_by_channel_name(temp_input,high_perf_sessions_indices,count_channel_labels,channel_labels_per_mouse,channel_label);

    outdegree_vector.Hit.bivar{1,count_channel_labels}{1,count_i_mice}{1,count_time_interval}=output_aligned;
    
    
    % --- multivar ---
    temp_input= network_adj_matrices.GO.multivar(:,:,:);
    %output_aligned=reorder_channels_by_name(temp_input,high_perf_sessions_indices,channel_sub_network,channel_labels_per_mouse,channel_label);
    output_aligned=reorder_adj_matrix_by_channel_name(temp_input,high_perf_sessions_indices,count_channel_labels,channel_labels_per_mouse,channel_label);

    outdegree_vector.Hit.multivar{1,count_channel_labels}{1,count_i_mice}{1,count_time_interval}=output_aligned;
    
    % ---------- CR ---------- 
    % --- NIfTy ---
%     temp_input= network_adj_matrices.NOGO.NIfTy(:,:,:);
%     output_aligned=reorder_channels_by_name(temp_input,high_perf_sessions_indices,channel_sub_network,channel_labels_per_mouse,channel_label);
%     outdegree_vector.CR.NIfTy{1,count_channel_labels}{1,count_i_mice}=output_aligned;
%     
%     network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
% 
%     outdegree.mean_NOGO.NIfTy(count_channel_labels,count_i_mice) = mean(network_outdegree,1)/scale_to_network;
%     outdegree.var_NOGO.NIfTy(count_channel_labels,count_i_mice) = std(network_outdegree,0,1)/scale_to_network;
    
    % --- bivar ---
    temp_input= network_adj_matrices.NOGO.bivar(:,:,:);
    %output_aligned=reorder_channels_by_name(temp_input,high_perf_sessions_indices,channel_sub_network,channel_labels_per_mouse,channel_label);
    output_aligned=reorder_adj_matrix_by_channel_name(temp_input,high_perf_sessions_indices,count_channel_labels,channel_labels_per_mouse,channel_label);

    outdegree_vector.CR.bivar{1,count_channel_labels}{1,count_i_mice}{1,count_time_interval}=output_aligned;
    
    % --- multivar ---
    temp_input= network_adj_matrices.NOGO.multivar(:,:,:);
    %output_aligned=reorder_channels_by_name(temp_input,high_perf_sessions_indices,channel_sub_network,channel_labels_per_mouse,channel_label);
    output_aligned=reorder_adj_matrix_by_channel_name(temp_input,high_perf_sessions_indices,count_channel_labels,channel_labels_per_mouse,channel_label);

    outdegree_vector.CR.multivar{1,count_channel_labels}{1,count_i_mice}{1,count_time_interval}=output_aligned;

    end
    
end
    
end % count_mice
end % count_time_interval

%% plot p-value of all links

% ----------------- across mice ---------------------
% ------ averaged across input_time_intervals--------
channel_label={'GP', 'DG', 'VL','Cpu', 'CA1_Py', 'S1_bf','VM', 'LDVL', 'Rt','M1','CA1_Mol','dCA1'};

for count_time_interval = 1:length(network_time_intervals)
                
    for count_source=1:length(channel_label)

    for count_target=1:length(channel_label)
        
%     input_Hit=[];
%     input_CR=[];
    input_Hit=double.empty(0,12);
    input_CR=double.empty(0,12);
    
        if count_source~=count_target

            for count_mice=1:mice
                if isempty(outdegree_vector.Hit.multivar{1,count_source}{1,count_mice})
                else
                    
                    input= mean(outdegree_vector.Hit.multivar{1,count_source}{1,count_mice}{1,count_time_interval},1,'omitnan');
                    % change the size of the empty input to concatenate
                    if isempty(input)
                        input=double.empty(0,12);
                    end
                    input_Hit=cat(1, input_Hit,input);
                    
                    input= mean(outdegree_vector.CR.multivar{1,count_source}{1,count_mice}{1,count_time_interval},1,'omitnan');
                    % change the size of the empty input to concatenate
                    if isempty(input)
                        input=double.empty(0,12);
                    end
                    input_CR=cat(1, input_CR, input);
                end

            end

%         input_Hit(isnan(input_Hit))=0;
%         input_CR(isnan(input_CR))=0; 

        predictor_Hit{1,count_time_interval}{1,count_source}= input_Hit;
        predictor_CR{1,count_time_interval}{1,count_source}= input_CR;

        end

    end

    end
end

%% plot appearance frequency of each source
% -----------------------
% --- tested by mouse ---
% -----------------------

time= cell2mat(network_time_intervals);
time=time(2:end,1);
mean_link_frequency_Hit=[];
mean_link_frequency_CR=[];

% dimentions link_frequency_Hit(:,count_source,count_target,count_mice)

for count_time_interval = 2:length(network_time_intervals)-1
            for count_source=1:12
            for count_target=1:12 
                        
            input=[];
            input=predictor_Hit{1,count_time_interval}{1,count_source}(:,count_target);
            % mean across mice
            mean_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');
            % sem across mice
            sem_link_frequency_Hit(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1)); 
            
            %mean_link_frequency_Hit(count_time_interval,count_source,count_target)= adjacency_matrix_Hit{1,count_time_interval}{1,count_source}(:,count_target);
            
            input=[];
            input=predictor_CR{1,count_time_interval}{1,count_source}(:,count_target);
            mean_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');
            sem_link_frequency_CR(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1));
              %mean_link_frequency_CR(count_time_interval,count_source,count_target)= adjacency_matrix_CR{1,count_time_interval}{1,count_source}(:,count_target);

            end
            end
end


%% plot 90th percentile of the strongest links

channel_label= channel_labels_all{1, 7}; 
% network_time_interval
time_interval=15:17;
percent_to_threshold=90;

subplot(2,2,1)

% Hit networks
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=squeeze(mean(mean_link_frequency_Hit(time_interval,:,:),1,'omitnan'));
sem_to_plot=squeeze(mean(sem_link_frequency_Hit(time_interval,:,:),1,'omitnan'));

% find the 95th percentile
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);

% use percentile threshold from the texture_presentation time
%percentile_threshold=0.0546;

G = digraph(input_to_plot>percentile_threshold,channel_label)
centrality_values = centrality(G,'pagerank') % 'pagerank'

p=plot(G,'Layout','force')
p.MarkerSize = 10*(centrality_values+0.01);

title('90th percentile of strongest Sources Hit texture presentation 5FPS')

subplot(2,2,2)

% Hit networks
%input_to_fit=log(sort(reshape(input_to_plot,1,[]),'descend'));
%semilogy(exp(input_to_fit),'b','LineWidth',2)
% ylim([1e-3 1]);
% xlim([0 132]);

[input_sorted, ind_sort]=sort(reshape(input_to_plot,1,[]),'descend');
variance_sorted= reshape(sem_to_plot,1,[]);
shadedErrorBar((1:144),input_sorted,variance_sorted(ind_sort),'k');
ylim([0 0.2]);
xlim([0 132]);

hold on

hline(percentile_threshold)
hold off

subplot(2,2,3)

% CR networks

dim1=size(mean_link_frequency_CR,1);
input_to_plot=squeeze(mean(mean_link_frequency_CR(time_interval,:,:),1,'omitnan'));
sem_to_plot=squeeze(mean(sem_link_frequency_CR(time_interval,:,:),1,'omitnan'));

% find the 95th percentile
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);

% use percentile threshold from the texture_presentation time
%percentile_threshold=0.0492;

G = digraph(input_to_plot>percentile_threshold,channel_label)
centrality_values = centrality(G,'pagerank') % 'pagerank'

p=plot(G,'Layout','force')
p.MarkerSize = 20*(centrality_values+0.01);

title('95th percentile of strongest Sources Hit texture presentation 5FPS')

subplot(2,2,4)

% Hit networks
% input_to_fit=log(sort(reshape(input_to_plot,1,[]),'descend'));
% semilogy(exp(input_to_fit),'r','LineWidth',2)
[input_sorted, ind_sort]=sort(reshape(input_to_plot,1,[]),'descend');
variance_sorted= reshape(sem_to_plot,1,[]);
shadedErrorBar((1:144),input_sorted,variance_sorted(ind_sort),'k');
ylim([0 0.2]);
xlim([0 132]);

%ylim([1e-3 1]);
xlim([0 132]);


hold on

hline(percentile_threshold)
hold off