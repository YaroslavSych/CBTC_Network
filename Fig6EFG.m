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
perf_def =1;

% cannel labels for pooling
% get selected channels
label=[];
outdegree_vector=[];

channel_label =channel_labels_all{1, 8}';
comb=length(channel_label);
        
pool_mice=[8 11 12 13];


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
                
                % select timesteps in the desired time interval
                network_time_steps = count_time_interval;
                
                if length(network_time_steps) < 1
                    error('ERROR: no data points found in the specified time interval.')
                end
                
                % store the new adjacency matrix (assumption: all training sessions have the same number of channels)
                if strcmp(trials_type{1,1},'Hit')
                    
                    temp_data.bivar=data_TE{1,count_mice}{1,count(ii)}.bivar.Hit.results;
                    network_adj_matrices.GO.bivar(:,:,count(ii)) = network_any_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.multivar=data_TE{1,count_mice}{1,count(ii)}.multivar.Hit.results;
                    network_adj_matrices.GO.multivar(:,:,count(ii)) = network_any_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);

                elseif strcmp(trials_type{1,1},'CR')
                    
                    temp_data.bivar=data_TE{1,count_mice}{1,count(ii)}.bivar.CR.results;
                    network_adj_matrices.NOGO.bivar(:,:,count(ii)) = network_any_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.multivar=data_TE{1,count_mice}{1,count(ii)}.multivar.CR.results;
                    network_adj_matrices.NOGO.multivar(:,:,count(ii)) = network_any_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);
                                        
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


%%
mice=length(pool_mice);
% ----------------- across mice ---------------------
% ------ averaged across input_time_intervals--------

for count_time_interval = 1:length(network_time_intervals)
                
    for count_source=1:length(channel_label)

    for count_target=1:length(channel_label)
        
%     input_Hit=[];
%     input_CR=[];
    input_Hit=double.empty(0,48);
    input_CR=double.empty(0,48);
    input_Hit_sessions=double.empty(0,48);
    input_CR_sessions=double.empty(0,48);
    
        if count_source~=count_target

            for count_mice=1:mice
                if isempty(outdegree_vector.Hit.bivar{1,count_source}{1,count_mice})
                else
                    
                    % does averaging for sessions
                    input= mean(outdegree_vector.Hit.multivar{1,count_source}{1,count_mice}{1,count_time_interval},1,'omitnan');
                    % change the size of the empty input to concatenate
                    if isempty(input)
                        input=double.empty(0,48);
                    end
                    % concatenates across mice -> dimentions: [mice, target_channels]
                    input_Hit=cat(1, input_Hit,input);
                    input_Hit_sessions= cat(1, input_Hit_sessions,outdegree_vector.Hit.multivar{1,count_source}{1,count_mice}{1,count_time_interval});
                    
                    input= mean(outdegree_vector.CR.multivar{1,count_source}{1,count_mice}{1,count_time_interval},1,'omitnan');
                    % change the size of the empty input to concatenate
                    if isempty(input)
                        input=double.empty(0,48);
                    end
                    input_CR=cat(1, input_CR, input);
                    input_CR_sessions= cat(1, input_CR_sessions,outdegree_vector.CR.multivar{1,count_source}{1,count_mice}{1,count_time_interval});
                end

            end

%         input_Hit(isnan(input_Hit))=0;
%         input_CR(isnan(input_CR))=0; 

        predictor_Hit{1,count_time_interval}{1,count_source}= input_Hit;
        predictor_CR{1,count_time_interval}{1,count_source}= input_CR;

        output_sessions_Hit(count_time_interval,count_source,:,:)=input_Hit_sessions;
        output_sessions_CR(count_time_interval,count_source,:,:)=input_CR_sessions;

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

for count_time_interval = 1:length(network_time_intervals)
            for count_source=1:48
            for count_target=1:48
                        
            input=[];
            input=predictor_Hit{1,count_time_interval}{1,count_source}(:,count_target);
           
            % mean across mice
            mean_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');
            
            % split mean acr0ss mice into train and test to validate the sequence
            % permute mouse indeces and split between train and test
            mouse_ind=randperm(size(input,1));
            train_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input(mouse_ind(3:end)),1,'omitnan');
            test_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input(mouse_ind(1:2)),1,'omitnan');
            
            % sem across mice
            sem_link_frequency_Hit(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1)); 
            
            all_links_Hit=input;
            %mean_link_frequency_Hit(count_time_interval,count_source,count_target)= adjacency_matrix_Hit{1,count_time_interval}{1,count_source}(:,count_target);
            
            input=[];
            input=predictor_CR{1,count_time_interval}{1,count_source}(:,count_target);
            mean_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');
            
            % split mean acr0ss mice into train and test to validate the sequence
            % permute mouse indeces and split between train and test
            mouse_ind=randperm(size(input,1));
            train_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(3:end)),1,'omitnan');
            test_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(1:2)),1,'omitnan');
            
            sem_link_frequency_CR(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1));
            %mean_link_frequency_CR(count_time_interval,count_source,count_target)= adjacency_matrix_CR{1,count_time_interval}{1,count_source}(:,count_target);

            end
            end
end

%% plot adjacency matrices

cmin=-0.05; % for links reducing their strength during stimulus
cmax=0.05; % for links increasing their strength during stimulus
link_strength_axis=0.2;

% re-order input to the matrix
channels_ordered=[1:12 25:48 13:24];

cue_interval=5:7;
stimulus_interval=15:17;

% use cue as the baseline time interval
baseline=squeeze(mean(mean_link_frequency_Hit(cue_interval,channels_ordered,channels_ordered),1,'omitnan'));
subplot(2,3,1)

mean_link_strength_expert=squeeze(mean(mean_link_frequency_Hit(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'));

imagesc(mean_link_strength_expert-baseline)
% clusters
%cgo_all = clustergram(mean_link_strength_expert-baseline,'Colormap',redbluecmap,'Standardize','Row')

colormap(redbluecmap)
axis('square')
caxis([cmin cmax])
title('Hit expert')

subplot(2,3,2)

% remove unique names in channel labels
ind=[];
% make node names nonunique
for i=1:length(channel_labels)
    if sum(strcmp(channel_labels{i,1},channel_labels))>1
        % find indeces
        [~,~,ind]=find(strcmp(channel_labels{i,1},channel_labels)==1)
            channel_labels{i,1}=[channel_labels{i,1} '1'];
            channel_labels{ind(2),1}=[channel_labels{ind,1} '2'];
    end
end

% bottom 5%
bottom_threshold=prctile(reshape(mean_link_strength_expert-baseline,1,[]), 2);
% top 95%
top_threshold=prctile(reshape(mean_link_strength_expert-baseline,1,[]), 98);

input_adj_matrix=(mean_link_strength_expert-baseline>=top_threshold);
G = digraph(input_adj_matrix, channel_labels);
p=plot(G,'Layout','force')

title('strong')

subplot(2,3,3)

input_adj_matrix=(mean_link_strength_expert-baseline<=bottom_threshold);
G = digraph(input_adj_matrix, channel_labels);
p=plot(G,'Layout','force')

title('weak')

subplot(2,3,4)
% EXPERT

% use cue as the baseline time interval
baseline=squeeze(mean(mean_link_frequency_CR(cue_interval,channels_ordered,channels_ordered),1,'omitnan'));

mean_link_strength_expert=squeeze(mean(mean_link_frequency_CR(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'));

%[~, index_source_marginal] =sort(mean(mean_link_strength_expert,2,'omitnan'),1,'descend');
% [~, index_target] =sort(mean_link_strength,2,'descend');
% [~, index_source] =sort(mean_link_strength,1,'descend');

imagesc(mean_link_strength_expert-baseline)
colormap(redbluecmap)
axis('square')
caxis([cmin cmax])
title('CR expert')

subplot(2,3,5)

% bottom 5%
bottom_threshold=prctile(reshape(mean_link_strength_expert-baseline,1,[]), 2);
% top 95%
top_threshold=prctile(reshape(mean_link_strength_expert-baseline,1,[]), 98);

input_adj_matrix=(mean_link_strength_expert-baseline>=top_threshold);
G = digraph(input_adj_matrix, channel_labels);
p=plot(G,'Layout','force')

title('strong')

subplot(2,3,6)

input_adj_matrix=(mean_link_strength_expert-baseline<=bottom_threshold);
G = digraph(input_adj_matrix, channel_labels);
p=plot(G,'Layout','force')

title('weak')

%% test inputs and outputs across all channels for mean_link_frequency adjacency

p_threshold=0.01/48;
channels_ordered=[1:12 25:48 13:24];

cue_interval=5:7;
tex_interval=15:17;
% use cue as the baseline time interval
baseline_Hit=squeeze(mean(mean_link_frequency_Hit(cue_interval,channels_ordered,channels_ordered),1,'omitnan'));
baseline_CR=squeeze(mean(mean_link_frequency_CR(cue_interval,channels_ordered,channels_ordered),1,'omitnan'));

t=0.2*(6:size(mean_link_frequency_Hit,1)-1);
count=1;
all_sources=(1:48);
count=1;
p_val=[];

for count_interval=6:size(mean_link_frequency_Hit,1)
    
    % test for the input
    for count_source=1:48
    %stimulus_interval=count_interval-1:count_interval+1;
    stimulus_interval=count_interval;

    input_Hit=squeeze(mean(mean_link_frequency_Hit(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))-baseline_Hit;
    input_CR=squeeze(mean(mean_link_frequency_CR(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))-baseline_CR;
    %input_Hit=input_CR;
    % remove source from the matrix
    sources_in = all_sources(all_sources~=count_source);
    [p,h] = ranksum(input_Hit(count_source,:),reshape(input_Hit(sources_in,:),1,[]));
    % test if mean is above all in the adj_matrix or below
    if and(mean(input_Hit(count_source,:),2,'omitnan') < mean(reshape(input_Hit(sources_in,:),1,[]),2,'omitnan'),...
            p<=p_threshold)
        p_val(count_source,count)=-1;
    elseif and(mean(input_Hit(count_source,:),2,'omitnan') > mean(reshape(input_Hit(sources_in,:),1,[]),2,'omitnan'),...
            p<=p_threshold)
        p_val(count_source,count)=1;
    else
        p_val(count_source,count)=0;
    end
    end % coutn_source
    count=count+1;

end % count_interval

subplot(2,1,1)
imagesc(p_val)
sum(abs(p_val(:,tex_interval-6)),'all')

mymap=[0 0 1; 
    1 1 1;
    1 0 0;];
colormap(mymap)
% cmin=-2; % for regions with mean input below the adj matrix average
% cmax=2; % for regions with mean input above the adj matrix average
% caxis([cmin cmax])
title('for sources = outputs')
count=1;
p_val=[];

for count_interval=6:size(mean_link_frequency_Hit,1)
    
    % test for the input
    for count_target=1:48
    %stimulus_interval=count_interval-1:count_interval+1;
    stimulus_interval=count_interval;

    input_Hit=squeeze(mean(mean_link_frequency_Hit(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))-baseline_Hit;
    input_CR=squeeze(mean(mean_link_frequency_CR(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))-baseline_CR;
    %input_Hit=input_CR;

    % remove source from the matrix
    sources_in = all_sources(all_sources~=count_target);
    [p,h] = ranksum(input_Hit(:,count_target),reshape(input_Hit(:,sources_in),1,[])');
    % test if mean is above all in the adj_matrix or below
    if and(mean(input_Hit(:,count_target),1,'omitnan') < mean(reshape(input_Hit(:,sources_in),1,[]),2,'omitnan'),...
            p<=p_threshold)
        p_val(count_target,count)=-1;
    elseif and(mean(input_Hit(:,count_target),1,'omitnan') > mean(reshape(input_Hit(:,sources_in),1,[]),2,'omitnan'),...
            p<=p_threshold)
        p_val(count_target,count)=1;
    else
        p_val(count_target,count)=0;
    end
    end % coutn_source
        count=count+1;

end % count_interval

subplot(2,1,2)
imagesc(p_val)
colormap(mymap)
% cmin=-2; % for regions with mean input below the adj matrix average
% cmax=2; % for regions with mean input above the adj matrix average
% caxis([cmin cmax])
title('for targets = inputs')

sum(abs(p_val(:,tex_interval-6)),'all')
%% plot distributions for the 12th  channel
stimulus_interval=10:15;%from 2 to 3 sec of the trial time
count_target=12;
input_Hit=squeeze(mean(mean_link_frequency_CR(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))-baseline_CR;
   sources_in = all_sources(all_sources~=count_target);
   test_distribution=input_Hit(:,count_target);
   all_distribution=reshape(input_Hit(:,sources_in),1,[]);
   
edges = [-0.15:0.01:0.15];
subplot(2,1,1)
hist(test_distribution,edges);
h1 = findobj(gca,'Type','patch');
h1.FaceColor = [0.5 0.5 0.5];
hold on
vline(median(test_distribution,1,'omitnan'))
hold off

subplot(2,1,2)
hist(all_distribution,edges);
hold on
vline(median(all_distribution,2,'omitnan'))

hold off

%% do Hit - CR and test

t=0.2*(6:size(mean_link_frequency_Hit,1)-1);
p_threshold=0.01/48;

count=1;
all_sources=(1:48);
count=1;
p_val=[];

for count_interval=6:size(mean_link_frequency_Hit,1)-1
    
    % test for the input
    for count_source=1:48
    %stimulus_interval=count_interval-1:count_interval+1;
    stimulus_interval=count_interval;

    input_Hit=squeeze(mean(mean_link_frequency_Hit(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))- ...
        squeeze(mean(mean_link_frequency_CR(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'));

    % remove source from the matrix
    sources_in = all_sources(all_sources~=count_source);
    [p,h] = ranksum(input_Hit(count_source,:),reshape(input_Hit(sources_in,:),1,[]));
    % test if mean is above all in the adj_matrix or below
    if and(mean(input_Hit(count_source,:),2,'omitnan') < mean(reshape(input_Hit(sources_in,:),1,[]),2,'omitnan'),...
            p<=p_threshold)
        p_val(count_source,count)=-1;
    elseif and(mean(input_Hit(count_source,:),2,'omitnan') > mean(reshape(input_Hit(sources_in,:),1,[]),2,'omitnan'),...
            p<=p_threshold)
        p_val(count_source,count)=1;
    else
        p_val(count_source,count)=0;
    end
    end % coutn_source
    count=count+1;

end % count_interval

imagesc(p_val)
mymap=[0 0 1; 
    1 1 1;
    1 0 0;];
colormap(mymap)
% cmin=-2; % for regions with mean input below the adj matrix average
% cmax=2; % for regions with mean input above the adj matrix average
% caxis([cmin cmax])
title('for Hit - CR')
count=1;
p_val=[];


%% ---make a video for averaged adjacency matrices across SESSIONS for Hit and CR ---

cmin=-0.05; % for links reducing their strength during stimulus
cmax=0.05; % for links increasing their strength during stimulus
link_strength_axis=0.2;
p_val=[];
t_val=[];
count_significant_links=[];
cd 'C:\Users\yaros\Documents\papers\MultiFiber_Learning\Manuscript\Figures\fig6'
% re-order input to the matrix
channels_ordered=[1:12 25:48 13:24];

cue_interval=5:7;
% average across sessions in the baseline-cue period
session_cue_Hit= squeeze(mean(output_sessions_Hit,3,'omitnan'));
session_cue_CR= squeeze(mean(output_sessions_CR,3,'omitnan'));

% use time samples of the cue period and reorder channels according to
% 'a','m' and 'p' notation

baseline_Hit=squeeze(mean(session_cue_Hit(cue_interval,channels_ordered,channels_ordered),1,'omitnan'));
baseline_CR=squeeze(mean(session_cue_CR(cue_interval,channels_ordered,channels_ordered),1,'omitnan'));

t=0.2*(6:size(mean_link_frequency_Hit,1)-1);
count=1;
for count_interval=6:size(mean_link_frequency_Hit,1)-1

figure(1)     

% ---Hit ---
subplot(2,2,1)

stimulus_interval=count_interval:count_interval+1;

% average acrosss sessions
session_av_Hit= squeeze(mean(output_sessions_Hit,3,'omitnan'));
session_av_CR= squeeze(mean(output_sessions_CR,3,'omitnan'));

% average within a time interval -1*0.2s : current_time_sample :+1*0.2s
input_Hit=squeeze(mean(session_av_Hit(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))-baseline_Hit;
input_CR=squeeze(mean(session_av_CR(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))-baseline_CR;

imagesc(input_Hit)
% clusters
%cgo_all = clustergram(mean_link_strength_expert-baseline,'Colormap',redbluecmap,'Standardize','Row')

colormap(redbluecmap)
axis('square')
caxis([cmin cmax])
title('Hit Expert')
set(gca,'XTick',[], 'YTick', [])

% --- P values ---
subplot(2,2,[2 4])

% reshape to channels
input1= permute(output_sessions_Hit,[1 2 4 3]);
input2= permute(output_sessions_CR,[1 2 4 3]);

input_to_test_Hit= squeeze(mean(input1(stimulus_interval,:,:,:)));
input_to_test_CR= squeeze(mean(input2(stimulus_interval,:,:,:)));

% --- test across links ---
for count_source=1:48
    for count_target=1:48
        if count_source~=count_target
            [p,h] = ttest(reshape(input_to_test_Hit(count_source,count_target,:),1,[]),reshape(input_to_test_CR(count_source,count_target,:),1,[]));
            p_val(count_source,count_target)=p;
        else
            p_val(count_source,count_target)=NaN;
        end
    end
end
% threshold p_val
p_val(p_val<=0.01)=0;
p_val(p_val>0.01)=1;
p_val_reshaped=reshape(p_val,1,[]);

count_significant_links(count)= sum(p_val_reshaped==0);
t_val=[t_val t(count)];

plot(t_val,count_significant_links,'--k');

% --- test across groups of links ---
% --- choose anterior, medial and posterior groups ---
% regions_in_groups= [(1:12); (13:24); (25:36); (37:48)];
% 
% for count_group=1:4
%     [p,h] = ranksum(reshape(input_to_test_Hit(regions_in_groups(count_group,:),regions_in_groups(count_group,:),:),1,[]),...
%         reshape(input_to_test_CR(regions_in_groups(count_group,:),regions_in_groups(count_group,:),:),1,[]));
%     p_val(count,count_group)=p;
% end

% for count_group=1:4
%     semilogy(t_val,p_val(:,count_group),'LineWidth',2);
%     hold on
% end
%semilogy(t,p_val_naive,'--k')

hold on
hline(0.01);
hold on
vline(3);
hold on
vline(3.5);
hold off
xlim([1 10]);
ylim([10^-8 100]);
title(['Time' num2str(t(count)) ' s'])
ylabel('Log p');
xlabel('Trial time');
%legend('Expert', 'Naive', 'Location','southwest');
legend('Anterior', 'Medial1', 'Medial2','Posterior', 'Location','southwest');

% --- CR ---
subplot(2,2,3)

imagesc(input_CR)
colormap(redbluecmap)
axis('square')
caxis([cmin cmax])
title('CR Expert')
set(gcf,'color','w');
set(gca,'XTick',[], 'YTick', [])

F(count) = getframe(gcf);
count=count+1;

drawnow

end
  % create the video writer with 1 fps
  writerObj = VideoWriter('adj_matrix.avi');
  writerObj.FrameRate = 5;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

%% ---make a video for averaged adjacency matrices across MICE Hit and CR ---

cmin=-0.05; % for links reducing their strength during stimulus
cmax=0.05; % for links increasing their strength during stimulus
link_strength_axis=0.2;
p_val=[];
t_val=[];
cd 'C:\Users\yaros\Documents\papers\MultiFiber_Learning\Manuscript\Figures\fig6'
% re-order input to the matrix
channels_ordered=[1:12 25:48 13:24];

cue_interval=5:7;
% use cue as the baseline time interval
baseline_Hit=squeeze(mean(mean_link_frequency_Hit(cue_interval,channels_ordered,channels_ordered),1,'omitnan'));
baseline_CR=squeeze(mean(mean_link_frequency_CR(cue_interval,channels_ordered,channels_ordered),1,'omitnan'));

t=0.2*(6:size(mean_link_frequency_Hit,1)-1);
count=1;
for count_interval=6:size(mean_link_frequency_Hit,1)-1

figure(1)     

% ---Hit ---
subplot(2,2,1)

stimulus_interval=count_interval-1:count_interval+1;

input_Hit=squeeze(mean(mean_link_frequency_Hit(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))-baseline_Hit;
input_CR=squeeze(mean(mean_link_frequency_CR(stimulus_interval,channels_ordered,channels_ordered),1,'omitnan'))-baseline_CR;

imagesc(input_Hit)
% clusters
%cgo_all = clustergram(input_Hit,'Colormap',redbluecmap,'Standardize','Row')

colormap(redbluecmap)
axis('square')
caxis([cmin cmax])
title('Hit Expert')
set(gca,'XTick',[], 'YTick', [])

% --- P values ---
subplot(2,2,[2 4])
[p,h] = ranksum(reshape(input_Hit,1,[]),reshape(input_CR,1,[]));
p_val(count)=p;
t_val=[t_val t(count)]
semilogy(t_val,p_val,'LineWidth',2);
hold on
%semilogy(t,p_val_naive,'--k')
hold on
hline(0.01);
hold on
vline(3);
hold on
vline(3.5);
hold off
xlim([1 10]);
ylim([10^-60 1]);
title(['Time' num2str(t(count)) ' s'])
ylabel('Log p');
xlabel('Trial time');
legend('Expert', 'Naive', 'Location','southwest');

% --- CR ---
subplot(2,2,3)

imagesc(input_CR)
colormap(redbluecmap)
axis('square')
caxis([cmin cmax])
title('CR Expert')
set(gcf,'color','w');
set(gca,'XTick',[], 'YTick', [])

F(count) = getframe(gcf);
count=count+1;

drawnow

end
  % create the video writer with 1 fps
  writerObj = VideoWriter('adj_matrix.avi');
  writerObj.FrameRate = 5;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


%% plot network of sources for the top strongest links
% top 22 sources are defined according to a 99th percentile of links from
% the plot above where mean links strentgh adj matrix was thresholded 
top_sources=22*1;
sample_range=15:17;

subplot(2,2,1)

% Hit networks
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=reshape(mean_link_frequency_Hit,dim1,[]);
[sorted_input_to_plot, index_to_plot_Hit] =sort(mean(input_to_plot(sample_range,:),1,'omitnan'),'descend');

% labels by source area
input_source=reshape(ones(1,48).*(1:48)',1,[])';
%findgroups
%hist(input_source(index_to_plot_Hit(1:14)));
%channel_label(input_source(index_to_plot_Hit(1:top_sources)))

% labels by target area
input_target=reshape((1:48).*ones(48,1),1,[])';
%channel_label(input_target(index_to_plot_Hit(1:top_sources)))

G = digraph(channel_label(input_source(index_to_plot_Hit(1:top_sources))), channel_label(input_target(index_to_plot_Hit(1:top_sources))))
% add node GP 31 between nodes 1:12
%G = addedge(G,31,1:12);

centrality_values = centrality(G,'pagerank') % 'pagerank'

p=plot(G,'Layout','force')
p.MarkerSize = 20*(centrality_values+0.01);

title('strongest 14 Sources Hit texture presentation')

subplot(2,2,2)

% cumulative sum of top to bottom sorted links to the texture
% discrimination period
semilogy(cumsum(sorted_input_to_plot))
hold on
vline(top_sources)
hold off
ylim([0 25])

subplot(2,2,3)

% CR networks

dim1=size(mean_link_frequency_CR,1);
input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);
[sorted_input_to_plot, index_to_plot_CR] =sort(mean(input_to_plot(sample_range,:),1,'omitnan'),'descend');

% labels by source area
input_source=reshape(ones(1,48).*(1:48)',1,[])';
%findgroups
%hist(input_source(index_to_plot_CR(1:14)));
%channel_label(input_source(index_to_plot_CR(1:top_sources)))

% labels by target area
input_target=reshape((1:48).*ones(48,1),1,[])';
%channel_label(input_target(index_to_plot_CR(1:top_sources)))

G = digraph(channel_label(input_source(index_to_plot_CR(1:top_sources))), channel_label(input_target(index_to_plot_CR(1:top_sources))));
centrality_values = centrality(G,'pagerank') % 'pagerank'
p=plot(G,'Layout','force')
p.MarkerSize = 20*centrality_values+0.01;
title('strongest 14 Sources CR texture presentation')

subplot(2,2,4)

semilogy(cumsum(sorted_input_to_plot))
hold on
vline(top_sources)
hold off
ylim([0 25])

%% plot unwrapped connection strength across the trial time 
cmax=0.2;

figure
% Hit trials sorted with own index to plot
subplot(1,8,1)
dim1=size(train_link_frequency_Hit,1);
input_to_plot=reshape(train_link_frequency_Hit,dim1,[]);

% sort by the time of the peak in the binned
% bin with 5 samples
bin=5;
count_bin=1;
for i=1:5:45
    input_to_plot_bin(count_bin,:)= mean(input_to_plot(i:i+bin,:),1,'omitnan');
    count_bin=count_bin+1;
end

[~,ind]=max(input_to_plot_bin,[],1)
[~, index_to_plot_Hit] =sort(ind,'ascend');

% sort by the max of mean
%[~, index_to_plot_Hit] =sort(mean(input_to_plot(14:17,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_Hit(12:end))')
shading interp;
% 'Colormap',redbluecmap
c = 1-gray(20);
colormap(c)
% save my colormap
% ax=gca;
% mymap=colormap(ax);
% save('mycolormap','mymap')
% colormap(mymap)
xlim([1 50])
caxis([0 cmax])
title('mean Hit')

% ax = gca;
% mymap = colormap(ax);
% save('MyColormap3','mymap')

% Hit test plot from the split of the test/train by mice
subplot(1,8,2)

dim1=size(train_link_frequency_Hit,1);
input_to_plot=reshape(test_link_frequency_Hit,dim1,[]);
imagesc(input_to_plot(:,index_to_plot_Hit(12:end))')
shading interp;
% 'Colormap',redbluecmap
c = 1-gray(20);
colormap(c)
% save my colormap
% ax=gca;
% mymap=colormap(ax);
% save('mycolormap','mymap')
% colormap(mymap)
xlim([1 50])
caxis([0 cmax])
title('mean test Hit')

% CR trials sorted with own index to plot
subplot(1,8,3)
% dim1=size(mean_link_frequency_CR,1);
% input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);

dim1=size(train_link_frequency_CR,1);
input_to_plot=reshape(train_link_frequency_CR,dim1,[]);

% sort by the time of the peak
% sort by the time of the peak in the binned
% bin with 5 samples
bin=5;
count_bin=1;
for i=1:5:45
    input_to_plot_bin(count_bin,:)= mean(input_to_plot(i:i+bin,:),1,'omitnan');
    count_bin=count_bin+1;
end

[~,ind]=max(input_to_plot_bin,[],1)
[~, index_to_plot_CR] =sort(ind,'ascend');

%[~, index_to_plot_CR] =sort(mean(input_to_plot(15:19,:),1,'omitnan'),'descend');
imagesc(input_to_plot(:,index_to_plot_CR(12:end))')
xlim([1 50])
caxis([0 cmax])
shading interp;
title('mean CR')

% test CR trials
subplot(1,8,4)

dim1=size(test_link_frequency_CR,1);
input_to_plot=reshape(test_link_frequency_CR,dim1,[]);

imagesc(input_to_plot(:,index_to_plot_CR(12:end))')
xlim([1 50])
caxis([0 cmax])
shading interp;
title('test CR')

% labels by source area
subplot(1,8,5)
input_to_plot=reshape(ones(1,48).*(1:48)',1,[])';
imagesc(input_to_plot(index_to_plot_Hit(12:end)))
xlim([0.5 1])

% plot Hit vs CR strength
% subplot(1,6,3)
% dim1=size(mean_link_frequency_CR,1);
% input_to_plot1=reshape(mean_link_frequency_Hit,dim1,[]);
% input_to_plot2=reshape(mean_link_frequency_CR,dim1,[]);
% input_to_plot=input_to_plot1-input_to_plot2;
% 
% [~, index_to_plot_CR] =sort(mean(input_to_plot(15:19,:),1,'omitnan'),'descend');
% imagesc(input_to_plot(:,index_to_plot_Hit(2:end))')
% xlim([1 50])
% caxis([0 cmax])

% labels by target area
subplot(1,8,6)
input_to_plot=reshape((1:48).*ones(48,1),1,[])';
imagesc(input_to_plot(index_to_plot_Hit(12:end)))
xlim([0.5 1])

% Hit trials sorted with own index to plot
subplot(1,8,7)
dim1=size(sem_link_frequency_Hit,1);
input_to_plot=reshape(sem_link_frequency_Hit,dim1,[]);
[~, index_to_plot_Hit] =sort(mean(input_to_plot(15:19,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_Hit(12:end))')
shading interp;
% 'Colormap',redbluecmap
colormap(c)
xlim([1 50])
caxis([0 cmax])
title('sem Hit')
% ax = gca;
% mymap = colormap(ax);
% save('MyColormap3','mymap')

% CR trials sorted with own index to plot
subplot(1,8,8)
dim1=size(sem_link_frequency_CR,1);
input_to_plot=reshape(sem_link_frequency_CR,dim1,[]);
[~, index_to_plot_CR] =sort(mean(input_to_plot(15:19,:),1,'omitnan'),'descend');
imagesc(input_to_plot(:,index_to_plot_CR(12:end))')
xlim([1 50])
caxis([0 cmax])
title('sem CR')

%% plot projection of links from the test dataset
number_of_permutations=4;

for count_permutations = 1:number_of_permutations
    
    % for each permutation calculate train and test_link_frequency_Hit and
    % train and test_link_frequency_CR
    
    train_link_frequency_Hit=[];
    test_link_frequency_Hit=[];
    
    train_link_frequency_CR=[];
    test_link_frequency_CR=[];
    % dimentions link_frequency_Hit(:,count_source,count_target,count_mice)

    for count_time_interval = 1:length(network_time_intervals)
                for count_source=1:48
                for count_target=1:48

                input=[];
                input=predictor_Hit{1,count_time_interval}{1,count_source}(:,count_target);

                % mean across mice
                mean_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');

                % split mean acr0ss mice into train and test to validate the sequence
                % permute mouse indeces and split between train and test
                mouse_ind=randperm(size(input,1));
                train_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input(mouse_ind(3:end)),1,'omitnan');
                test_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input(mouse_ind(1:2)),1,'omitnan');

                % sem across mice
                sem_link_frequency_Hit(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1)); 

                all_links_Hit=input;
                %mean_link_frequency_Hit(count_time_interval,count_source,count_target)= adjacency_matrix_Hit{1,count_time_interval}{1,count_source}(:,count_target);

                input=[];
                input=predictor_CR{1,count_time_interval}{1,count_source}(:,count_target);
                mean_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');

                % split mean acr0ss mice into train and test to validate the sequence
                % permute mouse indeces and split between train and test
                mouse_ind=randperm(size(input,1));
                train_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(3:end)),1,'omitnan');
                test_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(1:2)),1,'omitnan');

                sem_link_frequency_CR(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1));
                  %mean_link_frequency_CR(count_time_interval,count_source,count_target)= adjacency_matrix_CR{1,count_time_interval}{1,count_source}(:,count_target);

                end
                end
    end

    % unwrap the link frequency
    train_to_plot=reshape(train_link_frequency_Hit,dim1,[]);

    % sort by the time of the peak in the binned
    % bin with 5 samples
    bin=5;
    count_bin=1;
    for i=1:5:45
        train_to_plot_bin(count_bin,:)= mean(train_to_plot(i:i+bin,:),1,'omitnan');
        count_bin=count_bin+1;
    end

    [~,ind]=max(train_to_plot_bin,[],1)
    [~, index_to_plot_Hit] =sort(ind,'ascend');

    % sort by the max of mean
    %[~, index_to_plot_Hit] =sort(mean(input_to_plot(14:17,:),1,'omitnan'),'descend');
    test_to_plot=reshape(test_link_frequency_Hit,dim1,[]);

    % calculate mean and SEM for multiple shuffles of train and test across
    % mice for the links sequence for evre block of 40 links
    delta_links=750;
    t=(1:50);
    input=[];
    
    for count_block=1:3
        input=train_to_plot(:,index_to_plot_Hit(12+delta_links*count_block-delta_links:12+delta_links*count_block));
        train_to_plot_mean(:,count_block,count_permutations)= mean(input,2,'omitnan');
        
        input=test_to_plot(:,index_to_plot_Hit(12+delta_links*count_block-delta_links:12+delta_links*count_block));
        test_to_plot_mean(:,count_block,count_permutations)= mean(input,2,'omitnan');
    end % coutn blocks of links

end % count permutations

delta_links=40;
t=(1:50);

% --- train_to_plot_mean contains averages of sorted links in train and test daa-sets
% --- across blocks_of_links for multiple permutations

for count_block=1:3
    
subplot(1,3,count_block)

input=mean(train_to_plot_mean(:,count_block,:),3,'omitnan');
input_sem=std(train_to_plot_mean(:,count_block,:),[],3,'omitnan')/sqrt(number_of_permutations);

shadedErrorBar(t,input,input_sem,'k')
hold on

input=mean(test_to_plot_mean(:,count_block,:),3,'omitnan');
input_sem=std(test_to_plot_mean(:,count_block,:),[],3,'omitnan')/sqrt(number_of_permutations);

shadedErrorBar(t,input,input_sem,'m')
ylim([0 0.1])
xlim([6 50])
end
title('average for Hit')

%% plot appearance frequency across the trial - "network unwrapped"
% --- mice pooled as observations ---

clear out_deg_Hit out_deg_CR
input_Hit=[];
input_CR=[];
adjacency_matrix_Hit=[];
adjacency_matrix_CR=[];
% --- do moving average-like link probability ---
% --- by moving time_interval across network_time_intervals
% --- here averaged across mice, not sessions ---

for count_time_interval = 2:length(network_time_intervals)-1
    for count_source=1:48
        
        for count_mice=1:mice
            
            if isempty(outdegree_vector.Hit.multivar{1,count_source}{1,count_mice})

            else
                
            time_interval= count_time_interval-1:1:count_time_interval+1;
            input_to_reshape_Hit=[];
            input_to_reshape_CR=[];
            
            for step=1:length(time_interval)
                
                input=outdegree_vector.Hit.multivar{1,count_source}{1,count_mice}{1,time_interval(step)};
                % average across sessions and concatenate for 3 timesteps
                input_to_reshape_Hit=cat(3,input_to_reshape_Hit,sum(input,1,'omitnan')/size(input,1));
                
                input=outdegree_vector.CR.multivar{1,count_source}{1,count_mice}{1,time_interval(step)};
                input_to_reshape_CR=cat(3,input_to_reshape_CR,sum(input,1,'omitnan')/size(input,1));
            end
                % average across timesteps and concatenate across mice
                input_Hit=cat(1, input_Hit, sum(input_to_reshape_Hit,3,'omitnan')/length(time_interval));
                input_CR=cat(1, input_CR, sum(input_to_reshape_CR,3,'omitnan')/length(time_interval));
            end            
        end
        
    out_deg_Hit{1,count_time_interval}{1,count_source}=input_Hit;
    out_deg_CR{1,count_time_interval}{1,count_source}=input_CR;

    adjacency_matrix_Hit{1,count_time_interval}{1,count_source}=mean(input_Hit,1,'omitnan');
    adjacency_matrix_CR{1,count_time_interval}{1,count_source}=mean(input_CR,1,'omitnan');
    
    input_Hit=[];
    input_CR=[];
    
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
            for count_source=1:48
            for count_target=1:48
                        
            input=[];
            input=out_deg_Hit{1,count_time_interval}{1,count_source}(:,count_target);
            % mean across mice
            mean_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');
            % sem across mice
            sem_link_frequency_Hit(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1)); 
            
            %mean_link_frequency_Hit(count_time_interval,count_source,count_target)= adjacency_matrix_Hit{1,count_time_interval}{1,count_source}(:,count_target);
            
            input=[];
            input=out_deg_CR{1,count_time_interval}{1,count_source}(:,count_target);
            mean_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');
            sem_link_frequency_CR(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1));
            %mean_link_frequency_CR(count_time_interval,count_source,count_target)= adjacency_matrix_CR{1,count_time_interval}{1,count_source}(:,count_target);

            end
            end
end

%% plot 99th percentile of the strongest links

ind=[];
% make node names nonunique
for i=1:length(channel_labels)
    if sum(strcmp(channel_labels{i,1},channel_labels))>1
        % find indeces
        [~,~,ind]=find(strcmp(channel_labels{i,1},channel_labels)==1)
            channel_labels{i,1}=[channel_labels{i,1} '1'];
            channel_labels{ind(2),1}=[channel_labels{ind,1} '2'];
    end
end

% network_time_interval
time_interval=15:17;
percent_to_threshold=99;

subplot(2,2,1)

% Hit networks
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=squeeze(mean(mean_link_frequency_Hit(time_interval,:,:),1,'omitnan'));
% find the 95th percentile
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);

% use percentile threshold from texture_time_interval 15:17;
%percentile_threshold=0.0994;

G = digraph(input_to_plot>percentile_threshold,channel_labels)
centrality_values = centrality(G,'pagerank') % 'pagerank'

p=plot(G,'Layout','force')
p.MarkerSize = 10*(centrality_values+0.01);

title('95th percentile of strongest Sources Hit texture presentation')

subplot(2,2,2)

% cumulative sum of top to bottom sorted links to the texture
% discrimination period
% h1=histfit(reshape(input_to_plot,1,[]),50,'exponential')
% h1(1).FaceColor=[.8 .8 1];
% h1(2).Color=[.2 .2 1];

input_to_fit=log(sort(reshape(input_to_plot,1,[]),'descend'));
semilogy(exp(input_to_fit),'b','LineWidth',2)
%plot(input_to_fit,'b','LineWidth',2)

hold on
%fitexp and truncate at zeros

line_space=1:sum(input_to_fit~=-Inf);
line_fit=polyfit(line_space,input_to_fit(1:sum(input_to_fit~=-Inf)),1);

line_space_plot=1:length(input_to_fit);
f = polyval(line_fit,line_space_plot);
plot(line_space_plot,exp(f),'--');
xlim([0 132])
hold on

hline(percentile_threshold)
hold off

subplot(2,2,3)

% CR networks

dim1=size(mean_link_frequency_CR,1);
input_to_plot=squeeze(mean(mean_link_frequency_CR(time_interval,:,:),1,'omitnan'));
% find the 99th percentile
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);

% use percentile threshold from texture_time_interval 15:17;
%percentile_threshold=0.1111;

G = digraph(input_to_plot>percentile_threshold,channel_labels)
centrality_values = centrality(G,'pagerank') % 'pagerank'

p=plot(G,'Layout','force')
p.MarkerSize = 20*(centrality_values+0.01);

title('95th percentile of strongest Sources Hit texture presentation')

subplot(2,2,4)

% cumulative sum of top to bottom sorted links to the texture
% discrimination period
% h1=histfit(reshape(input_to_plot,1,[]),50,'exponential')
% h1(1).FaceColor=[.8 .8 1];
% h1(2).Color=[.2 .2 1];
input_to_fit=log(sort(reshape(input_to_plot,1,[]),'descend'));
semilogy(exp(input_to_fit),'b','LineWidth',2)
%plot(input_to_fit,'b','LineWidth',2)

hold on
%fitexp and truncate at zeros

line_space=1:sum(input_to_fit~=-Inf);
line_fit=polyfit(line_space,input_to_fit(1:sum(input_to_fit~=-Inf)),1);

line_space_plot=1:length(input_to_fit);
f = polyval(line_fit,line_space_plot);
plot(line_space_plot,exp(f),'--');
xlim([0 132])
hold on

hline(percentile_threshold)
hold off


%%
cmax=0.1;
% sort by mean in a time iterval
time_interval=15:17;

figure
% Hit trials sorted with own index to plot
subplot(1,6,1)
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=reshape(mean_link_frequency_Hit,dim1,[]);

% sort by the time of the peak
% [~,ind]=max(input_to_plot,[],1)
% [~, index_to_plot_Hit] =sort(ind,'ascend');

% sort by the max of mean
[~, index_to_plot_Hit] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_Hit)')
% 'Colormap',redbluecmap
c = 1-gray(20);
colormap(c)
% save my colormap
% ax=gca;
% mymap=colormap(ax);
% save('mycolormap','mymap')
% colormap(mymap)
xlim([1 50])
caxis([0 cmax])
title('mean Hit')
%colormapeditor

% ax = gca;
% mymap = colormap(ax);
% save('MyColormap3','mymap')

% CR trials sorted with own index to plot
subplot(1,6,2)
dim1=size(mean_link_frequency_CR,1);
input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);
[~, index_to_plot_CR] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');
imagesc(input_to_plot(:,index_to_plot_Hit)')
xlim([1 50])
caxis([0 cmax])

% CR trials sorted with own index to plot
subplot(1,6,3)
dim1=size(mean_link_frequency_CR,1);
input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);
[~, index_to_plot_CR] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');
imagesc(input_to_plot(:,index_to_plot_CR)')
xlim([1 50])
caxis([0 cmax])

% Hit trials sorted with own index to plot
subplot(1,6,5)
dim1=size(sem_link_frequency_Hit,1);
input_to_plot=reshape(sem_link_frequency_Hit,dim1,[]);
%[~, index_to_plot_Hit] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_Hit)')
shading interp;
% 'Colormap',redbluecmap
colormap(c)
xlim([1 50])
caxis([0 cmax])
title('sem Hit')
% ax = gca;
% mymap = colormap(ax);
% save('MyColormap3','mymap')

% CR trials sorted with own index to plot
subplot(1,6,6)
dim1=size(sem_link_frequency_CR,1);
input_to_plot=reshape(sem_link_frequency_CR,dim1,[]);
%[~, index_to_plot_CR] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');
imagesc(input_to_plot(:,index_to_plot_CR)')
xlim([1 50])
caxis([0 cmax])
title('sem CR')

%% plot network of sources for the top strongest links
% top 22 sources are defined according to a 99th percentile of links from
% the plot above where mean links strentgh adj matrix was thresholded 
top_sources=22*1;
sample_range=15:17;

subplot(2,2,1)

% Hit networks
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=reshape(mean_link_frequency_Hit,dim1,[]);
[sorted_input_to_plot, index_to_plot_Hit] =sort(mean(input_to_plot(sample_range,:),1,'omitnan'),'descend');

% labels by source area
input_source=reshape(ones(1,48).*(1:48)',1,[])';
%findgroups
%hist(input_source(index_to_plot_Hit(1:14)));
%channel_label(input_source(index_to_plot_Hit(1:top_sources)))

% labels by target area
input_target=reshape((1:48).*ones(48,1),1,[])';
%channel_label(input_target(index_to_plot_Hit(1:top_sources)))

G = digraph(channel_label(input_source(index_to_plot_Hit(1:top_sources))), channel_label(input_target(index_to_plot_Hit(1:top_sources))))
% add node GP 31 between nodes 1:12
%G = addedge(G,31,1:12);

centrality_values = centrality(G,'pagerank') % 'pagerank'

p=plot(G,'Layout','force')
p.MarkerSize = 20*(centrality_values+0.01);

title('strongest 14 Sources Hit texture presentation')

subplot(2,2,2)

% cumulative sum of top to bottom sorted links to the texture
% discrimination period
semilogy(cumsum(sorted_input_to_plot))
hold on
vline(top_sources)
hold off
ylim([0 25])

subplot(2,2,3)

% CR networks

dim1=size(mean_link_frequency_CR,1);
input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);
[sorted_input_to_plot, index_to_plot_CR] =sort(mean(input_to_plot(sample_range,:),1,'omitnan'),'descend');

% labels by source area
input_source=reshape(ones(1,48).*(1:48)',1,[])';
%findgroups
%hist(input_source(index_to_plot_CR(1:14)));
%channel_label(input_source(index_to_plot_CR(1:top_sources)))

% labels by target area
input_target=reshape((1:48).*ones(48,1),1,[])';
%channel_label(input_target(index_to_plot_CR(1:top_sources)))

G = digraph(channel_label(input_source(index_to_plot_CR(1:top_sources))), channel_label(input_target(index_to_plot_CR(1:top_sources))));
centrality_values = centrality(G,'pagerank') % 'pagerank'
p=plot(G,'Layout','force')
p.MarkerSize = 20*centrality_values+0.01;
title('strongest 14 Sources CR texture presentation')

subplot(2,2,4)

semilogy(cumsum(sorted_input_to_plot))
hold on
vline(top_sources)
hold off
ylim([0 25])

%% test with ranks if the Hit vs CR link frequency differences are not random

delta_links=[];

for count_time_interval = 2:length(network_time_intervals)-1
            for count_source=1:48
                %for count_mice=1:mice

                    delta_links(count_time_interval,count_source,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:)-out_deg_CR{1,count_time_interval}{1,count_source}(:,:);

                %end
            end
end
               
% 49x48x4x48 trial_time_interval x source x count_mice x target
trial_time=10:40;

% average link freq. across 3-3.6s time interval
input_to_rank1=squeeze(mean(delta_links(trial_time,:,:,:),1,'omitnan'));

% select links for the distribution as independent observations across trial time
% input_to_rank1=delta_links(trial_time,:,:,:);

% reorder as source x target x mice
input_to_rank2=permute(input_to_rank1,[2 1 3 ]);
input_to_rank3=reshape(input_to_rank2,mice,[]);

% input_to_rank2=permute(input_to_rank1,[3 2 4 1]);
% input_to_rank3=reshape(input_to_rank2,mice,[]);

% sort/ Rank 
[~,rank_input]=sort(input_to_rank3,1);
% find a sign for each Rank
sign_plus=[];
sign_minus=[];
for count_mice=1:mice
%     sign_plus(count_mice,:)=input_to_rank3(count_mice,rank_input(count_mice,:))>0;
%     sign_minus(count_mice,:)=input_to_rank3(count_mice,rank_input(count_mice,:))<0;
    sign_plus(count_mice,:)=input_to_rank3(count_mice,:)>0;
    sign_minus(count_mice,:)=input_to_rank3(count_mice,:)<0;
end
sign= sign_plus-sign_minus;
% modulus of the sum across mice for the product sign and rank
test_rank=sign.*rank_input;
%imagesc(test_rank) % for each mouse
%imagesc(abs(sum(sign.*rank_input,1))) % across mice
distrbution_of_ranks=(sum(sign.*rank_input,1));

% generate random distribution of ranks
rand_input=rand(mice,size(input_to_rank3,2)) - rand(mice,size(input_to_rank3,2));
[~,rand_rank_input]=sort(rand_input,1);
rand_sign_plus=rand_input>0;
rand_sign_minus=rand_input<0;
rand_sign= rand_sign_plus-rand_sign_minus;

rand_distrbution_of_ranks=(sum(rand_sign.*rand_rank_input,1));
[p,h] =ranksum(rand_distrbution_of_ranks,distrbution_of_ranks);

subplot(3,1,1)
h2=histfit(rand_distrbution_of_ranks,20,'normal')
h2(1).FaceColor=[.8 .8 .8];
h2(2).Color=[.2 .2 .2];
hold on

h1=histfit(distrbution_of_ranks,20,'normal')
h1(1).FaceColor=[.8 .8 1];
h1(2).Color=[.2 .2 1];
hold off

subplot(3,1,2)
plot(rand_distrbution_of_ranks)
subplot(3,1,3)
plot(distrbution_of_ranks)

%% test for the mean link strength

delta_links=[];

% true link delta
for count_time_interval = 2:length(network_time_intervals)-1
            for count_source=1:48
                %for count_mice=1:mice

                    delta_links(count_time_interval,count_source,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:)-out_deg_CR{1,count_time_interval}{1,count_source}(:,:);

                %end
            end
end

% 49x48x4x48 trial_time_interval x source x count_mice x target
trial_time=10:40;

% average link freq. across 3-3.6s time interval
input1=squeeze(mean(delta_links(trial_time,:,:,:),1,'omitnan'));

% select links for the distribution as independent observations across trial time
% input_to_rank1=delta_links(trial_time,:,:,:);

% reorder as source x target x mice
input2=permute(input1,[2 1 3 ]);
input_delta=reshape(input2,mice,[]);

distrbution_of_ranks=(sum(input_delta,1));

% generate random distribution of ranks
% generate shuffled sign
for count_mice=1:mice
    rand_var=rand(1,size(input_delta,2));
    sign_minus=rand_var<0.5
    sign_plus=rand_var>0.5

    rand_sign(count_mice,:)=sign_plus-sign_minus;
end

rand_input_delta= abs(input_delta).*rand_sign; 
rand_distrbution_of_ranks=(sum(rand_input_delta,1));
[p,h] =ranksum(rand_distrbution_of_ranks,distrbution_of_ranks);

subplot(3,1,1)
h2=histfit(rand_distrbution_of_ranks)
h2(1).FaceColor=[.8 .8 .8];
h2(2).Color=[.2 .2 .2];
hold on

h1=histfit(distrbution_of_ranks)
h1(1).FaceColor=[.8 .8 1];
h1(2).Color=[.2 .2 1];
hold off

subplot(3,1,2)
plot(rand_distrbution_of_ranks)
subplot(3,1,3)
plot(distrbution_of_ranks)

%% do a scatter plot for link strength

subplot(4,1,1)

linksHit=[];
linksCR=[];
input1=[];
input2=[];
trial_time=15:17;
%trial_time=20:22;

% all links
for count_time_interval = 2:length(network_time_intervals)-1
            for count_source=1:48
                %for count_mice=1:mice

                    linksHit(count_time_interval,count_source,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:);
                    linksCR(count_time_interval,count_source,:,:)=out_deg_CR{1,count_time_interval}{1,count_source}(:,:);

                %end
            end
end
% average link freq. across 3-3.6s time interval and mice
input1=squeeze(mean(mean(linksHit(trial_time,:,:,:),1,'omitnan'),3));
% average link freq. across 3-3.6s time interval
input2=squeeze(mean(mean(linksCR(trial_time,:,:,:),1,'omitnan'),3));
sz=3;
% grey ALL 
c=[0.8 0.8 0.8];
scatter(reshape(input2,1,[]),reshape(input1,1,[]),sz,c,'filled')
hold on

linksHit=[];
linksCR=[];
input1=[];
input2=[];
% all outgoing GP links
for count_time_interval = 2:length(network_time_intervals)-1
            for count_source=[31 42]
                %for count_mice=1:mice

                    linksHit(count_time_interval,count_source,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:);
                    linksCR(count_time_interval,count_source,:,:)=out_deg_CR{1,count_time_interval}{1,count_source}(:,:);

                %end
            end
end
% average link freq. across 3-3.6s time interval and mice
input1=squeeze(mean(mean(linksHit(trial_time,:,:,:),1,'omitnan'),3));
% average link freq. across 3-3.6s time interval
input2=squeeze(mean(mean(linksCR(trial_time,:,:,:),1,'omitnan'),3));
sz=20;
% light blue GP links
c=[0. 0.5 0.8];
scatter(reshape(input2,1,[]),reshape(input1,1,[]),sz,c,'filled')
hold on

linksHit=[];
linksCR=[];
input1=[];
input2=[];
% all outgoing VM links
for count_time_interval = 2:length(network_time_intervals)-1
            for count_source=34
                %for count_mice=1:mice

                    linksHit(count_time_interval,count_source,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:);
                    linksCR(count_time_interval,count_source,:,:)=out_deg_CR{1,count_time_interval}{1,count_source}(:,:);

                %end
            end
end
% average link freq. across 3-3.6s time interval and mice
input1=squeeze(mean(mean(linksHit(trial_time,:,:,:),1,'omitnan'),3));
% average link freq. across 3-3.6s time interval
input2=squeeze(mean(mean(linksCR(trial_time,:,:,:),1,'omitnan'),3));
sz=20;
% violet VM links
c=[0.8 0. 0.5];
scatter(reshape(input2,1,[]),reshape(input1,1,[]),sz,c,'filled')

linksHit=[];
linksCR=[];
% all outgoing PrL links
for count_time_interval = 2:length(network_time_intervals)-1
            for count_source=11
                %for count_mice=1:mice

                    linksHit(count_time_interval,count_source,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:);
                    linksCR(count_time_interval,count_source,:,:)=out_deg_CR{1,count_time_interval}{1,count_source}(:,:);

                %end
            end
end
% average link freq. across 3-3.6s time interval and mice
input1=squeeze(mean(mean(linksHit(trial_time,:,:,:),1,'omitnan'),3));
% average link freq. across 3-3.6s time interval
input2=squeeze(mean(mean(linksCR(trial_time,:,:,:),1,'omitnan'),3));
sz=20;
% light green PrL links
c=[0 0.8 0.];
scatter(reshape(input2,1,[]),reshape(input1,1,[]),sz,c,'filled')


hold on

linksHit=[];
linksCR=[];
input1=[];
input2=[];
% all outgoing MO links
for count_time_interval = 2:length(network_time_intervals)-1
            for count_source=12
                %for count_mice=1:mice

                    linksHit(count_time_interval,count_source,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:);
                    linksCR(count_time_interval,count_source,:,:)=out_deg_CR{1,count_time_interval}{1,count_source}(:,:);

                %end
            end
end
% average link freq. across 3-3.6s time interval and mice
input1=squeeze(mean(mean(linksHit(trial_time,:,:,:),1,'omitnan'),3));
% average link freq. across 3-3.6s time interval
input2=squeeze(mean(mean(linksCR(trial_time,:,:,:),1,'omitnan'),3));
sz=20;
% dark green MO links
c=[0 0.5 0.];
scatter(reshape(input2,1,[]),reshape(input1,1,[]),sz,c,'filled')
ylim([0 0.15])
xlim([0 0.15])
%axis('square')
hold on

% linksHit=[];
% linksCR=[];
% % all outgoing S1BF links
% for count_time_interval = 2:length(network_time_intervals)-1
%             for count_source=[25 28]
%                 %for count_mice=1:mice
% 
%                     linksHit(count_time_interval,count_source,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:);
%                     linksCR(count_time_interval,count_source,:,:)=out_deg_CR{1,count_time_interval}{1,count_source}(:,:);
% 
%                 %end
%             end
% end
% % average link freq. across 3-3.6s time interval and mice
% input1=squeeze(mean(mean(linksHit(trial_time,:,:,:),1,'omitnan'),3));
% % average link freq. across 3-3.6s time interval
% input2=squeeze(mean(mean(linksCR(trial_time,:,:,:),1,'omitnan'),3));
% sz=10;
% % pink S1BF links
% c=[1 0.5 1.];
% scatter(reshape(input2,1,[]),reshape(input1,1,[]),sz,c,'filled')
% ylim([0 0.15])
% xlim([0 0.15])
% axis('square')
hold off
title('mean link strength')


% GP link distribution
subplot(4,1,2)
delta_links=[];

% true link delta
count=1;
for count_source=[31 42]
            for count_time_interval = 2:length(network_time_intervals)-1
                %for count_mice=1:mice

                    delta_links(count_time_interval,count,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:)-out_deg_CR{1,count_time_interval}{1,count_source}(:,:);
                %end
            end
            count=count+1;

end

% average link freq. across 3-3.6s time interval
input1=squeeze(mean(delta_links(trial_time,:,:,:),1,'omitnan'));

% select links for the distribution as independent observations across trial time
% input_to_rank1=delta_links(trial_time,:,:,:);

% reorder as source x target x mice
input2=permute(input1,[2 1 3]);
input_delta=reshape(input2,mice,[]);

distrbution_of_GP=(sum(input_delta,1));
% compare random distribution of link strengths (called rand_distrbution_of_ranks)
% and a distrbution_of_link_strength
[p,h] =ranksum(rand_distrbution_of_ranks(1:length(distrbution_of_GP)),distrbution_of_GP);

h2=histfit(rand_distrbution_of_ranks(1:length(distrbution_of_GP)))
h2(1).FaceColor=[.8 .8 .8];
h2(2).Color=[.2 .2 .2];
hold on

h1=histfit(distrbution_of_GP)
h1(1).FaceColor=[.8 .8 1];
h1(2).Color=[.2 .2 1];

title(['mean link strength GP p=' num2str(p)])
xlim([-0.4 0.4])
hold off

% MO link distribution
subplot(4,1,3)
delta_links=[];

% true link delta
count=1;
for count_source=12
            for count_time_interval = 2:length(network_time_intervals)-1
                %for count_mice=1:mice

                    delta_links(count_time_interval,count,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:)-out_deg_CR{1,count_time_interval}{1,count_source}(:,:);
                %end
            end
            count=count+1;

end

% average link freq. across 3-3.6s time interval
input1=squeeze(mean(delta_links(trial_time,:,:,:),1,'omitnan'));

% select links for the distribution as independent observations across trial time
% input_to_rank1=delta_links(trial_time,:,:,:);

% reorder as source x target x mice
input2=permute(input1,[2 1 3 ]);
input_delta=reshape(input2,mice,[]);

distrbution_of_MO=(sum(input_delta,1));
% compare random distribution of link strengths (called rand_distrbution_of_ranks)
% and a distrbution_of_link_strength
[p,h] =ranksum(rand_distrbution_of_ranks(1:length(distrbution_of_MO)),distrbution_of_MO);

h2=histfit(rand_distrbution_of_ranks(1:length(distrbution_of_MO)))
h2(1).FaceColor=[.8 .8 .8];
h2(2).Color=[.2 .2 .2];
hold on

h1=histfit(distrbution_of_MO)
h1(1).FaceColor=[0 0.5 0.];
h1(2).Color=[0 0.5 0.];

title(['mean link strength MO p=' num2str(p)])
xlim([-0.4 0.4])
hold off

% VM link distribution
subplot(4,1,4)
delta_links=[];

% true link delta
count=1;
for count_source=34
            for count_time_interval = 2:length(network_time_intervals)-1
                %for count_mice=1:mice

                    delta_links(count_time_interval,count,:,:)=out_deg_Hit{1,count_time_interval}{1,count_source}(:,:)-out_deg_CR{1,count_time_interval}{1,count_source}(:,:);
                %end
            end
            count=count+1;

end

% average link freq. across 3-3.6s time interval
input1=squeeze(mean(delta_links(trial_time,:,:,:),1,'omitnan'));

% select links for the distribution as independent observations across trial time
% input_to_rank1=delta_links(trial_time,:,:,:);

% reorder as source x target x mice
input2=permute(input1,[2 1 3 ]);
input_delta=reshape(input2,mice,[]);

distrbution_of_PrL=(sum(input_delta,1));
% compare random distribution of link strengths (called rand_distrbution_of_ranks)
% and a distrbution_of_link_strength
[p,h] =ranksum(rand_distrbution_of_ranks(1:length(distrbution_of_PrL)),distrbution_of_PrL);

h2=histfit(rand_distrbution_of_ranks(1:length(distrbution_of_PrL)))
h2(1).FaceColor=[.8 .8 .8];
h2(2).Color=[.2 .2 .2];
hold on

h1=histfit(distrbution_of_PrL)
h1(1).FaceColor=[0.8 0. 0.5];
h1(2).Color=[0.8 0. 0.5];

title(['mean link strength VM p=' num2str(p)])
xlim([-0.4 0.4])
hold off

%% try PCA
linksHit=[];
linksCR=[];
number_of_shuffles=10;

% averaged across mice links
for count_time_interval = 2:length(network_time_intervals)-1
                for count_source=1:48
                    for count_target=1:48

                        linksHit(count_time_interval,count_source,count_target)=mean(out_deg_Hit{1,count_time_interval}{1,count_source}(:,count_target));
                        linksCR(count_time_interval,count_source,count_target)=mean(out_deg_CR{1,count_time_interval}{1,count_source}(:,count_target));
                        
                    end
                end
end

% first subtracted and then averaged across mice links
for count_time_interval = 2:length(network_time_intervals)-1
                for count_source=1:48
                    for count_target=1:48
                        for count_mice=1:size(out_deg_Hit{1,count_time_interval}{1,count_source}(:,count_target),1)

                        links_Hit_mice(count_time_interval,count_source,count_target,count_mice)=out_deg_Hit{1,count_time_interval}{1,count_source}(count_mice,count_target);
                        links_CR_mice(count_time_interval,count_source,count_target,count_mice)=out_deg_CR{1,count_time_interval}{1,count_source}(count_mice,count_target);
                        
                        delta_links_mice(count_time_interval,count_source,count_target,count_mice)=links_Hit_mice(count_time_interval,count_source,count_target,count_mice)-...
                            links_CR_mice(count_time_interval,count_source,count_target,count_mice);
                        
                            for count_shuffle=1:number_of_shuffles
                                delta_links_shuffle(count_time_interval,count_source,count_target,count_mice,count_shuffle)=(2*randi([0 1],1)-1)*links_Hit_mice(count_time_interval,count_source,count_target,count_mice)-...
                                    (2*randi([0 1],1)-1)*links_CR_mice(count_time_interval,count_source,count_target,count_mice);

                            end
                        end
                        
                    end
                end
end

%%
%plot a difference of Hit-CR connectivity matrices

subplot(4,1,1)
% plot connectivity matrix for Hit
input1=mean(links_Hit_mice,4);
trial_time=15:17;
input2=squeeze(mean(input1(trial_time,:,:),1));
imagesc(input2,[0 0.1]);
axis('square');
colormap(cool);
hold off

subplot(4,1,2)
% [plot connectivity matrix for CR
input1=mean(links_CR_mice,4);
trial_time=15:17;
input2=squeeze(mean(input1(trial_time,:,:),1));
imagesc(input2,[0 0.1]);
axis('square');
colormap(cool);

subplot(4,1,3)
% plot connectivity matrix for the difference 0- to + direction
input1=mean(delta_links_mice,4);
trial_time=15:17;
input2=squeeze(mean(input1(trial_time,:,:),1));
imagesc(input2,[0 0.1]);
axis('square');
colormap(cool);

subplot(4,1,4)
% plot connectivity matrix for the difference 0- to - direction
input1=mean(delta_links_mice,4);
trial_time=15:17;
input2=squeeze(mean(input1(trial_time,:,:),1));
imagesc(input2,[-0.1 0]);
axis('square');
colormap(flipud(cool))
%cmap = colormap(gca);

%% Plot trial dynamics with PCA 

% dimesions of input: count_time_interval x count_source x count_target x count_mice
% delta_links_mice
% links_Hit_mice
% links_CR_mice

% ---Do the PCA---
% for Hit
input_to_PCA=reshape(mean(links_Hit_mice,4),49,[]);
[coeff,score,latent,~,explained] = pca(input_to_PCA-mean(input_to_PCA));
score_hit=score;

% for CR
input_to_PCA=reshape(mean(links_CR_mice,4),49,[]);
[coeff,score,latent,~,explained] = pca(input_to_PCA-mean(input_to_PCA));
score_cr=score;

% covarianceMatrix = cov(input_to_PCA);
% [V,D] = eig(covarianceMatrix);

% top 3 principal component vectors
vectors_PCA=coeff(:,1:3);

input_mean_shuffled_scores=squeeze(mean(shuffled_scores,3));

% plot PCA1 vs PCA2
subplot(2,1,1)
%plot(score(:,1),score(:,2),'k','LineWidth',2)

% plot Hit
x = score_hit(:,1)';
y = score_hit(:,2)';
z = zeros(size(x));
col = (1:49);  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
colormap(jet)

hold on
    plot3(x,y,z,'o','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 .5 .5],...
    'MarkerIndices',[1 15 30])
text(x(1),y(1),z(1),'start')
text(x(15),y(15),z(15),'texture')
text(x(30),y(30),z(30),'reward')
hold on

% ---plot pca by mouse---

% for count_mice=1:length(pool_mice)
%     input_to_PCA=reshape(links_Hit_mice(:,:,:,count_mice),49,[]);
%     [coeff,score,latent,~,explained] = pca(input_to_PCA-mean(input_to_PCA));
%     %plot3(score(:,1),score(:,2),z)
%     x = score(:,1)';
%     y = score(:,2)';
%     z = zeros(size(x));
%     surface([x;x],[y;y],[z;z],[col;col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);
%     
%     hold on
%     plot3(score(:,1),score(:,2),z,'o','MarkerSize',3,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[0.5 .5 .5],...
%     'MarkerIndices',[1 15 30])
% text(x(1),y(1),z(1),'start')
% text(x(15),y(15),z(15),'texture')
% text(x(30),y(30),z(30),'reward')
% hold on
% end
hold off

% plot CR
x = score_cr(:,1)';
y = score_cr(:,2)';
z = zeros(size(x));
col = (1:49);  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
colormap(jet)

hold on
    plot3(x,y,z,'s','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 .5 .5],...
    'MarkerIndices',[1 15 30])
text(x(1),y(1),z(1),'start')
text(x(15),y(15),z(15),'texture')
text(x(30),y(30),z(30),'no-reward')
hold on

xlabel('pca 1')
ylabel('pca 2')
axis('square')



subplot(2,1,2)

%plot3(score(:,1),score(:,2),score(:,3),'k','LineWidth',2)
x = score_hit(:,1)';
y = score_hit(:,2)';
z =score_hit(:,3)';
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    colormap(jet)
    hold on
    
plot3(x,y,z,'o','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 .5 .5],...
    'MarkerIndices',[1 15 30])
text(x(1),y(1),z(1),'start')
text(x(15),y(15),z(15),'texture')
text(x(30),y(30),z(30),'reward')

x = score_cr(:,1)';
y = score_cr(:,2)';
z =score_cr(:,3)';
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    colormap(jet)
    hold on
    
plot3(x,y,z,'s','MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 .5 .5],...
    'MarkerIndices',[1 15 30])
text(x(1),y(1),z(1),'start')
text(x(15),y(15),z(15),'texture')
text(x(30),y(30),z(30),'no-reward')
    
xlabel('pca 1')
ylabel('pca 2')
zlabel('pca 3')

%--- plot pca by mouse---

% for count_mice=1:length(pool_mice)
%     input_to_PCA=reshape(links_CR_mice(:,:,:,count_mice),49,[]);
%     [coeff,score,latent,~,explained] = pca(input_to_PCA-mean(input_to_PCA));
%     %plot3(score(:,1),score(:,2),z)
%     x = score(:,1)';
%     y = score(:,2)';
%     z = score(:,3)';
%     surface([x;x],[y;y],[z;z],[col;col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);
%     
%     hold on
%     plot3(x,y,z,'o','MarkerSize',3,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[0.5 .5 .5],...
%     'MarkerIndices',[1 15 30])
% text(x(1),y(1),z(1),'start')
% text(x(15),y(15),z(15),'texture')
% text(x(30),y(30),z(30),'reward')
% hold on
% 
% end
hold off

%Alternatively, will plot the fraction of variance explained by each component.
% 100*latent./sum(latent) == explained

% figure
% bar(explained)

%% plot PCA for the diff of Hit vs CR vs time

% dimesions of input: observations x variables
input_to_PCA=reshape(mean(delta_links_mice,4),49,[]);

% ---Do the PCA---
[coeff,score,latent,~,explained] = pca(input_to_PCA-mean(input_to_PCA));

% covarianceMatrix = cov(input_to_PCA);
% [V,D] = eig(covarianceMatrix);

% top 3 principal component vectors
vectors_PCA=coeff(:,1:3);

% calculate averaged scores for PCA for shuffled sign
for count_shuffle=1:number_of_shuffles
    
    input_1shuffle=delta_links_shuffle(:,:,:,:,count_shuffle);
    input_shuffle=reshape(mean(input_1shuffle,4),49,[]);
    [coeff_shuffled,score_shuffled,latent_shuffled] = pca(input_shuffle-mean(input_shuffle));
    vectors_shuffle=coeff_shuffled(:,1:3);
    
    shuffled_scores(:,:,count_shuffle)=score_shuffled(:,1:5);
end

input_mean_shuffled_scores=squeeze(mean(shuffled_scores,3));
input_std_shuffled_scores=squeeze(std(shuffled_scores,[],3))/sqrt(number_of_shuffles);

% plot PCA1 vs PCA2
subplot(5,1,1)
shadedErrorBar((1:49),input_mean_shuffled_scores(:,1),input_std_shuffled_scores(:,1))
hold on
plot(score(:,1),'k','LineWidth',2)
hold off
%% plot average adjacency matrix with strongest source 

cmax=0.1;
texture_presentation=15:17;
subplot(2,3,1)
% EXPERT
mean_link_strength_expert=squeeze(mean(mean_link_frequency_Hit(texture_presentation,:,:),1,'omitnan'));
mean_link_strength_expert(mean_link_strength_expert==0)=nan;

[~, index_source_marginal] =sort(mean(mean_link_strength_expert,2,'omitnan'),1,'descend');
% [~, index_target] =sort(mean_link_strength,2,'descend');
% [~, index_source] =sort(mean_link_strength,1,'descend');

imagesc(mean_link_strength_expert(index_source_marginal,:))
colormap(redbluecmap)
axis('square')
caxis([0 cmax])
title('Hit expert')

% NAIVE
subplot(2,3,2)
mean_link_strength_naive=squeeze(mean(mean_link_frequency_Hit(texture_presentation,:,:),1,'omitnan'));
mean_link_strength_naive(mean_link_strength_naive==0)=nan;

imagesc(mean_link_strength_naive(index_source_marginal,:))
colormap(redbluecmap)
axis('square')
caxis([0 cmax])
title('Hit naive')

% marginal distribution
subplot(2,3,3)
input_to_plot_1=mean(mean_link_strength_naive,2,'omitnan');
input_to_plot_2=mean(mean_link_strength_expert,2,'omitnan');

bar([input_to_plot_1(index_source_marginal) input_to_plot_2(index_source_marginal)])
ylabel([channel_label{1,index_source_marginal}])
ylim([0 0.5])
camroll(-90)


subplot(2,3,4)
% EXPERT
%mean_link_frequency_CR(isnan(mean_link_frequency_CR))=0;
mean_link_strength_expert=squeeze(mean(mean_link_frequency_CR(texture_presentation,:,:),1,'omitnan'));
mean_link_strength_expert(mean_link_strength_expert==0)=nan;

[~, index_source_marginal] =sort(mean(mean_link_strength_expert,2,'omitnan'),1,'descend');
% [~, index_target] =sort(mean_link_strength,2,'descend');
% [~, index_source] =sort(mean_link_strength,1,'descend');

imagesc(mean_link_strength_expert(index_source_marginal,:))
colormap(redbluecmap)
axis('square')
caxis([0 cmax])
title('Hit expert')

% NAIVE
subplot(2,3,5)
mean_link_strength_naive=squeeze(mean(mean_link_frequency_CR(texture_presentation,:,:),1,'omitnan'));
mean_link_strength_naive(mean_link_strength_naive==0)=nan;

imagesc(mean_link_strength_naive(index_source_marginal,:))
colormap(redbluecmap)
axis('square')
caxis([0 cmax])
title('Hit naive')

% marginal distribution
subplot(2,3,6)
input_to_plot_1=mean(mean_link_strength_naive,2,'omitnan');
input_to_plot_2=mean(mean_link_strength_expert,2,'omitnan');

bar([input_to_plot_1(index_source_marginal) input_to_plot_2(index_source_marginal)])
ylabel([channel_label{1,index_source_marginal}])
ylim([0 0.5])
camroll(-90)