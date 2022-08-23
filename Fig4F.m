%% --------------------------------
% --- calculate for 12 channels ---
% ---------------------------------

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
perf_def = 1;

% cannel labels for pooling
% get selected channels
label=[];
channel_label = {'GP', 'Cpu',  'VL',  'VM', 'LDVL', 'Rt', 'M1',  'S1_bf', 'CA1_Mol', 'DG', 'CA1_Py', 'dCA1'};

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

for count_time_interval = 1:length(network_time_intervals)
            for count_source=1:12
            for count_target=1:12 
                        
            input=[];
            input=predictor_Hit{1,count_time_interval}{1,count_source}(:,count_target);
           
            % mean across mice
            mean_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');
            
            % split mean acr0ss mice into train and test to validate the sequence
            % permute mouse indeces and split between train and test
            mouse_ind=randperm(size(input,1));
            train_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input(mouse_ind(6:end)),1,'omitnan');
            test_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input(mouse_ind(1:6)),1,'omitnan');
            
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
            train_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(6:end)),1,'omitnan');
            test_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(1:6)),1,'omitnan');
            
            sem_link_frequency_CR(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1));
              %mean_link_frequency_CR(count_time_interval,count_source,count_target)= adjacency_matrix_CR{1,count_time_interval}{1,count_source}(:,count_target);

            end
            end
end

%% create averages across 3 time intervals across all mice

intervals=[5:7; 15:17; 30:32];
p_threshold=0.05;

p_cue=[];
p_stim=[];
p_act=[];

for count_source=1:12
    for count_target=1:12

        input_interval_Hit=[];
        input_interval_CR=[];
           if count_source~=count_target
               
           for count_intervals=1:3
           % define samples of cue, stimulus and action-reward window
           interval=intervals(count_intervals,:);    
              
           input_Hit=[];
           input_CR=[];
              for count_time_interval=1:length(interval)
                  
                  input_Hit=cat(1,input_Hit,...
                      predictor_Hit{1,interval(count_time_interval)}{1,count_source}(:,count_target));
                  input_CR=cat(1,input_CR,...
                      predictor_CR{1,interval(count_time_interval)}{1,count_source}(:,count_target));
              end
              % remove nan from input hit and cr
              input_Hit(isnan(input_Hit))=[];
              input_CR(isnan(input_CR))=[];
              
              input_interval_Hit(count_intervals,:)=input_Hit;
              input_interval_CR(count_intervals,:)=input_CR;
           end % count_intervals for cue, stimulus and action-reward window
           
           % test permutations of all windows
           % cue
           [p1,h]= signrank(input_interval_Hit(1,:),input_interval_Hit(2,:));
           % stimulus
           [p2,h]= signrank(input_interval_Hit(2,:),input_interval_Hit(3,:));
           % action reward
           [p3,h]= signrank(input_interval_Hit(1,:),input_interval_Hit(3,:));
           
           p_stim.Hit(count_source,count_target)=sum([p1 p2]< p_threshold);
           p_cue.Hit(count_source,count_target)=sum([p1 p3]< p_threshold);
           p_act.Hit(count_source,count_target)=sum([p2 p3]< p_threshold);
           
                      % test permutations of all windows
           % cue
           [p1,h]= signrank(input_interval_CR(1,:),input_interval_CR(2,:));
           % stimulus
           [p2,h]= signrank(input_interval_CR(2,:),input_interval_CR(3,:));
           % action reward
           [p3,h]= signrank(input_interval_CR(1,:),input_interval_CR(3,:));
           
           p_stim.CR(count_source,count_target)=sum([p1 p2]< p_threshold);
           p_cue.CR(count_source,count_target)=sum([p1 p3]< p_threshold);
           p_act.CR(count_source,count_target)=sum([p2 p3]< p_threshold);
           

    end 
end 

end

subplot(3,2,1)
%imagesc(p_cue.Hit)
% G=digraph(p_cue.Hit>0,channel_label);
% network_plot = plot(G,'EdgeColor','b','layout','force');
title('cue Hit')

subplot(3,2,2)
imagesc(p_cue.CR)
title('cue CR')

subplot(3,2,3)
imagesc(p_stim.Hit)
title('stim Hit')

subplot(3,2,4)
imagesc(p_stim.CR)
title('stim CR')

subplot(3,2,5)
imagesc(p_act.Hit)
title('act Hit')

subplot(3,2,6)
imagesc(p_act.CR)
title('act CR')

%% use indeces to do average vs trial time plot

mean_to_plot=[];
sem_to_plot=[];

t=(1:50)';
subplot(3,1,1)
[source,target]=find(p_cue.Hit>0);
% calculate averaged across classified groups of links
mean_to_plot= mean(mean(mean_link_frequency_Hit(:,source,target),3,'omitnan'),2,'omitnan');
sem_to_plot=mean(mean(sem_link_frequency_Hit(:,source,target),3,'omitnan'),2,'omitnan');
shadedErrorBar(t,mean_to_plot,sem_to_plot)
title('cue related')

subplot(3,1,2)
[source,target]=find(p_stim.Hit>0);
% calculate averaged across classified groups of links
mean_to_plot= mean(mean(mean_link_frequency_Hit(:,source,target),3,'omitnan'),2,'omitnan');
sem_to_plot=mean(mean(sem_link_frequency_Hit(:,source,target),3,'omitnan'),2,'omitnan');
shadedErrorBar(t,mean_to_plot,sem_to_plot)
title('stimulus related')

subplot(3,1,3)
[source,target]=find(p_act.Hit>0);
% calculate averaged across classified groups of links
mean_to_plot= mean(mean(mean_link_frequency_Hit(:,source,target),3,'omitnan'),2,'omitnan');
sem_to_plot=mean(mean(sem_link_frequency_Hit(:,source,target),3,'omitnan'),2,'omitnan');
shadedErrorBar(t,mean_to_plot,sem_to_plot)
title('action related')

%% ---plot top 75th percentile of the strongest links---

% ind=[];
% % make node names nonunique
% for i=1:length(channel_labels)
%     if sum(strcmp(channel_labels{i,1},channel_labels))>1
%         % find indeces
%         [~,~,ind]=find(strcmp(channel_labels{i,1},channel_labels)==1)
%             channel_labels{i,1}=[channel_labels{i,1} '1'];
%             channel_labels{ind,1}=[channel_labels{ind,1} '2'];
%     end
% end
% network_time_interval
time_interval=15:17;
percent_to_threshold=94;

subplot(2,2,1)

% Hit networks
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=squeeze(mean(mean_link_frequency_Hit(time_interval,:,:),1,'omitnan'));
% find the 95th percentile
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);

% use percentile threshold from the texture_presentation time
%percentile_threshold=0.0546;

G = digraph(input_to_plot>percentile_threshold,channel_label)
centrality_values = centrality(G,'pagerank') % 'pagerank'
top_channels.Hit= G.Edges;

p=plot(G,'Layout','force')
p.MarkerSize = 10*(centrality_values+0.01);

title('90th percentile of strongest Sources Hit texture presentation 5FPS')

subplot(2,2,2)

% Hit networks
input_to_fit=log(sort(reshape(input_to_plot,1,[]),'descend'));
semilogy(exp(input_to_fit),'b','LineWidth',2)
ylim([1e-3 1]);
xlim([0 132]);
% cumulative sum of top to bottom sorted links to the texture
% discrimination period
% h1=histfit(reshape(input_to_plot,1,[]),50,'exponential')
% h1(1).FaceColor=[.8 .8 1];
% h1(2).Color=[.2 .2 1];

%input_to_fit=sort(reshape(input_to_plot,1,[]),'descend');
%semilogy(input_to_fit,'b','LineWidth',2)


%plot(input_to_fit,'b','LineWidth',2)

hold on
%fitexp and truncate at zeros

% line_space=1:sum(input_to_fit~=-Inf);
% line_fit=polyfit(line_space,input_to_fit(1:sum(input_to_fit~=-Inf)),1);
% 
% line_space_plot=1:length(input_to_fit);
% f = polyval(line_fit,line_space_plot);
% plot(line_space_plot,exp(f),'--');
% xlim([0 132])
% hold on

hline(percentile_threshold)
hold off

subplot(2,2,3)

% CR networks

dim1=size(mean_link_frequency_CR,1);
input_to_plot=squeeze(mean(mean_link_frequency_CR(time_interval,:,:),1,'omitnan'));
% find the 95th percentile
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);

% use percentile threshold from the texture_presentation time
%percentile_threshold=0.0492;

G = digraph(input_to_plot>percentile_threshold,channel_label)
centrality_values = centrality(G,'pagerank') % 'pagerank'
top_channels.CR= G.Edges;

p=plot(G,'Layout','force')
p.MarkerSize = 20*(centrality_values+0.01);

title('95th percentile of strongest Sources Hit texture presentation 5FPS')

subplot(2,2,4)

% Hit networks
input_to_fit=log(sort(reshape(input_to_plot,1,[]),'descend'));
semilogy(exp(input_to_fit),'r','LineWidth',2)
ylim([1e-3 1]);
xlim([0 132]);

% cumulative sum of top to bottom sorted links to the texture
% discrimination period
% h1=histfit(reshape(input_to_plot,1,[]),50,'exponential')
% h1(1).FaceColor=[.8 .8 1];
% h1(2).Color=[.2 .2 1];

% input_to_fit=sort(reshape(input_to_plot,1,[]),'descend');
% semilogy(input_to_fit,'r','LineWidth',2)

hold on
%fitexp and truncate at zeros

% line_space=1:sum(input_to_fit~=-Inf);
% line_fit=polyfit(line_space,input_to_fit(1:sum(input_to_fit~=-Inf)),1);
% 
% line_space_plot=1:length(input_to_fit);
% f = polyval(line_fit,line_space_plot);
% plot(line_space_plot,exp(f),'--');
% xlim([0 132])
% hold on

hline(percentile_threshold)
hold off

%% ---plot bottom 25th percentile of links---

% ind=[];
% % make node names nonunique
% for i=1:length(channel_labels)
%     if sum(strcmp(channel_labels{i,1},channel_labels))>1
%         % find indeces
%         [~,~,ind]=find(strcmp(channel_labels{i,1},channel_labels)==1)
%             channel_labels{i,1}=[channel_labels{i,1} '1'];
%             channel_labels{ind,1}=[channel_labels{ind,1} '2'];
%     end
% end
% network_time_interval
time_interval=15:17;
percent_to_threshold=25;

subplot(2,2,1)

% Hit networks
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=squeeze(mean(mean_link_frequency_Hit(time_interval,:,:),1,'omitnan'));
% find the 95th percentile
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);

% use percentile threshold from the texture_presentation time
%percentile_threshold=0.0546;

G = digraph(input_to_plot<percentile_threshold,channel_label)
centrality_values = centrality(G,'pagerank') % 'pagerank'
bottom_channels.Hit= G.Edges;

p=plot(G,'Layout','force')
p.MarkerSize = 10*(centrality_values+0.01);

title('90th percentile of strongest Sources Hit texture presentation 5FPS')

subplot(2,2,2)

% Hit networks
input_to_fit=log(sort(reshape(input_to_plot,1,[]),'descend'));
semilogy(exp(input_to_fit),'b','LineWidth',2)
ylim([1e-3 1]);
xlim([0 132]);
% cumulative sum of top to bottom sorted links to the texture
% discrimination period
% h1=histfit(reshape(input_to_plot,1,[]),50,'exponential')
% h1(1).FaceColor=[.8 .8 1];
% h1(2).Color=[.2 .2 1];

%input_to_fit=sort(reshape(input_to_plot,1,[]),'descend');
%semilogy(input_to_fit,'b','LineWidth',2)


%plot(input_to_fit,'b','LineWidth',2)

hold on
%fitexp and truncate at zeros

% line_space=1:sum(input_to_fit~=-Inf);
% line_fit=polyfit(line_space,input_to_fit(1:sum(input_to_fit~=-Inf)),1);
% 
% line_space_plot=1:length(input_to_fit);
% f = polyval(line_fit,line_space_plot);
% plot(line_space_plot,exp(f),'--');
% xlim([0 132])
% hold on

hline(percentile_threshold)
hold off

subplot(2,2,3)

% CR networks

dim1=size(mean_link_frequency_CR,1);
input_to_plot=squeeze(mean(mean_link_frequency_CR(time_interval,:,:),1,'omitnan'));
% find the 95th percentile
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);

% use percentile threshold from the texture_presentation time
%percentile_threshold=0.0492;

G = digraph(input_to_plot<percentile_threshold,channel_label)
centrality_values = centrality(G,'pagerank') % 'pagerank'
bottom_channels.CR= G.Edges;

p=plot(G,'Layout','force')
p.MarkerSize = 20*(centrality_values+0.01);

title('95th percentile of strongest Sources Hit texture presentation 5FPS')

subplot(2,2,4)

% Hit networks
input_to_fit=log(sort(reshape(input_to_plot,1,[]),'descend'));
semilogy(exp(input_to_fit),'r','LineWidth',2)
ylim([1e-3 1]);
xlim([0 132]);

% cumulative sum of top to bottom sorted links to the texture
% discrimination period
% h1=histfit(reshape(input_to_plot,1,[]),50,'exponential')
% h1(1).FaceColor=[.8 .8 1];
% h1(2).Color=[.2 .2 1];

% input_to_fit=sort(reshape(input_to_plot,1,[]),'descend');
% semilogy(input_to_fit,'r','LineWidth',2)

hold on
%fitexp and truncate at zeros

% line_space=1:sum(input_to_fit~=-Inf);
% line_fit=polyfit(line_space,input_to_fit(1:sum(input_to_fit~=-Inf)),1);
% 
% line_space_plot=1:length(input_to_fit);
% f = polyval(line_fit,line_space_plot);
% plot(line_space_plot,exp(f),'--');
% xlim([0 132])
% hold on

hline(percentile_threshold)
hold off

%%
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

for i=1:size(input_to_plot,2)
    input_to_plot_smooth(:,i) = smooth(input_to_plot(:,i),'moving');
end
% use here input_to_plot_smooth instead of input_to_plot_bin
[~,ind]=max(input_to_plot_smooth,[],1)
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
for i=1:size(input_to_plot,2)
    input_to_plot_smooth(:,i) = smooth(input_to_plot(:,i),'moving');
end
% use here input_to_plot_smooth instead of input_to_plot_bin
[~,ind]=max(input_to_plot_smooth,[],1)
[~, index_to_plot_CR] =sort(ind,'ascend');

% sort by the max of mean
%[~, index_to_plot_Hit] =sort(mean(input_to_plot(14:17,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_CR(12:end))')
shading interp;

xlim([1 50])
caxis([0 cmax])
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
input_to_plot=reshape(ones(1,12).*(1:12)',1,[])';
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
input_to_plot=reshape((1:12).*ones(12,1),1,[])';
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
% ----------- for Hit trials -----------  
number_of_permutations=20;

for count_permutations = 1:number_of_permutations
    
    % for each permutation calculate train and test_link_frequency_Hit and
    % train and test_link_frequency_CR
    
    train_link_frequency_Hit=[];
    test_link_frequency_Hit=[];
    
    train_link_frequency_CR=[];
    test_link_frequency_CR=[];
    % dimentions link_frequency_Hit(:,count_source,count_target,count_mice)

    for count_time_interval = 1:length(network_time_intervals)
                for count_source=1:12
                for count_target=1:12 

                input=[];
                input=predictor_Hit{1,count_time_interval}{1,count_source}(:,count_target);

                % mean across mice
                mean_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');

                % split mean acr0ss mice into train and test to validate the sequence
                % permute mouse indeces and split between train and test
                mouse_ind=randperm(size(input,1));
                train_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input(mouse_ind(6:end)),1,'omitnan');
                test_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input(mouse_ind(1:6)),1,'omitnan');

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
                train_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(6:end)),1,'omitnan');
                test_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(1:6)),1,'omitnan');

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

    index_Hit_sorted(:,count_permutations)=index_to_plot_Hit;

    % sort by the max of mean
    %[~, index_to_plot_Hit] =sort(mean(input_to_plot(14:17,:),1,'omitnan'),'descend');
    test_to_plot=reshape(test_link_frequency_Hit,dim1,[]);

    % calculate mean and SEM for multiple shuffles of train and test across
    % mice for the links sequence for evre block of 40 links
    delta_links=40;
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


%% plot projection of links from the test dataset
% ----------- for Hit and CR trials -----------  
number_of_permutations=20;

for count_permutations = 1:number_of_permutations
    
    % for each permutation calculate train and test_link_frequency_Hit and
    % train and test_link_frequency_CR
    
    train_link_frequency_Hit=[];
    test_link_frequency_Hit=[];
    
    train_link_frequency_CR=[];
    test_link_frequency_CR=[];
    % dimentions link_frequency_Hit(:,count_source,count_target,count_mice)

for count_time_interval = 1:length(network_time_intervals)
                for count_source=1:12
                for count_target=1:12 

                input=[];
                input=predictor_Hit{1,count_time_interval}{1,count_source}(:,count_target);

                % mean across mice
                mean_link_frequency_Hit(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');

                % split mean acr0ss mice into train and test to validate the sequence
                % permute mouse indeces and split between train and test
                mouse_ind=randperm(size(input,1));
                train_link_frequency_Hit(count_time_interval,count_source,count_target,count_permutations)=mean(input(mouse_ind(6:end)),1,'omitnan');
                test_link_frequency_Hit(count_time_interval,count_source,count_target,count_permutations)=mean(input(mouse_ind(1:6)),1,'omitnan');

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
                train_link_frequency_CR(count_time_interval,count_source,count_target,count_permutations)=mean(input(mouse_ind(6:end)),1,'omitnan');
                test_link_frequency_CR(count_time_interval,count_source,count_target,count_permutations)=mean(input(mouse_ind(1:6)),1,'omitnan');

                sem_link_frequency_CR(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1));
                  %mean_link_frequency_CR(count_time_interval,count_source,count_target)= adjacency_matrix_CR{1,count_time_interval}{1,count_source}(:,count_target);

                end
                end
    end % count time interval
end % count permutations

% -------- use Hit as the input -----------------------
    % take mean across bootstrap (mouse permutations)and unwrap the link frequency
    train_to_plot=reshape(mean(train_link_frequency_Hit,4,'omitnan'),dim1,[]);

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
    test_to_plot=reshape(mean(test_link_frequency_Hit,4,'omitnan'),dim1,[]);

    cmax=0.01;

    figure
    % Hit trials sorted with own index to plot
    subplot(1,4,1)

    imagesc(train_to_plot')
    shading interp;
    % 'Colormap',redbluecmap
    c = 1-gray(20);
    colormap(c)
    xlim([1 50])
    caxis([0 cmax])
    title('train Hit')

    subplot(1,4,2)

    imagesc(test_to_plot(:,index_to_plot_Hit(12:end))')
    shading interp;
    c = 1-gray(20);
    colormap(c)
    xlim([1 50])
    caxis([0 cmax])
    title('mean test Hit')

% -------- use CR as the input -----------------------
    % take mean across bootstrap (mouse permutations)and unwrap the link frequency
    train_to_plot=reshape(mean(train_link_frequency_CR,4,'omitnan'),dim1,[]);

    % sort by the time of the peak in the binned
    % bin with 5 samples
    bin=5;
    count_bin=1;
    for i=1:5:45
        train_to_plot_bin(count_bin,:)= mean(train_to_plot(i:i+bin,:),1,'omitnan');
        count_bin=count_bin+1;
    end

    [~,ind]=max(train_to_plot_bin,[],1)
    [~, index_to_plot_CR] =sort(ind,'ascend');

    % sort by the max of mean
    %[~, index_to_plot_Hit] =sort(mean(input_to_plot(14:17,:),1,'omitnan'),'descend');
    test_to_plot=reshape(mean(test_link_frequency_CR,4,'omitnan'),dim1,[]);
    
    subplot(1,4,3)

    imagesc(train_to_plot')
    shading interp;
    % 'Colormap',redbluecmap
    c = 1-gray(20);
    colormap(c)
    xlim([1 50])
    caxis([0 cmax])
    title('train CR')

    subplot(1,4,4)

    imagesc(test_to_plot(:,index_to_plot_Hit(12:end))')
    shading interp;
    c = 1-gray(20);
    colormap(c)
    xlim([1 50])
    caxis([0 cmax])
    title('mean test CR')


%% plot projection of links from the test dataset
% ----------- for CR trials -----------  
number_of_permutations=20;

for count_permutations = 1:number_of_permutations
    
    % for each permutation calculate train and test_link_frequency_Hit and
    % train and test_link_frequency_CR
    
    train_link_frequency_CR=[];
    test_link_frequency_CR=[];
    % dimentions link_frequency_Hit(:,count_source,count_target,count_mice)

    for count_time_interval = 1:length(network_time_intervals)
                for count_source=1:12
                for count_target=1:12 

                input=[];
                input=predictor_CR{1,count_time_interval}{1,count_source}(:,count_target);
                mean_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input,1,'omitnan');

                % split mean acr0ss mice into train and test to validate the sequence
                % permute mouse indeces and split between train and test
                mouse_ind=randperm(size(input,1));
                train_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(6:end)),1,'omitnan');
                test_link_frequency_CR(count_time_interval,count_source,count_target)=mean(input(mouse_ind(1:6)),1,'omitnan');

                sem_link_frequency_CR(count_time_interval,count_source,count_target)=std(input,[],1,'omitnan')/sqrt(size(input,1));
                  %mean_link_frequency_CR(count_time_interval,count_source,count_target)= adjacency_matrix_CR{1,count_time_interval}{1,count_source}(:,count_target);

                end
                end
    end

    % unwrap the link frequency
    train_to_plot=reshape(train_link_frequency_CR,dim1,[]);

    % sort by the time of the peak in the binned
    % bin with 5 samples
    bin=5;
    count_bin=1;
    for i=1:5:45
        train_to_plot_bin(count_bin,:)= mean(train_to_plot(i:i+bin,:),1,'omitnan');
        count_bin=count_bin+1;
    end

    [~,ind]=max(train_to_plot_bin,[],1)
    [~, index_to_plot_CR] =sort(ind,'ascend');

    index_CR_sorted(:,count_permutations)=index_to_plot_CR;
    
    % sort by the max of mean
    %[~, index_to_plot_Hit] =sort(mean(input_to_plot(14:17,:),1,'omitnan'),'descend');
    test_to_plot=reshape(test_link_frequency_CR,dim1,[]);

    % calculate mean and SEM for multiple shuffles of train and test across
    % mice for the links sequence for evre block of 40 links
    delta_links=40;
    t=(1:50);
    input=[];
    
    for count_block=1:3
        input=train_to_plot(:,index_to_plot_CR(12+delta_links*count_block-delta_links:12+delta_links*count_block));
        train_to_plot_mean(:,count_block,count_permutations)= mean(input,2,'omitnan');
        
        input=test_to_plot(:,index_to_plot_CR(12+delta_links*count_block-delta_links:12+delta_links*count_block));
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

%% plot for train_to_plot and test_to_plot averaged across block of links 
% Here train_to_plot and test_to_plot are generated from one shuffle acros mice
delta_links=40;
t=(1:50);

for i=1:3
subplot(1,3,i)
input=train_to_plot(:,index_to_plot_Hit(12+delta_links*i-delta_links:12+delta_links*i));
shadedErrorBar(t,mean(input,2,'omitnan'),std(input,[],2,'omitnan')/sqrt(40),'k')
hold on
input=test_to_plot(:,index_to_plot_Hit(12+delta_links*i-delta_links:12+delta_links*i));
shadedErrorBar(t,mean(input,2,'omitnan'),std(input,[],2,'omitnan')/sqrt(40),'m')
ylim([0 0.1])
xlim([0 50])
end
title('average for Hit')

%% calculate probability that cue, stimulus and test links stay in the same block
delta_links=40;
input=[];
    
for count_block=1:3
    for count_links=1:144
        links_in_block= index_Hit_sorted(13+delta_links*count_block-delta_links:13+delta_links*count_block,:);
        block_stability_Hit(count_links,count_block)= sum(sum(links_in_block==count_links,2),1)/number_of_permutations;
        
        links_in_block= index_CR_sorted(13+delta_links*count_block-delta_links:13+delta_links*count_block,:);
        block_stability_CR(count_links,count_block)= sum(sum(links_in_block==count_links,2),1)/number_of_permutations;
    end
end

%imagesc(block_stability_CR)

stable_links.Hit=block_stability_Hit>0.75;
%sum(stable_links.Hit,1)
stable_links.CR=block_stability_CR>0.75;
%sum(stable_links.CR,1)

% rename ---stable_links_indices--- if the comparison for naive to expert is required 
stable_links_indices.Hit.cue=find(stable_links.Hit(:,1)>0);
stable_links_indices.Hit.tex=find(stable_links.Hit(:,2)>0);
stable_links_indices.Hit.rew=find(stable_links.Hit(:,3)>0);
stable_links_indices.CR.cue=find(stable_links.CR(:,1)>0);
stable_links_indices.CR.tex=find(stable_links.CR(:,2)>0);
stable_links_indices.CR.rew=find(stable_links.CR(:,3)>0);

%% compare stable links for naive vs expert conditions

% only 1 link per block is overlapping
cue_links=intersect(stable_links_indices.Hit.cue,stable_links_indices_expert.Hit.cue)
tex_links=intersect(stable_links_indices.Hit.tex,stable_links_indices_expert.Hit.tex)
rew_links=intersect(stable_links_indices.Hit.rew,stable_links_indices_expert.Hit.rew)

intersect(stable_links_indices.CR.cue,stable_links_indices_expert.CR.cue)
intersect(stable_links_indices.CR.tex,stable_links_indices_expert.CR.tex)
intersect(stable_links_indices.CR.rew,stable_links_indices_expert.CR.rew)

intersect(stable_links_indices_expert.Hit.cue,stable_links_indices_expert.CR.cue)
intersect(stable_links_indices_expert.Hit.tex,stable_links_indices_expert.CR.tex)
intersect(stable_links_indices_expert.Hit.rew,stable_links_indices_expert.CR.rew)

% recalculate from links to matrix element
% labels by source area
input_source=reshape(ones(1,12).*(1:12)',1,[]);
%findgroups
%hist(input_source(index_to_plot_Hit(1:14)));
%channel_label(input_source(index_to_plot_Hit(1:top_sources)))

% labels by target area
input_target=reshape((1:12).*ones(12,1),1,[]);
%channel_label(input_target(index_to_plot_Hit(1:top_sources)))
template_source_to_target=[input_source;input_target];
% cue
subplot(3,2,1)
%G = digraph(channel_label(input_source(cue_links)), channel_label(input_target(cue_links)))
input_links=template_source_to_target(:,stable_links_indices_expert.Hit.cue);
G = digraph(channel_label(input_links(1,:)), channel_label(input_links(2,:)));

plot(G,'Layout','force')
title('Hit cue')

subplot(3,2,2)
input_links=template_source_to_target(:,cue_links);
G = digraph(channel_label(input_links(1,:)), channel_label(input_links(2,:)));

plot(G,'Layout','force')
title('Hit Shared naive to expert cue')

subplot(3,2,3)
% cue
%G = digraph(channel_label(input_source(tex_links)), channel_label(input_target(tex_links)))
input_links=template_source_to_target(:,stable_links_indices_expert.Hit.tex);
G = digraph(channel_label(input_links(1,:)), channel_label(input_links(2,:)));
plot(G,'Layout','force')

title('Hit tex')

subplot(3,2,4)
input_links=template_source_to_target(:,tex_links);
G = digraph(channel_label(input_links(1,:)), channel_label(input_links(2,:)));

plot(G,'Layout','force')
title('Hit Shared naive to expert tex')

subplot(3,2,5)
% cue
%G = digraph(channel_label(input_source(rew_links)), channel_label(input_target(rew_links)))
input_links=template_source_to_target(:,stable_links_indices_expert.Hit.rew);
G = digraph(channel_label(input_links(1,:)), channel_label(input_links(2,:)));
plot(G,'Layout','force')
title('Hit reward')

subplot(3,2,6)
input_links=template_source_to_target(:,rew_links);
G = digraph(channel_label(input_links(1,:)), channel_label(input_links(2,:)));

plot(G,'Layout','force')
title('Hit Shared naive to expert rew')

%% stable links indeces Hit

subplot(3,1,1)
y=[40 length(stable_links_indices.Hit.cue); 40 length(stable_links_indices_expert.Hit.cue); 40 length(cue_links)];
bar(y,'stacked')
title('cue related')
ylim([0 60])

subplot(3,1,2)
y=[40 length(stable_links_indices.Hit.tex); 40 length(stable_links_indices_expert.Hit.tex); 40 length(tex_links)];
bar(y,'stacked')
title('stimulus related')
ylim([0 60])

subplot(3,1,3)
y=[40 length(stable_links_indices.Hit.rew); 40 length(stable_links_indices_expert.Hit.rew); 40 length(rew_links)];
bar(y,'stacked')
title('action-reward related')
ylim([0 60])
%% --------------------------------
% --- calculate for 48 channels ---
% ---------------------------------

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

%% plot p-value of all links

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
    
        if count_source~=count_target

            for count_mice=1:mice
                if isempty(outdegree_vector.Hit.bivar{1,count_source}{1,count_mice})
                else
                    
                    input= mean(outdegree_vector.Hit.multivar{1,count_source}{1,count_mice}{1,count_time_interval},1,'omitnan');
                    % change the size of the empty input to concatenate
                    if isempty(input)
                        input=double.empty(0,48);
                    end
                    input_Hit=cat(1, input_Hit,input);
                    
                    input= mean(outdegree_vector.CR.multivar{1,count_source}{1,count_mice}{1,count_time_interval},1,'omitnan');
                    % change the size of the empty input to concatenate
                    if isempty(input)
                        input=double.empty(0,48);
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

%% --- plot top 95th percentile of the strongest links ---

ind=[];
%make node names nonunique
for i=1:length(channel_label)
    if sum(strcmp(channel_label{1,i},channel_label))>1
        % find indeces
        [r,c]=find(strcmp(channel_label{1,i},channel_label)==1);
            channel_label{1,i}=channel_label{1,i};
            channel_label{1,c(2)}=[channel_label{1,i} '2'];
    end
end

% network_time_interval
time_interval=25:27;
percent_to_threshold=99;

subplot(2,2,1)

% Hit networks
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=squeeze(mean(mean_link_frequency_Hit(time_interval,:,:),1,'omitnan'));
% find the 95th percentile
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);

% use percentile threshold from texture_time_interval 15:17;
%percentile_threshold=0.0994;

G = digraph(input_to_plot>percentile_threshold,channel_label)
centrality_values = centrality(G,'pagerank') % 'pagerank'

%test_links=G.Edges;

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

G = digraph(input_to_plot>percentile_threshold,channel_label)
centrality_values = centrality(G,'pagerank') % 'pagerank'

% test_links=G.Edges;
% template_links=top_channels.Hit;
% 
% count=1;
% ind_template_CR=[];
% ind_test_CR=[];
% for count_template=1:size(template_links,1)
%     for count_test=1:size(test_links,1)
%         
%         input_test= [test_links{count_test,:}];
%         input_template= [template_links{count_template,:}];
% 
%         if strcmp(input_test ,input_template)>0
%             ind_template_CR(count)=count_template;
%             ind_test_CR(count)=count_test;
%             count=count+1;
%         end
%     end
% end

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

%% calculate indeces for top and bottom links as compared to 12 ch

% test across all possible links in 48ch
% input_source=reshape(ones(1,48).*(1:48)',1,[])';
% input_target=reshape((1:48).*ones(48,1),1,[])';
% test_links=[channel_label(input_source); channel_label(input_target)]';
% network_time_interval
time_interval=15:17;
percent_to_threshold=10;

subplot(2,2,1)

% Hit networks
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=squeeze(mean(mean_link_frequency_Hit(time_interval,:,:),1,'omitnan'));
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);
G = digraph(input_to_plot>percentile_threshold,channel_label);
test_links=G.Edges;

% template from the top 25% of 12ch
template_links=top_channels.Hit;

% --- Hit ---

count=1;
top_template_Hit=[];
top_test_Hit=[];
for count_template=1:size(template_links,1)
    for count_test=1:size(test_links,1)
        
        %input_test= {test_links{count_test,1},test_links{count_test,2}};
        %input_template= template_links{count_template,1};

        input_test= [test_links{count_test,:}];
        input_template= [template_links{count_template,:}];

        if strcmp(input_test ,input_template)>0
            top_template_Hit(count)=count_template;
            top_test_Hit(count)=count_test;
            count=count+1;
        end
    end
end

% template from the top 25% of 12ch
template_links=bottom_channels.Hit;

count=1;
bottom_template_Hit=[];
bottom_test_Hit=[];

for count_template=1:size(template_links,1)
    for count_test=1:size(test_links,1)
        
%         input_test= {test_links{count_test,1},test_links{count_test,2}};
%         input_template= template_links{count_template,1};
        input_test= [test_links{count_test,:}];
        input_template= [template_links{count_template,:}];

        if strcmp(input_test ,input_template)>0
            bottom_template_Hit(count)=count_template;
            bottom_test_Hit(count)=count_test;
            count=count+1;
        end
    end
end

% Hit networks
dim1=size(mean_link_frequency_CR,1);
input_to_plot=squeeze(mean(mean_link_frequency_CR(time_interval,:,:),1,'omitnan'));
percentile_threshold=prctile(reshape(input_to_plot,1,[]),percent_to_threshold);
G = digraph(input_to_plot>percentile_threshold,channel_label);
test_links=G.Edges;

% template from the top 25% of 12ch
template_links=top_channels.Hit;

% --- CR ---

% template from the top 25% of 12ch
template_links=top_channels.CR;

count=1;
top_template_CR=[];
top_test_CR=[];
for count_template=1:size(template_links,1)
    for count_test=1:size(test_links,1)
        
%         input_test= {test_links{count_test,1},test_links{count_test,2}};
%         input_template= template_links{count_template,1};
        input_test= [test_links{count_test,:}];
        input_template= [template_links{count_template,:}];

        if strcmp(input_test ,input_template)>0
            top_template_CR(count)=count_template;
            top_test_CR(count)=count_test;
            count=count+1;
        end
    end
end

% template from the top 25% of 12ch
template_links=bottom_channels.CR;

count=1;
bottom_template_CR=[];
bottom_test_CR=[];

for count_template=1:size(template_links,1)
    for count_test=1:size(test_links,1)
        
%         input_test= {test_links{count_test,1},test_links{count_test,2}};
%         input_template= template_links{count_template,1};
        input_test= [test_links{count_test,:}];
        input_template= [template_links{count_template,:}];

        if strcmp(input_test ,input_template)>0
            bottom_template_CR(count)=count_template;
            bottom_test_CR(count)=count_test;
            count=count+1;
        end
    end
end

%%
cmax=0.1;
% sort by mean in a time iterval
time_interval=15:17;

figure
% Hit trials sorted with own index to plot
subplot(1,8,1)
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=reshape(mean_link_frequency_Hit,dim1,[]);

% sort by the time of the peak
% [~,ind]=max(input_to_plot,[],1)
% [~, index_to_plot_Hit] =sort(ind,'ascend');

% sort by the max of mean
[~, index_to_plot_Hit] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_Hit)')
shading interp;
% 'Colormap',redbluecmap
%colormap(redbluecmap)
% ax=gca;
% mymap=colormap(ax);
% save('mycolormap','mymap')

%colormap('jet')
c = 1-gray(20);
colormap(c)

xlim([1 50])
caxis([0 cmax])
title('mean Hit')

%colormapeditor

% ax = gca;
% mymap = colormap(ax);
% save('MyColormap3','mymap')

% --- do a plot for the Hit link strength sidebar
subplot(1,8,2)
input_linespace1=zeros(1,48*48);
input_linespace1(top_test_Hit)=1;

bar(input_linespace1','FaceColor',[1 0.5 0.5],'EdgeColor',[1 0 0],'LineWidth',0.5)
hold on

input_linespace2=zeros(1,48*48);
input_linespace2(bottom_test_Hit)=1;
bar(input_linespace2','FaceColor',[0.5 0.5 1],'EdgeColor',[0 0 1],'LineWidth',0.5)
camroll(-90)
hold off

% CR trials sorted with own index to plot
subplot(1,8,3)
dim1=size(mean_link_frequency_CR,1);
input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);
[~, ~] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');
imagesc(input_to_plot(:,index_to_plot_Hit)')
xlim([1 50])
caxis([0 cmax])

% CR trials sorted with own index to plot
subplot(1,8,4)
dim1=size(mean_link_frequency_CR,1);
input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);
[~, index_to_plot_CR] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');
imagesc(input_to_plot(:,index_to_plot_CR)')
xlim([1 50])
caxis([0 cmax])

% --- do a plot for the Hit link strength sidebar
subplot(1,8,5)
input_linespace1=zeros(1,48*48);
input_linespace1(top_test_CR)=1;

bar(input_linespace1','FaceColor',[1 0.5 0.5],'EdgeColor',[1 0 0],'LineWidth',0.5)
hold on

input_linespace2=zeros(1,48*48);
input_linespace2(bottom_test_CR)=1;
bar(input_linespace2','FaceColor',[0.5 0.5 1],'EdgeColor',[0 0 1],'LineWidth',0.5)
camroll(-90)
hold off

% labels by source area
% subplot(1,6,3)
% input_to_plot=reshape(ones(1,12).*(1:12)',1,[])';
% imagesc(input_to_plot(index_to_plot_Hit))
% xlim([0.5 1])

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
% subplot(1,6,4)
% input_to_plot=reshape((1:12).*ones(12,1),1,[])';
% imagesc(input_to_plot(index_to_plot_Hit))
% xlim([0.5 1])

% Hit trials sorted with own index to plot
subplot(1,6,5)
dim1=size(sem_link_frequency_Hit,1);
input_to_plot=reshape(sem_link_frequency_Hit,dim1,[]);
[~, index_to_plot_Hit] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_Hit)')
shading interp;
% 'Colormap',redbluecmap
%colormap(mymap)
xlim([1 50])
% caxis([0 cmax])

colormap(c)
title('sem Hit')
% ax = gca;
% mymap = colormap(ax);
% save('MyColormap3','mymap')

% CR trials sorted with own index to plot
subplot(1,6,6)
dim1=size(sem_link_frequency_CR,1);
input_to_plot=reshape(sem_link_frequency_CR,dim1,[]);
[~, index_to_plot_CR] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');
imagesc(input_to_plot(:,index_to_plot_CR)')
xlim([1 50])
caxis([0 cmax])
title('sem CR')


%% plot sequence like link activation and sorted by the mean in texture window

figure
cmax=0.1;
input_to_plot_bin=[];

% ---- Hit ----
ax(1) = subplot(1,5,1)
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=reshape(mean_link_frequency_Hit,dim1,[]);

% sort by the time of the peak in the binned
% bin with 5 samples
bin=5;
count_bin=1;
for i=1:5:45
    input_to_plot_bin(count_bin,:)= mean(input_to_plot(i:i+bin,:),1,'omitnan');
    count_bin=count_bin+1;
end
for i=1:size(input_to_plot,2)
    input_to_plot_smooth(:,i) = smooth(input_to_plot(:,i),'moving');
end
% use here input_to_plot_smooth instead of input_to_plot_bin
[~,ind]=max(input_to_plot_smooth,[],1)
[~, index_to_plot_Hit] =sort(ind,'ascend');

% sort by the max of mean
%[~, index_to_plot_Hit] =sort(mean(input_to_plot(14:17,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_Hit(48:end))')
shading interp;
% 'Colormap',redbluecmap
c = 1-gray(20);
colormap(ax(1), c)
xlim([1 50])
caxis([0 cmax])
title('mean Hit')

% Hit trials sorted with own index to plot
ax(2) = subplot(1,5,2)
dim1=size(mean_link_frequency_Hit,1);
input_to_plot=reshape(mean_link_frequency_Hit,dim1,[]);

% sort by the time of the peak
% [~,ind]=max(input_to_plot,[],1)
% [~, index_to_plot_Hit] =sort(ind,'ascend');

% sort by the max of mean
[~, index_to_plot_Hit] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');

cue_baseline= mean(input_to_plot(1:10,:),'omitnan');

imagesc((input_to_plot(:,index_to_plot_Hit) - cue_baseline(index_to_plot_Hit))')
shading interp;
colormap(ax(2), redbluecmap)
xlim([1 50])
caxis([-0.05 0.05])
title('mean Hit')

input_to_plot_bin=[];
% ---- CR ----
ax(3) = subplot(1,5,3)
dim1=size(mean_link_frequency_CR,1);
input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);

% sort by the time of the peak in the binned
% bin with 5 samples
bin=5;
count_bin=1;
for i=1:5:45
    input_to_plot_bin(count_bin,:)= mean(input_to_plot(i:i+bin,:),1,'omitnan');
    count_bin=count_bin+1;
end

for i=1:size(input_to_plot,2)
    input_to_plot_smooth(:,i) = smooth(input_to_plot(:,i),'moving');
end
% use here input_to_plot_smooth instead of input_to_plot_bin
[~,ind]=max(input_to_plot_smooth,[],1)
[~, index_to_plot_CR] =sort(ind,'ascend');

% sort by the max of mean
%[~, index_to_plot_Hit] =sort(mean(input_to_plot(14:17,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_Hit(:))')
shading interp;
% 'Colormap',redbluecmap
c = 1-gray(20);
colormap(ax(3), c)
xlim([1 50])
caxis([0 cmax])
title('mean CR sorted by Hit')

ax(3) = subplot(1,5,4)
dim1=size(mean_link_frequency_CR,1);
input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);

% sort by the time of the peak in the binned
% bin with 5 samples
bin=5;
count_bin=1;
for i=1:5:45
    input_to_plot_bin(count_bin,:)= mean(input_to_plot(i:i+bin,:),1,'omitnan');
    count_bin=count_bin+1;
end

for i=1:size(input_to_plot,2)
    input_to_plot_smooth(:,i) = smooth(input_to_plot(:,i),'moving');
end
% use here input_to_plot_smooth instead of input_to_plot_bin
[~,ind]=max(input_to_plot_smooth,[],1)
[~, index_to_plot_CR] =sort(ind,'ascend');

% sort by the max of mean
%[~, index_to_plot_Hit] =sort(mean(input_to_plot(14:17,:),1,'omitnan'),'descend');

imagesc(input_to_plot(:,index_to_plot_CR(48:end))')
shading interp;
% 'Colormap',redbluecmap
c = 1-gray(20);
colormap(ax(3), c)
xlim([1 50])
caxis([0 cmax])
title('mean CR sorted by CR')

% Hit trials sorted with own index to plot
ax(4) = subplot(1,5,5)
dim1=size(mean_link_frequency_CR,1);
input_to_plot=reshape(mean_link_frequency_CR,dim1,[]);
cue_baseline= mean(input_to_plot(1:10,:),'omitnan');
% sort by the time of the peak
% [~,ind]=max(input_to_plot,[],1)
% [~, index_to_plot_Hit] =sort(ind,'ascend');

% sort by the max of mean
[~, index_to_plot_CR] =sort(mean(input_to_plot(time_interval,:),1,'omitnan'),'descend');

imagesc((input_to_plot(:,index_to_plot_CR) - cue_baseline(index_to_plot_CR))')
shading interp;
colormap(ax(4), redbluecmap)
xlim([1 50])
caxis([-0.05 0.05])
title('mean CR')


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
p=plot(G,'Layout','force');
p.MarkerSize = 20*centrality_values+0.01;
title('strongest 14 Sources CR texture presentation')

subplot(2,2,4)

semilogy(cumsum(sorted_input_to_plot))
hold on
vline(top_sources)
hold off
ylim([0 25])