
mice=length(mouse_IDs);

% cannel labels for pooling
% get selected channels
label=[];
channel_label = {'GP', 'DG', 'VL','Cpu', 'CA1_Py', 'S1_bf','VM', 'LDVL', 'Rt','M1','CA1_Mol','dCA1'};

comb=length(channel_label);
        
pool_mice=[1 2 3 4 5 6 7 8 9 10 11 12 13];
for count_i_mice=1:length(pool_mice)

count_mice=pool_mice(count_i_mice);
session_names=sessions_IDs{1,count_mice};
sessions=length(session_names);

% define thresholds
network_p_threshold = 0.01;
network_performance_threshold = 0.7;
pagerank_constant = 0.85;
network_frequency_threshold = 0.5;
% define as 0 for naive animals i.e. pool sessions below the performance threshold
% define as 1 for expert animals i.e. pool sessions above the performance threshold
perf_def = 1;

% --- import channel labels ---
channel_labels=channel_labels_all{1,count_mice};
count_time_step=1;

% 3:3:16
for count_time_intervals =1:3:50-2
    
% --- Choose time intervals (in seconds) ---
% CUE=[1,1.5], TEX=[2,3.5], LIK=[3.5,6]
%network_time_interval = {[1,1.5],[3,3.5],[6,6.5]};
%network_time_interval = {[3,3.5]};

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
                
                %temp_data.parameteres=data_TE{1,count_mice}{1,count(ii)}.NIfTy.Hit.results.parameters;
                % select timesteps in the desired time interval
                %network_time_steps =
                %find(temp_data.parameteres.samples_timesteps >= min(network_time_interval{1,1}) & temp_data.parameteres.samples_timesteps < max(network_time_interval{1,1}))
                network_time_steps =(count_time_intervals:count_time_intervals+2);

                if length(network_time_steps) < 1
                    error('ERROR: no data points found in the specified time interval.')
                end
                
                % store the new adjacency matrix (assumption: all training sessions have the same number of channels)
                if strcmp(trials_type{1,1},'Hit')
                    
                    temp_data.bivar=data_TE{1,count_mice}{1,count(ii)}.bivar.Hit.results;
                    network_adj_matrices.GO.bivar(:,:,count(ii)) = network_bivar_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.multivar=data_TE{1,count_mice}{1,count(ii)}.multivar.Hit.results;
                    network_adj_matrices.GO.multivar(:,:,count(ii)) = network_bivar_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.NIfTy=data_TE{1,count_mice}{1,count(ii)}.NIfTy.Hit.results;
                    network_adj_matrices.GO.NIfTy(:,:,count(ii)) = network_threshold(temp_data.NIfTy.p_table,network_time_steps,network_p_threshold);
                    
                elseif strcmp(trials_type{1,1},'CR')
                    
                    temp_data.bivar=data_TE{1,count_mice}{1,count(ii)}.bivar.CR.results;
                    network_adj_matrices.NOGO.bivar(:,:,count(ii)) = network_bivar_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                    temp_data.multivar=data_TE{1,count_mice}{1,count(ii)}.multivar.CR.results;
                    network_adj_matrices.NOGO.multivar(:,:,count(ii)) = network_bivar_threshold(temp_data.multivar.p_table,network_time_steps,network_p_threshold);
                    
                    temp_data.NIfTy=data_TE{1,count_mice}{1,count(ii)}.NIfTy.CR.results;
                    network_adj_matrices.NOGO.NIfTy(:,:,count(ii)) = network_threshold(temp_data.NIfTy.p_table,network_time_steps,network_p_threshold);
                    
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
        outdegree.mean_GO.NIfTy(count_channel_labels,count_i_mice,count_time_step) =NaN;
        outdegree.var_GO.NIfTy(count_channel_labels,count_i_mice,count_time_step) =NaN;
        outdegree.mean_GO.bivar(count_channel_labels,count_i_mice,count_time_step)=NaN;
        outdegree.var_GO.bivar(count_channel_labels,count_i_mice,count_time_step)=NaN;
        outdegree.mean_GO.multivar(count_channel_labels,count_i_mice,count_time_step)=NaN;
        outdegree.var_GO.multivar(count_channel_labels,count_i_mice,count_time_step)=NaN;
        
        outdegree.mean_NOGO.NIfTy(count_channel_labels,count_i_mice,count_time_step)=NaN;
        outdegree.var_NOGO.NIfTy(count_channel_labels,count_i_mice,count_time_step)=NaN;
        outdegree.mean_NOGO.bivar(count_channel_labels,count_i_mice,count_time_step)=NaN;
        outdegree.var_NOGO.bivar(count_channel_labels,count_i_mice,count_time_step)=NaN;
        outdegree.mean_NOGO.multivar(count_channel_labels,count_i_mice,count_time_step)=NaN;
        outdegree.var_NOGO.multivar(count_channel_labels,count_i_mice,count_time_step)=NaN;
        
        clustering_coeff.mean_GO.NIfTy(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.var_GO.NIfTy(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.mean_GO.bivar(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.var_GO.bivar(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.mean_GO.multivar(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.var_GO.multivar(count_channel_labels,count_i_mice,count_time_step) =NaN;
        
        clustering_coeff.mean_NOGO.NIfTy(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.var_NOGO.NIfTy(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.mean_NOGO.bivar(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.var_NOGO.bivar(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.mean_NOGO.multivar(count_channel_labels,count_i_mice,count_time_step) =NaN;
        clustering_coeff.var_NOGO.multivar(count_channel_labels,count_i_mice,count_time_step) =NaN;
    else
    
    current_channel_labels=channel_labels_all{count_mice};
    
    % ------------------------
    % ---------- GO ----------
    % --------- bivar---------
    
    temp_input= network_adj_matrices.GO.bivar(:,:,:);
    scale_to_network=size(temp_input,1)*(size(temp_input,1)-1);
    
    network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
    outdegree.mean_GO.bivar(count_channel_labels,count_i_mice,count_time_step) = mean(network_outdegree,1)/scale_to_network;
    outdegree.var_GO.bivar(count_channel_labels,count_i_mice,count_time_step) = std(network_outdegree,0,1)/scale_to_network;
    
    network_cluster_GO = NaN([size(temp_input,1),length(high_perf_sessions_indices)]);
    
    clustered_nodes_GO=[];
    for count_sessions = 1:length(high_perf_sessions_indices)
        network_cluster_GO(:,count_sessions) = weighted_clust_coeff(temp_input(:,:,high_perf_sessions_indices(count_sessions)));
        %network_cluster_GO(:,count_sessions) = clustering_coef_wd(temp_input(:,:,high_perf_sessions_indices(count_sessions)));

        % find connected nodes from the source
        neigh=kneighbors(temp_input(:,:,high_perf_sessions_indices(count_sessions)),channel_sub_network,1);
        % find if neighbors are connected
        neigh_cat=[];
        if length(neigh)>1
            for jj=1:length(neigh)
                % neighbors of neighbors
                neigh_2=kneighbors(temp_input(:,:,high_perf_sessions_indices(count_sessions)),neigh(jj),1);
                % concatenate all neighbors of neighbors
                neigh_cat=[neigh_cat neigh_2];
            end
        end
        % intersect neighbors and all neigh_2
        connected_neigh=intersect(neigh,neigh_cat);
        if connected_neigh>0
            for jj=1:length(connected_neigh)
                clustered_nodes_GO=[clustered_nodes_GO current_channel_labels(connected_neigh(jj))];
            end
        end
        
    end
    
    clustering_coeff.mean_GO.bivar(count_channel_labels,count_i_mice,count_time_step) = mean(network_cluster_GO(channel_sub_network,:),2);
    clustering_coeff.var_GO.bivar(count_channel_labels,count_i_mice,count_time_step) = var(network_cluster_GO(channel_sub_network,:),0,2);
    
    % nodes conecting when count_channel_labels is connected to both
    connected_nodes_bivar_GO{1,count_channel_labels}{1,count_i_mice}=clustered_nodes_GO;
    
    % ------------------------
    % ---------- GO ----------
    % ------- multivar -------

    temp_input= network_adj_matrices.GO.multivar(:,:,:);
    network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
    outdegree.mean_GO.multivar(count_channel_labels,count_i_mice,count_time_step) = mean(network_outdegree,1)/scale_to_network;
    outdegree.var_GO.multivar(count_channel_labels,count_i_mice,count_time_step) = std(network_outdegree,0,1)/scale_to_network;
    
    network_cluster_GO = NaN([size(temp_input,1),length(high_perf_sessions_indices)]);

    clustered_nodes_GO=[];
    for count_sessions = 1:length(high_perf_sessions_indices)
        network_cluster_GO(:,count_sessions) = weighted_clust_coeff(temp_input(:,:,high_perf_sessions_indices(count_sessions)));
        %network_cluster_GO(:,count_sessions) = clustering_coef_wd(temp_input(:,:,high_perf_sessions_indices(count_sessions)));

        % find connected nodes from the source
        neigh=kneighbors(temp_input(:,:,high_perf_sessions_indices(count_sessions)),channel_sub_network,1);
        % find if neighbors are connected
        neigh_cat=[];
        if length(neigh)>1
            for jj=1:length(neigh)
                % neighbors of neighbors
                neigh_2=kneighbors(temp_input(:,:,high_perf_sessions_indices(count_sessions)),neigh(jj),1);
                % concatenate all neighbors of neighbors
                neigh_cat=[neigh_cat neigh_2];
            end
        end
        % intersect neighbors and all neigh_2
        connected_neigh=intersect(neigh,neigh_cat);
        if connected_neigh>0
            for jj=1:length(connected_neigh)
                clustered_nodes_GO=[clustered_nodes_GO current_channel_labels(connected_neigh(jj))];
            end
        end
        
    end % count sessions
    
    clustering_coeff.mean_GO.multivar(count_channel_labels,count_i_mice,count_time_step) = mean(network_cluster_GO(channel_sub_network,:),2);
    clustering_coeff.var_GO.multivar(count_channel_labels,count_i_mice,count_time_step) = var(network_cluster_GO(channel_sub_network,:),0,2);
        
    % nodes conecting when count_channel_labels is connected to both
    connected_nodes_multivar_GO{1,count_channel_labels}{1,count_i_mice}=clustered_nodes_GO;
 
    % --------------------------
    % ---------- NOGO ---------- 
    % ---------- bivar ---------
    
    temp_input= network_adj_matrices.NOGO.bivar(:,:,:);
    network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
    outdegree.mean_NOGO.bivar(count_channel_labels,count_i_mice,count_time_step) = mean(network_outdegree,1)/scale_to_network;
    outdegree.var_NOGO.bivar(count_channel_labels,count_i_mice,count_time_step) = std(network_outdegree,0,1)/scale_to_network;
    
    network_cluster_NOGO = NaN([size(temp_input,1),length(high_perf_sessions_indices)]);

    clustered_nodes_NOGO=[];
    for count_sessions = 1:length(high_perf_sessions_indices)
        network_cluster_NOGO(:,count_sessions) = weighted_clust_coeff(temp_input(:,:,high_perf_sessions_indices(count_sessions)));
        %network_cluster_NOGO(:,count_sessions) = clustering_coef_wd(temp_input(:,:,high_perf_sessions_indices(count_sessions)));
        
        % find connected nodes from the source
        neigh=kneighbors(temp_input(:,:,high_perf_sessions_indices(count_sessions)),channel_sub_network,1);
        % find if neighbors are connected
        neigh_cat=[];
        if length(neigh)>1
            for jj=1:length(neigh)
                % neighbors of neighbors
                neigh_2=kneighbors(temp_input(:,:,high_perf_sessions_indices(count_sessions)),neigh(jj),1);
                % concatenate all neighbors of neighbors
                neigh_cat=[neigh_cat neigh_2];
            end
        end
        % intersect neighbors and all neigh_2
        connected_neigh=intersect(neigh,neigh_cat);
        if connected_neigh>0
            for jj=1:length(connected_neigh)
                clustered_nodes_NOGO=[clustered_nodes_NOGO current_channel_labels(connected_neigh(jj))];
            end
        end
    end % count sessions
    
    % save clustering coeff and labels of connected nodes
    clustering_coeff.mean_NOGO.bivar(count_channel_labels,count_i_mice,count_time_step) = mean(network_cluster_NOGO(channel_sub_network,:),2);
    clustering_coeff.var_NOGO.bivar(count_channel_labels,count_i_mice,count_time_step) = var(network_cluster_NOGO(channel_sub_network,:),0,2);
    
    % nodes conecting when count_channel_labels is connected to both
    connected_nodes_bivar_NOGO{1,count_channel_labels}{1,count_i_mice}=clustered_nodes_NOGO;
    
    % --------------------------
    % ---------- NOGO ---------- 
    % -------- multivar --------
    temp_input= network_adj_matrices.NOGO.multivar(:,:,:);
    network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
    outdegree.mean_NOGO.multivar(count_channel_labels,count_i_mice,count_time_step) = mean(network_outdegree,1)/scale_to_network;
    outdegree.var_NOGO.multivar(count_channel_labels,count_i_mice,count_time_step) = std(network_outdegree,0,1)/scale_to_network;
        
    network_cluster_NOGO = NaN([size(temp_input,1),length(high_perf_sessions_indices)]);

    clustered_nodes_NOGO=[];
    for count_sessions = 1:length(high_perf_sessions_indices)
        network_cluster_NOGO(:,count_sessions) = weighted_clust_coeff(temp_input(:,:,high_perf_sessions_indices(count_sessions)));
        %network_cluster_NOGO(:,count_sessions) = clustering_coef_wd(temp_input(:,:,high_perf_sessions_indices(count_sessions)));
        
        % find connected nodes from the source
        neigh=kneighbors(temp_input(:,:,high_perf_sessions_indices(count_sessions)),channel_sub_network,1);
        % find if neighbors are connected
        neigh_cat=[];
        if length(neigh)>1
            for jj=1:length(neigh)
                % neighbors of neighbors
                neigh_2=kneighbors(temp_input(:,:,high_perf_sessions_indices(count_sessions)),neigh(jj),1);
                % concatenate all neighbors of neighbors
                neigh_cat=[neigh_cat neigh_2];
            end
        end
        % intersect neighbors and all neigh_2
        connected_neigh=intersect(neigh,neigh_cat);
        if connected_neigh>0
            for jj=1:length(connected_neigh)
                clustered_nodes_NOGO=[clustered_nodes_NOGO current_channel_labels(connected_neigh(jj))];
            end
        end
    end
    
    clustering_coeff.mean_NOGO.multivar(count_channel_labels,count_i_mice,count_time_step) = mean(network_cluster_NOGO(channel_sub_network,:),2);
    clustering_coeff.var_NOGO.multivar(count_channel_labels,count_i_mice,count_time_step) = var(network_cluster_NOGO(channel_sub_network,:),0,2);
    % nodes conecting when count_channel_labels is connected to both 
    connected_nodes_multivar_NOGO{1,count_channel_labels}{1,count_i_mice}=clustered_nodes_NOGO;

    end
    
end

count_time_step=count_time_step+1;

end % count_time_intervals   
end % count_mice

%% plot compare out degree for Hit vs CR

ch=1;
t=0.2:0.6:10-0.4;
%plot for Hit trials
subplot(4,1,1)
input_to_plot= squeeze(mean(outdegree.mean_GO.multivar(ch,:,:),2,'omitnan'));
input_to_sem=squeeze(std(outdegree.mean_GO.multivar(ch,:,:)/sqrt(count_mice)));
shadedErrorBar(t,input_to_plot,input_to_sem,'-b',0)
title(channel_label{1,ch})

% plot for CR trials
subplot(4,1,2)
input_to_plot= squeeze(mean(outdegree.mean_NOGO.multivar(ch,:,:),2,'omitnan'));
input_to_sem=squeeze(std(outdegree.mean_NOGO.multivar(ch,:,:)/sqrt(count_mice)));
shadedErrorBar(t,input_to_plot,input_to_sem,'-r',0)
title(channel_label{1,ch})


% plot for CR trials
subplot(4,1,3)
input_to_plot= squeeze(mean(outdegree.mean_GO.multivar(ch,:,:) - outdegree.mean_NOGO.multivar(ch,:,:),2,'omitnan'));
input_to_sem=squeeze(std(outdegree.mean_GO.multivar(ch,:,:) - outdegree.mean_NOGO.multivar(ch,:,:),[],2,'omitnan')/sqrt(count_mice))
shadedErrorBar(t,input_to_plot,input_to_sem,'-b',0)
title(channel_label{1,ch})

subplot(4,1,4)
for i=1:length(t)
    [p,h]=signrank(outdegree.mean_GO.multivar(ch,:,i),outdegree.mean_NOGO.multivar(ch,:,i));
    p_val(i)=p;
end
p_val(p_val>0.01)=1;
p_val(p_val<=0.01)=0;
map=[0, 0, 0
    1, 1, 1];
imagesc(p_val)
colormap(map)

%% outdegree map of all regions

input_to_plot= squeeze(mean(outdegree.mean_GO.multivar(:,:,:) - outdegree.mean_NOGO.multivar(:,:,:),2,'omitnan'));
imagesc(input_to_plot)
colormap(mymap)

%% plot compare bivar clustering coeff for Hit vs CR
%mouse_pool=[3 5 7 9 10 11 12 13];
mouse_pool=(1:13);
count_mice=length(mouse_pool);
ch=1;
t=0.2:0.6:10-0.4;
%plot for Hit trials
subplot(4,1,1)
input_to_plot= squeeze(mean(clustering_coeff.mean_GO.bivar(ch,mouse_pool,:),2,'omitnan'));
input_to_sem=squeeze(std(clustering_coeff.mean_GO.bivar(ch,mouse_pool,:)/sqrt(count_mice)));
shadedErrorBar(t,input_to_plot,input_to_sem,'-b',0)
title(channel_label{1,ch})

% plot for CR trials
subplot(4,1,2)
input_to_plot= squeeze(mean(clustering_coeff.mean_NOGO.bivar(ch,mouse_pool,:),2,'omitnan'));
input_to_sem=squeeze(std(clustering_coeff.mean_NOGO.bivar(ch,mouse_pool,:)/sqrt(count_mice)));
shadedErrorBar(t,input_to_plot,input_to_sem,'-r',0)
title(channel_label{1,ch})


% plot for CR trials
subplot(4,1,3)
input_to_plot= squeeze(mean(clustering_coeff.mean_GO.bivar(ch,mouse_pool,:) - clustering_coeff.mean_NOGO.bivar(ch,mouse_pool,:),2,'omitnan'));
input_to_sem=squeeze(std(clustering_coeff.mean_GO.bivar(ch,mouse_pool,:) - clustering_coeff.mean_NOGO.bivar(ch,mouse_pool,:),[],2,'omitnan')/sqrt(count_mice))
shadedErrorBar(t,input_to_plot,input_to_sem,'-k',0)
title(channel_label{1,ch})
ylim([-0.1 0.1])

subplot(4,1,4)
for i=1:length(t)
    [p,h]=signrank(clustering_coeff.mean_GO.bivar(ch,mouse_pool,i),clustering_coeff.mean_NOGO.bivar(ch,mouse_pool,i));
    p_val(i)=p;
end
p_val(p_val>0.01)=1;
p_val(p_val<=0.01)=0;
map=[0, 0, 0
    1, 1, 1];
imagesc(p_val)
colormap(map)

%% 

subplot(2,1,1)
input_to_plot= squeeze(mean(clustering_coeff.mean_GO.bivar(:,mouse_pool,:) - clustering_coeff.mean_NOGO.bivar(:,mouse_pool,:),2,'omitnan'));
imagesc(input_to_plot)
colormap(mymap)
caxis([-0.1 0.15])

% calculate p-values for all channels
for ch=1:12
    for i=1:length(t)
        [p,h]=signrank(clustering_coeff.mean_GO.bivar(ch,mouse_pool,i),clustering_coeff.mean_NOGO.bivar(ch,mouse_pool,i));
        p_val(ch,i)=p;
    end
end
p_val(p_val>0.05/12)=1;
p_val(p_val<=0.05/12)=0;
subplot(2,1,2)
imagesc(input_to_plot)
map=[0, 0, 0
    1, 1, 1];
imagesc(p_val)
colormap(map)

%% plot compare multivar clustering coeff for Hit vs CR
%mouse_pool=[3 5 7 9 10 11 12 13];
mouse_pool=(1:13);
count_mice=length(mouse_pool);
ch=5;
t=0.2:0.6:10-0.4;
%plot for Hit trials
subplot(4,1,1)
input_to_plot= squeeze(mean(clustering_coeff.mean_GO.multivar(ch,mouse_pool,:),2,'omitnan'));
input_to_sem=squeeze(std(clustering_coeff.mean_GO.multivar(ch,mouse_pool,:)/sqrt(count_mice)));
shadedErrorBar(t,input_to_plot,input_to_sem,'-b',0)
title(channel_label{1,ch})

% plot for CR trials
subplot(4,1,2)
input_to_plot= squeeze(mean(clustering_coeff.mean_NOGO.multivar(ch,mouse_pool,:),2,'omitnan'));
input_to_sem=squeeze(std(clustering_coeff.mean_NOGO.multivar(ch,mouse_pool,:)/sqrt(count_mice)));
shadedErrorBar(t,input_to_plot,input_to_sem,'-r',0)
title(channel_label{1,ch})


% plot for CR trials
subplot(4,1,3)
input_to_plot= squeeze(mean(clustering_coeff.mean_GO.multivar(ch,mouse_pool,:) - clustering_coeff.mean_NOGO.multivar(ch,mouse_pool,:),2,'omitnan'));
input_to_sem=squeeze(std(clustering_coeff.mean_GO.multivar(ch,mouse_pool,:) - clustering_coeff.mean_NOGO.multivar(ch,mouse_pool,:),[],2,'omitnan')/sqrt(count_mice))
shadedErrorBar(t,input_to_plot,input_to_sem,'-b',0)
title(channel_label{1,ch})
%ylim([-0.05 0.1])

subplot(4,1,4)
for i=1:length(t)
    [p,h]=signrank(clustering_coeff.mean_GO.multivar(ch,mouse_pool,i),clustering_coeff.mean_NOGO.multivar(ch,mouse_pool,i));
     %[h,p]=ttest(clustering_coeff.mean_GO.multivar(ch,mouse_pool,i),clustering_coeff.mean_NOGO.multivar(ch,mouse_pool,i));

    p_val(i)=p;
end
p_val(p_val>0.01)=1;
p_val(p_val<=0.01)=0;
map=[0, 0, 0
    1, 1, 1];
imagesc(p_val)
colormap(map)

%% multivariate TE map of all regions
subplot(2,1,1)
input_to_plot= squeeze(mean(clustering_coeff.mean_GO.multivar(:,mouse_pool,:) - clustering_coeff.mean_NOGO.multivar(:,mouse_pool,:),2,'omitnan'));
imagesc(input_to_plot)
colormap(mymap)

% calculate p-values for all channels
for ch=1:12
    for i=1:length(t)
        [p,h]=signrank(clustering_coeff.mean_GO.multivar(ch,mouse_pool,i),clustering_coeff.mean_NOGO.multivar(ch,mouse_pool,i));
        p_val(ch,i)=p;
    end
end
p_val(p_val>0.05/12)=1;
p_val(p_val<=0.05/12)=0;
subplot(2,1,2)
imagesc(input_to_plot)
map=[0, 0, 0
    1, 1, 1];
imagesc(p_val)
colormap(map)

%% plot diff of clustering across all regions
t=0.2:0.2:10-0.4;

for ch=1:12

subplot(4,3,ch)
input_to_plot= squeeze(mean(clustering_coeff.mean_GO.multivar(ch,:,:) - clustering_coeff.mean_NOGO.multivar(ch,:,:),2,'omitnan'));
input_to_sem=squeeze(std(clustering_coeff.mean_GO.multivar(ch,:,:) - clustering_coeff.mean_NOGO.multivar(ch,:,:),[],2,'omitnan')/sqrt(count_mice));
shadedErrorBar(t,input_to_plot,input_to_sem,'-b',1)
title(channel_label{1,ch})
ylim([-0.05 0.1])
hold on
hline(0.05)
hold off
end
%% --- boxplot plus lines ---
number_of_mice= length(pool_mice);
% input_field_1=clustering_coeff.mean_GO.NIfTy;
% input_field_2=clustering_coeff.mean_NOGO.NIfTy;

input_field_1=reshape(clustering_coeff.mean_GO.multivar,12,[]);
input_field_2=reshape(clustering_coeff.mean_NOGO.multivar,12,[]);

% input_field_1=outdegree.mean_GO.multivar;
% input_field_2=outdegree.mean_NOGO.multivar;

% [p,tbl,stats]=anova1([input_field_1(1,:); input_field_2(1,:)]');
% multcompare(stats)
for i=1:12
    %[h,p] =ttest2(input_field_1(i,:), input_field_2(i,:));
    [p,h] =signrank(input_field_1(i,:), input_field_2(i,:));
    p_val(i)=p;
end

input_1= std(input_field_1,0,2,'omitnan')./sqrt(number_of_mice);
input_2= mean(input_field_1,2,'omitnan');

input_3= std(input_field_2,0,2,'omitnan')./sqrt(number_of_mice);
input_4= mean(input_field_2,2,'omitnan');

% [p,tbl,stats]=anova1([input_field_1 input_field_2]');
% multcompare(stats)

h1=barwitherr([input_1 input_3],[input_2 input_4]);
set(gca,'XTickLabel',channel_label)
set(h1(1),'FaceColor',[0.0 0.0 0.6]);
%ylim([0 1]);
hold on

% for i=1:12
%
% plot([i-0.1:0.2:0.1+i],[input_field_1(i,:); input_field_2(i,:)],'-o',...
%     'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
%     'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);
%
% end


hold off
title(num2str(p));

%% re-name input fields

% input_field_1_e=input_field_1;
% input_field_2_e=input_field_2;

input_field_1_n=input_field_1;
input_field_2_n=input_field_2;
%% plot for one channel
ch=5;
figure
h1=barwitherr([input_1(ch,1) input_3(ch,1)],[input_2(ch,1) input_4(ch,1)]);
set(gca,'XTickLabel',{'Hit','CR'})
set(h1(1),'FaceColor',[0.0 0.0 0.6]);
%ylim([0 1]);
hold on

plot(1:1:2,[input_field_1(ch,:); input_field_2(ch,:)],'-o',...
    'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
    'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);
ylim([0 0.3])
title(channel_label(ch))
%% compare multivar on naive vs expert clstering coeff

% Hit
input_1_e= std(input_field_1_e,0,2,'omitnan')./sqrt(number_of_mice);
input_2_e= mean(input_field_1_e,2,'omitnan');
% CR
input_3_e= std(input_field_2_e,0,2,'omitnan')./sqrt(number_of_mice);
input_4_e= mean(input_field_2_e,2,'omitnan');

% Hit
input_1_n= std(input_field_1_n,0,2,'omitnan')./sqrt(number_of_mice);
input_2_n= mean(input_field_1_n,2,'omitnan');
% CR
input_3_n= std(input_field_2_n,0,2,'omitnan')./sqrt(number_of_mice);
input_4_n= mean(input_field_2_n,2,'omitnan');

% test Hit trials
for i=1:12
    [h,p] =ttest2(input_field_1_e(i,:), input_field_1_n(i,:));
    p_val_Hit(i)=p<0.05
end

for i=1:12
    [h,p] =ttest2(input_field_2_e(i,:), input_field_2_n(i,:));
    p_val_CR(i)=p<0.05
end
%channel_label{1,p_val}

% plot channels significantly changing clustering coeff for Hit trials
%figure
subplot(2,1,1)
h1=barwitherr([input_1_n(p_val_Hit) input_1_e(p_val_Hit)],[input_2_n(p_val_Hit) input_2_e(p_val_Hit)]);
set(gca,'XTickLabel',channel_label(p_val_Hit))
set(h1(1),'FaceColor',[0.0 0.0 0.6]);

hold on
count_plot=1;
for i=1:12
 if p_val_Hit(i)>0
    plot([count_plot-0.1:0.2:0.1+count_plot],[input_field_1_n(i,:); input_field_1_e(i,:)],'-o',...
    'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
    'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);
    count_plot=count_plot+1;
 end
end
xlim([0 8]);
ylim([0 0.3]);
title('Hit');
hold off

%figure
subplot(2,1,2)
h2=barwitherr([input_3_n(p_val_CR) input_3_e(p_val_CR)],[input_4_n(p_val_CR) input_4_e(p_val_CR)]);
set(gca,'XTickLabel',channel_label(p_val_CR))
set(h2(1),'FaceColor',[0.0 0.0 0.6]);

hold on
count_plot=1;
for i=1:12
 if p_val_CR(i)>0
    plot([count_plot-0.1:0.2:0.1+count_plot],[input_field_2_n(i,:); input_field_2_e(i,:)],'-o',...
    'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
    'MarkerEdgeColor','r','MarkerFaceColor',[1,1,1]);
    count_plot=count_plot+1;
 end
end
xlim([0 8]);
ylim([0 0.4]);
title('CR');
hold off

%% test
for i=1:12
    [h,p] =ttest2(input_field_1_e(i,:), input_field_1_n(i,:));
    p_val(i)=p<0.05
end
channel_label{1,p_val}