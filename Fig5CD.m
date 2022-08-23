%----- permute network? ----- 
permute_network=1;

% ----- define thresholds ----- 
network_p_threshold = 0.01;
network_performance_threshold = 0.7;
pagerank_constant = 0.85;
network_frequency_threshold = 0.5;

% ----- define as 0 for naive animals i.e. pool sessions below the performance threshold ----- 
% ----- define as 1 for expert animals i.e. pool sessions above the
% performance threshold ----- 
perf_def = 1;

% --- Choose time intervals (in seconds) ---
% CUE=[1,1.5], TEX=[2,3.5], LIK=[3.5,6]
%network_time_interval = {[1,1.5],[3,3.5],[6,6.5]};
network_time_interval = {[3,3.5]};

mice=length(mouse_IDs);

% cannel labels for pooling
% get selected channels
label=[];
channel_label = {'GP', 'DG', 'VL','Cpu', 'CA1_Py', 'S1_bf','VM', 'LDVL', 'Rt','M1','CA1_Mol','dCA1'};

comb=length(channel_label);
pool_mice=[1 2 3 4 5 6 7 8 9 10 11 12 13];

% --- add a loop for counting shuffles ---
number_of_shuffles= 30;

for count_shuffles = 1:number_of_shuffles
    
for count_i_mice=1:length(pool_mice)

count_mice=pool_mice(count_i_mice);
session_names=sessions_IDs{1,count_mice};
sessions=length(session_names);

% --- import channel labels ---
channel_labels=channel_labels_all{1,count_mice};

% --- choose NIFTy, bivariate and multivariate measures ---
% --- Generate adjacency matrices for GO and NOGO trials ---
    
clear network_adj_matrices

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
                network_time_steps = find(temp_data.parameteres.samples_timesteps >= min(network_time_interval{1,1}) & temp_data.parameteres.samples_timesteps < max(network_time_interval{1,1}))-1

                if length(network_time_steps) < 1
                    error('ERROR: no data points found in the specified time interval.')
                end
                
                temp_adj_matrices_GO =[];
                temp_adj_matrices_NOGO =[];
                for trials_type = {'Hit','CR'}
                
                % ---do a loop over network_time_steps---
                for time_frames=1:length(network_time_steps)
                    
                % store the new adjacency matrix (assumption: all training sessions have the same number of channels)
                if strcmp(trials_type{1,1},'Hit')

        %                     temp_data.NIfTy=data_TE{1,count_mice}{1,count(ii)}.NIfTy.Hit.results;
        %                     network_adj_matrices.GO.NIfTy(:,:,count(ii)) = network_threshold(temp_data.NIfTy.p_table,network_time_steps,network_p_threshold);
        %                     temp_data.bivar=data_TE{1,count_mice}{1,count(ii)}.bivar.Hit.results;
        %                     network_adj_matrices.GO.bivar(:,:,count(ii)) = network_bivar_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                            temp_data.multivar=data_TE{1,count_mice}{1,count(ii)}.multivar.Hit.results;
                            temp_input= network_bivar_threshold(temp_data.multivar.p_table,network_time_steps(time_frames),network_p_threshold);
                            adj_dim=size(temp_input,1);
                            
                            if permute_network==1
                                % randomize netwrok while preserving the degree distribution
                                input_to_permute=reshape(temp_input,1,[]);
                                indeces_to_permute=randperm(length(input_to_permute));
                                adj_permute=reshape(input_to_permute(indeces_to_permute), adj_dim, adj_dim);
                                
                                temp_adj_matrices_GO = cat(3, temp_adj_matrices_GO, adj_permute);
                                
                            else
                                temp_adj_matrices_GO = cat(3, temp_adj_matrices_GO, temp_input);
                            end
                    disp(['session ' num2str(count_sessions) ' Hit time frame ' num2str(time_frames)])
                    
                    elseif strcmp(trials_type{1,1},'CR')
                        
        %                     temp_data.NIfTy=data_TE{1,count_mice}{1,count(ii)}.NIfTy.CR.results;
        %                     network_adj_matrices.NOGO.NIfTy(:,:,count(ii)) = network_threshold(temp_data.NIfTy.p_table,network_time_steps,network_p_threshold);
        %                     temp_data.bivar=data_TE{1,count_mice}{1,count(ii)}.bivar.CR.results;
        %                     network_adj_matrices.NOGO.bivar(:,:,count(ii)) = network_bivar_threshold(temp_data.bivar.p_table,network_time_steps,network_p_threshold);

                        temp_data.multivar=data_TE{1,count_mice}{1,count(ii)}.multivar.CR.results;
                        temp_input= network_bivar_threshold(temp_data.multivar.p_table,network_time_steps(time_frames),network_p_threshold);

                            if permute_network==1
                                % randomize netwrok while preserving the degree distribution
                                input_to_permute=reshape(temp_input,1,[]);
                                indeces_to_permute=randperm(length(input_to_permute));
                                adj_permute=reshape(input_to_permute(indeces_to_permute), adj_dim, adj_dim);
                                
                                temp_adj_matrices_NOGO = cat(3, temp_adj_matrices_NOGO, adj_permute);
                            else
                                temp_adj_matrices_NOGO = cat(3, temp_adj_matrices_NOGO, temp_input);
                            end

                    disp(['session ' num2str(count_sessions) ' Hit time frame ' num2str(time_frames)])
                    end % if strcmp for Hit vs CR
                    
                end % time frames
                end % count trial type

                network_adj_matrices.GO.multivar(:,:,count(ii)) = mean(temp_adj_matrices_GO,3,'omitnan');
                network_adj_matrices.NOGO.multivar(:,:,count(ii)) = mean(temp_adj_matrices_NOGO,3,'omitnan');
                ii=ii+1;

                end % if both NifTy and bivar TE exist
                
                end
    end % count sessions
   
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
        outdegree.mean_GO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        outdegree.var_GO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        outdegree.mean_GO.bivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_GO.bivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.mean_GO.multivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_GO.multivar(count_channel_labels,count_i_mice)=NaN;
        
        outdegree.mean_NOGO.NIfTy(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_NOGO.NIfTy(count_channel_labels,count_i_mice)=NaN;
        outdegree.mean_NOGO.bivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_NOGO.bivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.mean_NOGO.multivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_NOGO.multivar(count_channel_labels,count_i_mice)=NaN;
        
        clustering_coeff.mean_GO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_GO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.mean_GO.bivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_GO.bivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.mean_GO.multivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_GO.multivar(count_channel_labels,count_i_mice) =NaN;
        
        clustering_coeff.mean_NOGO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_NOGO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.mean_NOGO.bivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_NOGO.bivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.mean_NOGO.multivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_NOGO.multivar(count_channel_labels,count_i_mice) =NaN;
    else
    
    current_channel_labels=channel_labels_all{count_mice};
    
    % ------------------------
    % ---------- GO ----------
    % ------- multivar -------
    temp_input=[];
    disp(['size of the network ' num2str(size(network_adj_matrices.GO.multivar))])

    temp_input= network_adj_matrices.GO.multivar(:,:,:);
%     if permute_network==1
%         % permute each session separately and concatenate 
%         for count_sessions=1:size(network_adj_matrices.GO.multivar,3)
%             temp_input= cat(3,temp_input,...
%                 randmio_dir(network_adj_matrices.GO.multivar(:,:,count_sessions),20));
%         end
%     else
%         temp_input= network_adj_matrices.GO.multivar(:,:,:);
%     end
    
%     network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
%     outdegree.mean_GO.multivar(count_channel_labels,count_i_mice) = mean(network_outdegree,1)/scale_to_network;
%     outdegree.var_GO.multivar(count_channel_labels,count_i_mice) = std(network_outdegree,0,1)/scale_to_network;
    
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
    
    clustering_coeff.mean_GO.multivar(count_channel_labels,count_i_mice) = mean(network_cluster_GO(channel_sub_network,:),2);
    clustering_coeff.var_GO.multivar(count_channel_labels,count_i_mice) = var(network_cluster_GO(channel_sub_network,:),0,2);
        
    % nodes conecting when count_channel_labels is connected to both
    connected_nodes_multivar_shuffle_GO{1,count_channel_labels}{1,count_i_mice}{1,count_shuffles}=clustered_nodes_GO;
 
    
    % --------------------------
    % ---------- NOGO ---------- 
    % -------- multivar --------
    
    temp_input=[];
    disp(['size of the network ' num2str(size(network_adj_matrices.GO.multivar))])
    temp_input= network_adj_matrices.NOGO.multivar(:,:,:);
%     if permute_network==1
%         % permute each session separately and concatenate 
%         for count_sessions=1:size(network_adj_matrices.NOGO.multivar,3)
%             temp_input= cat(3,temp_input,...
%                 randmio_dir(network_adj_matrices.NOGO.multivar(:,:,count_sessions),20));
%         end
%     else
%         temp_input= network_adj_matrices.NOGO.multivar(:,:,:);
%     end
    
%     network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
%     outdegree.mean_NOGO.multivar(count_channel_labels,count_i_mice) = mean(network_outdegree,1)/scale_to_network;
%     outdegree.var_NOGO.multivar(count_channel_labels,count_i_mice) = std(network_outdegree,0,1)/scale_to_network;
        
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
    
    clustering_coeff.mean_NOGO.multivar(count_channel_labels,count_i_mice) = mean(network_cluster_NOGO(channel_sub_network,:),2);
    clustering_coeff.var_NOGO.multivar(count_channel_labels,count_i_mice) = var(network_cluster_NOGO(channel_sub_network,:),0,2);
    % nodes conecting when count_channel_labels is connected to both 
    connected_nodes_multivar_shuffle_NOGO{1,count_channel_labels}{1,count_i_mice}{1,count_shuffles}=clustered_nodes_NOGO;

    end
    
end
    
end % count_mice

end % count_shuffles

%% calculate connected nodes without permutation

connected_nodes_multivar_GO =cell(1,comb);
connected_nodes_multivar_NOGO =cell(1,comb);

% ----- permute network? --- 
permute_network=0;

mice=length(mouse_IDs);

% cannel labels for pooling
% get selected channels
label=[];
channel_label = {'GP', 'DG', 'VL','Cpu', 'CA1_Py', 'S1_bf','VM', 'LDVL', 'Rt','M1','CA1_Mol','dCA1'};

comb=length(channel_label);
pool_mice=[1 2 3 4 5 6 7 8 9 10 11 12 13];

for count_i_mice=1:length(pool_mice)

count_mice=pool_mice(count_i_mice)
session_names=sessions_IDs{1,count_mice};
sessions=length(session_names);

% define thresholds
network_p_threshold = 0.01;
network_performance_threshold = 0.7;
pagerank_constant = 0.85;
network_frequency_threshold = 0.5;
% define as 0 for naive animals i.e. pool sessions below the performance threshold
% define as 1 for expert animals i.e. pool sessions above the performance threshold

%---Performance to pool is defined in the shuffled clustering above---
%perf_def = 1;

% --- import channel labels ---
channel_labels=channel_labels_all{1,count_mice};

% --- Choose time intervals (in seconds) ---
%--- Network_time_interval is defined in the shuffled clustering above---
% network_time_interval = {[3 3.5]};

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
                network_time_steps = find(temp_data.parameteres.samples_timesteps >= min(network_time_interval{1,1}) & temp_data.parameteres.samples_timesteps < max(network_time_interval{1,1}))-1;

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
        outdegree.mean_GO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        outdegree.var_GO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        outdegree.mean_GO.bivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_GO.bivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.mean_GO.multivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_GO.multivar(count_channel_labels,count_i_mice)=NaN;
        
        outdegree.mean_NOGO.NIfTy(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_NOGO.NIfTy(count_channel_labels,count_i_mice)=NaN;
        outdegree.mean_NOGO.bivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_NOGO.bivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.mean_NOGO.multivar(count_channel_labels,count_i_mice)=NaN;
        outdegree.var_NOGO.multivar(count_channel_labels,count_i_mice)=NaN;
        
        clustering_coeff.mean_GO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_GO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.mean_GO.bivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_GO.bivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.mean_GO.multivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_GO.multivar(count_channel_labels,count_i_mice) =NaN;
        
        clustering_coeff.mean_NOGO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_NOGO.NIfTy(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.mean_NOGO.bivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_NOGO.bivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.mean_NOGO.multivar(count_channel_labels,count_i_mice) =NaN;
        clustering_coeff.var_NOGO.multivar(count_channel_labels,count_i_mice) =NaN;
    else
    
    current_channel_labels=channel_labels_all{count_mice};
    
    % ------------------------
    % ---------- GO ----------
    % --------- bivar---------
    temp_input=[];
    %temp_input= network_adj_matrices.GO.bivar(:,:,:);
    if permute_network==1
        % permute each session separately and concatenate 
        for count_sessions=1:size(network_adj_matrices.GO.bivar,3)
            temp_input= cat(3,temp_input,...
                randmio_dir(network_adj_matrices.GO.bivar(:,:,count_sessions),1));
        end
    else
        temp_input= network_adj_matrices.GO.bivar(:,:,:);
    end
    
    scale_to_network=size(temp_input,1)*(size(temp_input,1)-1);
    
    network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
    outdegree.mean_GO.bivar(count_channel_labels,count_i_mice) = mean(network_outdegree,1)/scale_to_network;
    outdegree.var_GO.bivar(count_channel_labels,count_i_mice) = std(network_outdegree,0,1)/scale_to_network;
    
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
    
    clustering_coeff.mean_GO.bivar(count_channel_labels,count_i_mice) = mean(network_cluster_GO(channel_sub_network,:),2);
    clustering_coeff.var_GO.bivar(count_channel_labels,count_i_mice) = var(network_cluster_GO(channel_sub_network,:),0,2);
    
    % nodes conecting when count_channel_labels is connected to both
    connected_nodes_bivar_GO{1,count_channel_labels}{1,count_i_mice}=clustered_nodes_GO;
    connected_nodes_sessions{1,count_channel_labels}{1,count_i_mice}=length(high_perf_sessions_indices);
    
    % ------------------------
    % ---------- GO ----------
    % ------- multivar -------

    %temp_input= network_adj_matrices.GO.multivar(:,:,:);
    if permute_network==1
        % permute each session separately and concatenate 
        for count_sessions=1:size(network_adj_matrices.GO.multivar,3)
            temp_input= cat(3,temp_input,...
                randmio_dir(network_adj_matrices.GO.multivar(:,:,count_sessions),1));
        end
    else
        temp_input= network_adj_matrices.GO.multivar(:,:,:);
    end
    
    network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
    outdegree.mean_GO.multivar(count_channel_labels,count_i_mice) = mean(network_outdegree,1)/scale_to_network;
    outdegree.var_GO.multivar(count_channel_labels,count_i_mice) = std(network_outdegree,0,1)/scale_to_network;
    
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
    
    clustering_coeff.mean_GO.multivar(count_channel_labels,count_i_mice) = mean(network_cluster_GO(channel_sub_network,:),2);
    clustering_coeff.var_GO.multivar(count_channel_labels,count_i_mice) = var(network_cluster_GO(channel_sub_network,:),0,2);
        
    % nodes conecting when count_channel_labels is connected to both
    connected_nodes_multivar_GO{1,count_channel_labels}{1,count_i_mice}=clustered_nodes_GO;
 
    % --------------------------
    % ---------- NOGO ---------- 
    % ---------- bivar ---------
    
    %temp_input= network_adj_matrices.NOGO.bivar(:,:,:);
    if permute_network==1
        % permute each session separately and concatenate 
        for count_sessions=1:size(network_adj_matrices.NOGO.bivar,3)
            temp_input= cat(3,temp_input,...
                randmio_dir(network_adj_matrices.NOGO.bivar(:,:,count_sessions),1));
        end
    else
        temp_input= network_adj_matrices.NOGO.bivar(:,:,:);
    end
    
    network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
    outdegree.mean_NOGO.bivar(count_channel_labels,count_i_mice) = mean(network_outdegree,1)/scale_to_network;
    outdegree.var_NOGO.bivar(count_channel_labels,count_i_mice) = std(network_outdegree,0,1)/scale_to_network;
    
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
    clustering_coeff.mean_NOGO.bivar(count_channel_labels,count_i_mice) = mean(network_cluster_NOGO(channel_sub_network,:),2);
    clustering_coeff.var_NOGO.bivar(count_channel_labels,count_i_mice) = var(network_cluster_NOGO(channel_sub_network,:),0,2);
    
    % nodes conecting when count_channel_labels is connected to both
    connected_nodes_bivar_NOGO{1,count_channel_labels}{1,count_i_mice}=clustered_nodes_NOGO;
    
    % --------------------------
    % ---------- NOGO ---------- 
    % -------- multivar --------
    %temp_input= network_adj_matrices.NOGO.multivar(:,:,:);
    if permute_network==1
        % permute each session separately and concatenate 
        for count_sessions=1:size(network_adj_matrices.NOGO.multivar,3)
            temp_input= cat(3,temp_input,...
                randmio_dir(network_adj_matrices.NOGO.multivar(:,:,count_sessions),1));
        end
    else
        temp_input= network_adj_matrices.NOGO.multivar(:,:,:);
    end
    
    network_outdegree = squeeze(sum(temp_input(channel_sub_network,:,high_perf_sessions)>0,2));
    outdegree.mean_NOGO.multivar(count_channel_labels,count_i_mice) = mean(network_outdegree,1)/scale_to_network;
    outdegree.var_NOGO.multivar(count_channel_labels,count_i_mice) = std(network_outdegree,0,1)/scale_to_network;
        
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
    
    clustering_coeff.mean_NOGO.multivar(count_channel_labels,count_i_mice) = mean(network_cluster_NOGO(channel_sub_network,:),2);
    clustering_coeff.var_NOGO.multivar(count_channel_labels,count_i_mice) = var(network_cluster_NOGO(channel_sub_network,:),0,2);
    % nodes conecting when count_channel_labels is connected to both 
    connected_nodes_multivar_NOGO{1,count_channel_labels}{1,count_i_mice}=clustered_nodes_NOGO;

    end
    
end
    
end % count_mice

%% --- plot connected nodes pool across sessions ---
 %{'GP'}    {'DG'}    {'VL'}    {'Cpu'}    {'CA1_Py'}    {'S1_bf'}    {'VM'}    {'LDVL'}    {'Rt'}    {'M1'}    {'CA1_Mol'}    {'dCA1'}
channel=7;

pool_of_regions_Hit=[];
pool_of_regions_CR=[];

for count_mice=1:size(connected_nodes_multivar_NOGO{1,1},2)
    
    pool_of_regions_Hit=[pool_of_regions_Hit connected_nodes_multivar_GO{1,channel}{1,count_mice}];
    pool_of_regions_CR=[pool_of_regions_CR connected_nodes_multivar_NOGO{1,channel}{1,count_mice}];

end % count_mice

%outkey='match';
%out = regexp(pool_of_regions_Hit,channel_label{1,2},outkey) 

subplot(4,1,1)
%[C,ia,ib] =intersect(channel_label,pool_of_regions_Hit, 'stable');
% find uniuqe elements and count them
clear region_labels ia ib region_name regional_histogram;
[region_labels,ia,ib] = unique(pool_of_regions_Hit, 'stable');

for jj=1:length(ia)
    regional_histogram(jj)=sum(ib==jj);  
end

% sort by the largest occurence
input_to_plot_Hit=[];
ind_sorted=[];
[input_to_plot_Hit, ind_sorted]=sort(regional_histogram,'descend');

% sort names and make categorical
region_name_Hit=[];
for jj=1:length(ia)
    region_name_Hit=[region_name_Hit region_labels(ind_sorted(jj))];
end
X= categorical(region_name_Hit);
X = reordercats(X,region_name_Hit);

% ---to normalize the bar plot by the number of observations---
observations=[];
observations_targets_Hit=[];
% ---calculate number of sessions across targets ---
i= channel;
%for i = 1:length(channel_label)
    
    % ---calculate presence of each channel across mice---

    input_count= clustering_coeff.mean_GO.multivar(i,:);
    input_count(isnan(input_count))=[];
    
    count_mice_in_channel(i)= size(input_count,2);
            
    % ---calculate if both source and target are present across mice---
        for count_targets=1:length(region_name_Hit)
            
            for count_mice =1:size(connected_nodes_sessions{1,i},2)

                if and(sum(strcmp(channel_label{1,i}, channel_labels_all{1,count_mice}))>0, ...
                    sum(strcmp(region_name_Hit{1,count_targets}, channel_labels_all{1,count_mice}))>0);

                    observations= [observations connected_nodes_sessions{1,i}{1,count_mice}];
                end
            
            end % count_mice
            % --- and sum up sessions for each target ---
            observations_targets_Hit(count_targets)= sum(observations);
            observations=[];
            
        end % count_targets
        
%end % count_sources

% ---and divide by the total number of expert sessions in each mouse---
bar(X,input_to_plot_Hit)
title('Connected to cluster Hit')
ylim([0 10])

subplot(4,1,2)
bar(X,input_to_plot_Hit./observations_targets_Hit)
title('Probability to observe connection Hit')
ylim([0 0.2])

subplot(4,1,3)
%[C,ia,ib] =intersect(channel_label,pool_of_regions_Hit, 'stable');
% find uniuqe elements and count them
clear region_labels ia ib region_name regional_histogram;
[region_labels,ia,ib] = unique(pool_of_regions_CR, 'stable');

for jj=1:length(ia)
    regional_histogram(jj)=sum(ib==jj);  
end

% sort by the largest occurence
input_to_plot_CR=[];
ind_sorted=[];
[input_to_plot_CR, ind_sorted]=sort(regional_histogram,'descend');

% sort names and make categorical
region_name_CR=[];
for jj=1:length(ia)
    region_name_CR=[region_name_CR region_labels(ind_sorted(jj))];
end
X= categorical(region_name_CR);
X = reordercats(X,region_name_CR);

bar(X,input_to_plot_CR)
title('Connected to cluster CR')
ylim([0 10])
% connecting_channels.labels.CR=X(input_to_plot_CR>=2);
% connecting_channels.strength.CR=input_to_plot(input_to_plot_CR>=2);

% ---to normalize the bar plot by the number of observations---
observations=[];
observations_targets_CR=[];
% ---calculate number of sessions across targets ---
i= channel;
%for i = 1:length(channel_label)
    
    % ---calculate presence of each channel across mice---

    input_count= clustering_coeff.mean_GO.multivar(i,:);
    input_count(isnan(input_count))=[];
    
    count_mice_in_channel(i)= size(input_count,2);
            
    % ---calculate if both source and target are present across mice---
        for count_targets=1:length(region_name_CR)
            
            for count_mice =1:size(connected_nodes_sessions{1,i},2)

                if and(sum(strcmp(channel_label{1,i}, channel_labels_all{1,count_mice}))>0, ...
                    sum(strcmp(region_name_CR{1,count_targets}, channel_labels_all{1,count_mice}))>0);

                    observations= [observations connected_nodes_sessions{1,i}{1,count_mice}];
                end
            
            end % count_mice
            % --- and sum up sessions for each target ---
            observations_targets_CR(count_targets)= sum(observations);
            observations=[];
            
        end % count_targets
        
%end % count_sources

subplot(4,1,4)
bar(X,input_to_plot_CR./observations_targets_CR)
title('Probability to observe connection CR')
ylim([0 0.2])

%% --- use frequent targets in region_name_Hit and region_name_CR to estimate them across mice ---

clear mean_in_mice_Hit mean_in_mice_CR
% repeat the loop across
for count_mice=1:size(connected_nodes_multivar_NOGO{1,1},2)
    
    % to quantify the standard deviation
    for jj=1:length(region_name_Hit)
        mean_in_mice_Hit(count_mice,jj)=mean(strcmp(region_name_Hit(jj), connected_nodes_multivar_GO{1,channel}{1,count_mice}));
    end
    
    for jj=1:length(region_name_CR)
        mean_in_mice_CR(count_mice,jj)=mean(strcmp(region_name_CR(jj), connected_nodes_multivar_NOGO{1,channel}{1,count_mice}));
    end

end % count_mice

input_count= mean_in_mice_Hit>0;
count_mice_in_channel_Hit= sum(input_count,1);

input_count= mean_in_mice_CR>0;
count_mice_in_channel_CR= sum(input_count,1);
%% --- concatenate for selected region across count_shuffles plot connected nodes ---
 
%{'GP'}    {'DG'}    {'VL'}    {'Cpu'}    {'CA1_Py'}    {'S1_bf'}    {'VM'}    {'LDVL'}    {'Rt'}    {'M1'}    {'CA1_Mol'}    {'dCA1'}

pool_of_shuffled_regions_Hit=[];
pool_of_shuffled_regions_CR=[];
shuffles_Hit = 0;
shuffles_CR = 0;

clear mean_in_shuffle_Hit mean_in_shuffle_CR

for count_mice=1:length(pool_mice)
        
        if isempty(connected_nodes_multivar_GO{1,channel}{1,count_mice})
        else
            for count_shuffles =1:number_of_shuffles
                pool_of_shuffled_regions_Hit=[pool_of_shuffled_regions_Hit connected_nodes_multivar_shuffle_GO{1,channel}{1,count_mice}{1,count_shuffles}];
                
                % increment for every iteration
                shuffles_Hit = shuffles_Hit+size(connected_nodes_multivar_shuffle_GO{1,channel}{1,count_mice},2);
                
                % to quantify the standard deviation
                for jj=1:length(region_name_Hit)
                        mean_in_shuffle_Hit(count_mice,count_shuffles,jj)=mean(strcmp(region_name_Hit(jj), connected_nodes_multivar_shuffle_GO{1,channel}{1,count_mice}{1,count_shuffles}));
                end
            end
        end
        
        if isempty(connected_nodes_multivar_NOGO{1,channel}{1,count_mice})
        else
            for count_shuffles =1:number_of_shuffles
                pool_of_shuffled_regions_CR=[pool_of_shuffled_regions_CR connected_nodes_multivar_shuffle_NOGO{1,channel}{1,count_mice}{1,count_shuffles}];
                shuffles_CR = shuffles_CR+size(connected_nodes_multivar_shuffle_NOGO{1,channel}{1,count_mice},2);
                
                % to quantify the standard deviation
                for jj=1:length(region_name_CR)
                        mean_in_shuffle_CR(count_mice,count_shuffles,jj)=mean(strcmp(region_name_CR(jj), connected_nodes_multivar_shuffle_NOGO{1,channel}{1,count_mice}{1,count_shuffles}));
                end
                
            end
        end
end

%%
% mean and var normalized by the number of observations (shuffles and mice)

subplot(2,1,1)
%[C,ia,ib] =intersect(channel_label,pool_of_regions_Hit, 'stable');
% find uniuqe elements and count them
clear region_labels ia ib region_name regional_histogram error_mice error_shuffle;

%choose only those present in region_name_Hit
mean_of_shuffle=[];
for jj=1:length(region_name_Hit)   
    mean_of_shuffle(jj)=sum(strcmp(region_name_Hit(jj), pool_of_shuffled_regions_Hit));
end

X= categorical(region_name_Hit);
X = reordercats(X,region_name_Hit);

%bar(X,[input_to_plot_Hit; mean_of_shuffle/number_of_shuffles])

error_mice=std(mean_in_mice_Hit,[],1);
error_shuffle=std(squeeze(mean(mean_in_shuffle_Hit,2)),[],1);
%errorbar(mean_of_shuffle/number_of_shuffles,error_shuffle);
barwitherr([((3*error_mice)./count_mice_in_channel_Hit)' ((3*error_shuffle)./count_mice_in_channel_Hit)'], ...
    [(input_to_plot_Hit./observations_targets_Hit)' (mean_of_shuffle./(number_of_shuffles*observations_targets_Hit))']);
xticks(1:length(input_to_plot_Hit))
set(gca,'XTickLabel',X)

ylim([0 0.5]);
title('Connected to cluster Hit')

% connecting_channels.labels.Hit=X(input_to_plot>=2);
% connecting_channels.strength.Hit=input_to_plot(input_to_plot>=2);

%subplot(4,1,2)
%test for significance
% for jj=1:length(region_name_Hit) 
%     p=ranksum(mean_in_mice_Hit(:,jj),reshape(mean_in_shuffle_Hit(:,:,jj),1,[]));
%     p_val(jj)=p;
% end
% p_val(p_val>0.05)=1;
% p_val(p_val<0.05)=0;
% imagesc(p_val);

subplot(2,1,2)

clear region_labels ia ib region_name regional_histogram error_mice error_shuffle;
%choose only those present in region_name_CR
mean_of_shuffle=[];
for jj=1:length(region_name_CR)   
    mean_of_shuffle(jj)=sum(strcmp(region_name_CR(jj), pool_of_shuffled_regions_CR));
end

X= categorical(region_name_CR);
X = reordercats(X,region_name_CR);

error_mice=std(mean_in_mice_CR,[],1);
mean_of_mice=sum(mean_in_mice_CR,1);
error_shuffle=std(squeeze(mean(mean_in_shuffle_CR,2)),[],1);
%errorbar(mean_of_shuffle/number_of_shuffles,error_shuffle);

barwitherr([((3*error_mice)./count_mice_in_channel_CR)' ((3*error_shuffle)./count_mice_in_channel_CR)'], ...
    [(input_to_plot_CR./observations_targets_CR)' (mean_of_shuffle./(number_of_shuffles*observations_targets_CR))']);
xticks(1:length(input_to_plot_CR))
set(gca,'XTickLabel',X)

ylim([0 0.5]);
title('Connected to cluster CR')


% subplot(4,1,4)
% %test for significance
% for jj=1:length(region_name_CR) 
%     p=ranksum(mean_in_mice_CR(:,jj),reshape(mean_in_shuffle_CR(:,:,jj),1,[]));
%     p_val(jj)=p;
% end
% p_val(p_val>0.05)=1;
% p_val(p_val<0.05)=0;
% imagesc(p_val);

%% plot shuffle in shaded error bar
error_shuffle=std(squeeze(mean(mean_in_shuffle_Hit,2)),[],1);
%choose only those present in region_name_Hit
mean_of_shuffle=[];
for jj=1:length(region_name_Hit)   
    mean_of_shuffle(jj)=sum(strcmp(region_name_Hit(jj), pool_of_shuffled_regions_Hit));
end
shadedErrorBar((1:length(region_name_Hit)),mean_of_shuffle/number_of_shuffles,3*error_shuffle);
ylim([0 10])
%% -----------------------
% --- Do circular plot ---
% ----------Hit-----------
input_size=length([0 connecting_channels.strength.Hit]);
input_adj_matrix= [[0 connecting_channels.strength.Hit]; zeros(input_size-1,input_size)];

input_labels=['RT' connecting_channels.labels.Hit];

myLabel = cell(length(input_labels));
for i = 1:length(input_labels)
  myLabel{i} = input_labels(i);
end

figure;
myColorMap = lines(length(input_adj_matrix));

circularGraph(input_adj_matrix,'Colormap',myColorMap,'Label',myLabel);

%% ------------------------
% --- Do circular plot ----
% -----------CR------------
input_size=length([0 connecting_channels.strength.CR]);
input_adj_matrix= [[0 connecting_channels.strength.CR]; zeros(input_size-1,input_size)];

input_labels=['RT' connecting_channels.labels.CR];

myLabel = cell(length(input_labels));
for i = 1:length(input_labels)
  myLabel{i} = input_labels(i);
end

figure;
myColorMap = lines(length(input_adj_matrix));

circularGraph(input_adj_matrix,'Colormap',myColorMap,'Label',myLabel);


%% -----------------------------------
% --- Hit vs CR boxplot plus lines ---
% ------------------------------------

number_of_mice= length(pool_mice);
% input_field_1=clustering_coeff.mean_GO.NIfTy;
% input_field_2=clustering_coeff.mean_NOGO.NIfTy;

input_field_1=clustering_coeff.mean_GO.multivar;
input_field_2=clustering_coeff.mean_NOGO.multivar;

% input_field_1=outdegree.mean_GO.multivar;
% input_field_2=outdegree.mean_NOGO.multivar;

% [p,tbl,stats]=anova1([input_field_1(1,:); input_field_2(1,:)]');
% multcompare(stats)
p_val=[];
for i=1:12
    [h,p] =ttest2(input_field_1(i,:), input_field_2(i,:));
    %[p,h] =signrank(input_field_1(i,:), input_field_2(i,:));
    p_val(i)=p;
end

input_1= std(input_field_1,0,2,'omitnan')./sqrt(number_of_mice);
input_2= mean(input_field_1,2,'omitnan');

input_3= std(input_field_2,0,2,'omitnan')./sqrt(number_of_mice);
input_4= mean(input_field_2,2,'omitnan');

% [p,tbl,stats]=anova1([input_field_1 input_field_2]');
% multcompare(stats)

h1=barwitherr([input_1(p_val<0.05) input_3(p_val<0.05)],[input_2(p_val<0.05) input_4(p_val<0.05)]);
set(gca,'XTickLabel',channel_label(p_val<0.05))
set(h1(1),'FaceColor',[0.0 0.0 0.6]);
%ylim([0 1]);
hold on
count_i=1;

for i=1:12

    if p_val(i)<0.05
        plot([count_i-0.15:0.3:count_i+0.15],[input_field_1(i,:); input_field_2(i,:)],'-o',...
            'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
            'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);
        
        count_i=count_i+1;
    end
end
% plot((1:12)-0.1,input_field_1','-o',...
%     'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
%     'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);
% 
% hold on
% plot((1:12)+0.1,input_field_2','-o',...
%     'LineWidth',1,'Color',[0.9,0.9,0.9], 'MarkerSize',5,...
%     'MarkerEdgeColor','r','MarkerFaceColor',[1,1,1]);
hold off
title(num2str(p));

% i=1;
% plot([i-0.1:0.2:0.1+i],[input_field_1(i,:); input_field_2(i,:)],'-o',...
%     'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
%     'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);

%% re-name input fields to save naive and expert clustering

% input_field_1_e=input_field_1; % go
% input_field_2_e=input_field_2; % nogo

% input_field_1_n=input_field_1;
% input_field_2_n=input_field_2;

%% ----------------------------------------------
% --- Hit: naive vs expert boxplot plus lines ---
% -----------------------------------------------

number_of_mice= length(pool_mice);
% plot only significant
plot_significant=1;

% --- Hit ---
input_field_1=input_field_1_n;
input_field_2=input_field_1_e;

p_val=[];
for i=1:12
    %[h,p] =ttest(input_field_1(i,:), input_field_2(i,:));
    [p,h] =ranksum(input_field_1(i,:), input_field_2(i,:));
    p_val(i)=p;
end

if plot_significant==0
    p_val=zeros(1,12);
end

input_1= std(input_field_1,0,2,'omitnan')./sqrt(number_of_mice);
input_2= mean(input_field_1,2,'omitnan');

input_3= std(input_field_2,0,2,'omitnan')./sqrt(number_of_mice);
input_4= mean(input_field_2,2,'omitnan');

subplot(2,1,1)
title('Hit')
h1=barwitherr([input_1(p_val<0.05) input_3(p_val<0.05)],[input_2(p_val<0.05) input_4(p_val<0.05)]);
set(gca,'XTickLabel',channel_label(p_val<0.05))
set(h1(1),'FaceColor',[0.0 0.0 0.6]);
%ylim([0 1]);
hold on
count_i=1;

for i=1:12

    if p_val(i)<0.05
        plot([count_i-0.15:0.3:count_i+0.15],[input_field_1(i,:); input_field_2(i,:)],'-o',...
            'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
            'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);
        
        count_i=count_i+1;
    end
end

hold off
title('Hit');

% --- CR ---
input_field_1=input_field_2_n;
input_field_2=input_field_2_e;

p_val=[];
for i=1:12
    [p,h] =ranksum(input_field_1(i,:), input_field_2(i,:));
    p_val(i)=p;
end

if plot_significant==0
    p_val=zeros(1,12);
end

input_1= std(input_field_1,0,2,'omitnan')./sqrt(number_of_mice);
input_2= mean(input_field_1,2,'omitnan');

input_3= std(input_field_2,0,2,'omitnan')./sqrt(number_of_mice);
input_4= mean(input_field_2,2,'omitnan');

subplot(2,1,2)
h2=barwitherr([input_1(p_val<0.05) input_3(p_val<0.05)],[input_2(p_val<0.05) input_4(p_val<0.05)]);
set(gca,'XTickLabel',channel_label(p_val<0.05))
set(h2(1),'FaceColor',[0.0 0.0 0.6]);
%ylim([0 1]);
hold on
count_i=1;

for i=1:12

    if p_val(i)<0.05
        plot([count_i-0.15:0.3:count_i+0.15],[input_field_1(i,:); input_field_2(i,:)],'-o',...
            'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
            'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);
        
        count_i=count_i+1;
    end
end

hold off
title('CR');

%% plot for all channels Hit vs CR

for i=1:12
ch=i;
subplot(4,3,i)
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
axis('square')
hold on
end
hold off

%% plot all channels naive vs expert
% input_field_1_e - expert go
% input_field_2_e - expert nogo

input_1= std(input_field_1_n,0,2,'omitnan')./sqrt(number_of_mice);
input_2= mean(input_field_1_n,2,'omitnan');

input_3= std(input_field_1_e,0,2,'omitnan')./sqrt(number_of_mice);
input_4= mean(input_field_1_e,2,'omitnan');


for i=1:12
ch=i;
subplot(4,3,i)
h1=barwitherr([input_1(ch,1) input_3(ch,1)],[input_2(ch,1) input_4(ch,1)]);
set(gca,'XTickLabel',{'naive','expert'})
set(h1(1),'FaceColor',[0.0 0.0 0.6]);
%ylim([0 1]);
hold on

plot(1:1:2,[input_field_1(ch,:); input_field_2(ch,:)],'-o',...
    'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
    'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);
ylim([0 0.3])
title(channel_label(ch))
axis('square')
hold on
end
hold off

%% plot for all channels naive vs expert

for i=1:12
ch=i;
subplot(4,3,i)
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
axis('square')
hold on
end
hold off

%% compare multivar on naive vs expert clustering coeff

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
    %[h,p] =ttest2(input_field_1_e(i,:), input_field_1_n(i,:));
    [p,h] =signrank(input_field_1_e(i,:), input_field_1_n(i,:));
    p_val_Hit(i)=p<0.05
end

for i=1:12
    %[h,p] =ttest2(input_field_2_e(i,:), input_field_2_n(i,:));
    [p,h] =signrank(input_field_2_e(i,:), input_field_2_n(i,:));
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

