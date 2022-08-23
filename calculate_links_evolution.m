function network_links_GO=calculate_links_evolution(network_adj_matrices_GO)
    
% Trial type GO
    count_network_time_interval=1;
    % find number of nodes and links
    network_size = size(network_adj_matrices_GO,1);
    network_links_total = network_size^2 - network_size;
    % find indices of weights above threshold
    links_indices = network_adj_matrices_GO > 0;
    % compute number of links
    network_links = squeeze(sum(sum(links_indices,1),2));
    
    % compute number of shared links between consecutive sessions
    indices_sum = links_indices(:,:,1:end-1) + links_indices(:,:,2:end);
    network_links_shared = squeeze(sum(sum(indices_sum == 2,1),2));
    % compute p-value of shared links
    network_links_shared_p_value = NaN([length(network_links_shared),1]);
    warning('off','MATLAB:nchoosek:LargeCoefficient');
    sessions=size(network_adj_matrices_GO,3);
    for count = 1:sessions
        network_links1 = network_links(count);
        
        if sessions==1
            network_links2 =NaN;
            network_links_shared_p_value(count) =NaN;
        elseif (count+1)<=sessions
            network_links2 = network_links(count+1);
            network_links_shared_p_value(count) = nchoosek(network_links_total,network_links_shared(count))*nchoosek(network_links_total-network_links_shared(count),network_links1-network_links_shared(count))*nchoosek(network_links_total-network_links1,network_links2-network_links_shared(count))/(nchoosek(network_links_total,network_links1)*nchoosek(network_links_total,network_links2));
        end
        
%         transitivity(count)=transitivity_wd(network_adj_matrices_GO(:,:,count));
%         assortativity(count)=assortativity_wei(network_adj_matrices_GO(:,:,count),3);
%         [Ci Q] = modularity_dir(network_adj_matrices_GO(:,:,count));
%         Q_max(count)=Q;
%         agreement_coeff(count)=agreement(Ci);
    end
    
network_links_GO{count_network_time_interval,1}=network_links;
network_links_GO{count_network_time_interval,2}=network_links_shared;
network_links_GO{count_network_time_interval,3}=network_links_shared_p_value;
% network_links_GO{count_network_time_interval,4}=transitivity;
% network_links_GO{count_network_time_interval,5}=assortativity;
% network_links_GO{count_network_time_interval,6}=Q_max;
% network_links_GO{count_network_time_interval,7}=agreement_coeff;

end