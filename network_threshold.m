function [ network_adj_mat ] = network_threshold(p_table,time_steps,p_threshold)
%create weighted adjacency matrix by thresholding p-values, averaging over time steps, and cutting low-frequency links
%   
    
    network_size = size(p_table,1);
    
    network_adj_mat = zeros(network_size);
     
    for ch1 = 1:network_size
        for ch2 = 1:network_size
            if ch1 ~= ch2
                network_adj_mat(ch1,ch2) = mean(p_table{ch1,ch2}{1}(time_steps) < p_threshold);
            end
        end
    end 
            
end

