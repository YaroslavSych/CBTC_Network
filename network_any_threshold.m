function [ network_adj_mat ] = network_any_threshold(p_table,time_steps,p_threshold)
%create unweighted adjacency matrix by thresholding p-values, taking OR over time steps
%   
                    network_size = size(p_table,2);

                    network_adj_mat = zeros(network_size);

                    for ch1 = 1:network_size
                        for ch2 = 1:network_size
                            if ch1 ~= ch2
                                network_adj_mat(ch1,ch2) = any(p_table(time_steps,ch1,ch2) < p_threshold, 1);
                            end
                        end
                    end
            
end

