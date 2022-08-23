function reorder_input=reorder_adj_matrix_by_channel_name(temp_input,all_sessions,count_source_label,channel_labels_per_mouse,channel_label) 
% input: 1. input adjacency matrix, 2. all sessions provides session
% indexes after reading across bivar and multivar adjacency matrices
% 3. channel index of the source= from which the out-going links should be aligned, 
% 4. channel labels in the current mouse,
% 5. template of channel lables as provided by the
% channel_labels variable

% output: output_aligned a vector of weighted by appearance frequency links

% generate empty template
reorder_input=[];

% ---find a channel label from channel_label template---
% and call its index channel_index
for count_channel_labels=1:length(channel_label)
                    
    label = strncmp(channel_labels_per_mouse,channel_label{1,count_channel_labels},5);
    channel_sub_network = find(label>0);
    
    if length(channel_sub_network)>1
        channel_index(count_channel_labels)=channel_sub_network(1,1);
    elseif isempty(channel_sub_network)
        % --- give here NaN to index ----
         channel_index(count_channel_labels)=NaN;
    else
        channel_index(count_channel_labels)=channel_sub_network;
    end
end

for count_target=1:length(channel_index)
    
    if count_source_label==count_target
        
    for count_sessions = 1:length(all_sessions)    
        reorder_input(count_sessions,count_target)=NaN;
    end

    else

    for count_sessions = 1:length(all_sessions)
        
        % channel_index(count_source_label) is the index for the source
        % channel_index(count_target) is the index for the target
        if or(isnan(channel_index(count_source_label)),isnan(channel_index(count_target)))
        reorder_input(count_sessions,count_target)=NaN;
        else
        reorder_input(count_sessions,count_target)=temp_input(channel_index(count_source_label),channel_index(count_target),all_sessions(count_sessions));
    end
    
    end
    
end

end 

end