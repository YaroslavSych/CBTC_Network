mice=length(mouse_IDs);

network_total_links_GO=[];
network_total_links_NOGO=[];    
network_shared_links_GO=[];
network_shared_links_NOGO=[];
network_shared_links_p_GO=[];
network_shared_links_p_NOGO=[];

performance_mice=cell(1,mice);
total_links=cell(1,mice);
shared_links=cell(1,mice);
network_links_shared_p_value=cell(1,mice);

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

pool_mice=[8 11 12 13];
for count_i_mice=1:length(pool_mice)

count_mice=pool_mice(count_i_mice);
session_names=sessions_IDs{1,count_mice};
sessions=length(session_names);

% --- import channel labels ---
channel_labels=channel_labels_all{1,count_mice};

% --- Choose time intervals (in seconds) ---
% CUE=[1,1.5], TEX=[2,3.5], LIK=[3.5,6]
%network_time_interval = {[1,1.5],[3,3.5],[6,6.5]};
network_time_interval = {[3,3.5]};

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
                network_time_steps = find(temp_data.parameteres.samples_timesteps >= min(network_time_interval{1,1}) & temp_data.parameteres.samples_timesteps < max(network_time_interval{1,1}));
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

% --- input to calculate network properties ---
adj_input_GO=network_adj_matrices.GO.bivar;
adj_input_NOGO=network_adj_matrices.NOGO.bivar;

% --- align all to learning session (first above 70% performance) --- 
high_perf_sessions = cell2mat(performance{1,count_mice}) > network_performance_threshold;
learning_ind =find(high_perf_sessions, 1, 'first');

% --- do alignment to learning onset ---
% create a vector +/-20 sessions around learning
length_to_pad=20;

% --- calculate network links and p-values of intersection ---
% --- Hit trials ---
network_links_GO=calculate_links_evolution(adj_input_GO);
total_links{1,count_i_mice}.Hit=...
    align_vectors_to_learning(network_links_GO{1,1}, learning_ind, length_to_pad);
shared_links{1,count_i_mice}.Hit=...
    align_vectors_to_learning(network_links_GO{1,2}, learning_ind, length_to_pad);
network_links_shared_p_value{1,count_i_mice}.Hit=...
    align_vectors_to_learning(network_links_GO{1,3}, learning_ind, length_to_pad);
% network_transitivity{1,count_i_mice}.Hit=...
%     align_vectors_to_learning(network_links_GO{1,4}, learning_ind, length_to_pad);
% network_assortativity{1,count_i_mice}.Hit=...
%     align_vectors_to_learning(network_links_GO{1,5}, learning_ind, length_to_pad);
% network_modularity{1,count_i_mice}.Hit=...
%     align_vectors_to_learning(network_links_NOGO{1,6}, learning_ind, length_to_pad);

% --- CR trials ---
network_links_NOGO=calculate_links_evolution(adj_input_NOGO);
total_links{1,count_i_mice}.CR=...
    align_vectors_to_learning(network_links_NOGO{1,1}, learning_ind, length_to_pad);
shared_links{1,count_i_mice}.CR=...
    align_vectors_to_learning(network_links_NOGO{1,2}, learning_ind, length_to_pad);
network_links_shared_p_value{1,count_i_mice}.CR=...
    align_vectors_to_learning(network_links_NOGO{1,3}, learning_ind, length_to_pad);
% network_transitivity{1,count_i_mice}.CR=...
%     align_vectors_to_learning(network_links_NOGO{1,4}, learning_ind, length_to_pad);
% network_assortativity{1,count_i_mice}.CR=...
%     align_vectors_to_learning(network_links_NOGO{1,5}, learning_ind, length_to_pad);
% network_modularity{1,count_i_mice}.CR=...
%     align_vectors_to_learning(network_links_NOGO{1,6}, learning_ind, length_to_pad);

% --- performance ---
performance_mice{1,count_i_mice}= ...
    align_vectors_to_learning(cell2mat(performance{1,count_mice}), learning_ind, length_to_pad);
end


%% plot Total links shaded error bar

mice=length(pool_mice);
input_Hit=[];
for i=1:mice
input_Hit= [input_Hit total_links{1,i}.Hit];
end
input_CR=[];
for i=1:mice
input_CR= [input_CR total_links{1,i}.CR];
end

mean_GO=mean(input_Hit,2,'omitnan');
std_GO=std(input_Hit,[],2,'omitnan')./sqrt(mice);

mean_NOGO=mean(input_CR,2,'omitnan');
std_NOGO=std(input_CR,[],2,'omitnan')./sqrt(mice);

mean_performance=mean(cell2mat(performance_mice),2,'omitnan');
std_performance=std(cell2mat(performance_mice),[],2,'omitnan')./sqrt(mice);

days=1:41;
subplot(2,1,1)
shadedErrorBar(days,mean_performance,std_performance,'b',0)
%ylim([0 1])
hold on
subplot(2,1,2)
shadedErrorBar(days,mean_GO,std_GO,'b',0);
%ylim([0 20])
hold on
shadedErrorBar(days,mean_NOGO,std_NOGO,'r',0);
%ylim([0 20])
hold off

%% plot Shared links shaded error bar

input_Hit_shared=[];
for i=1:mice
input_Hit_shared= [input_Hit_shared shared_links{1,i}.Hit];
end
input_CR_shared=[];
for i=1:mice
input_CR_shared= [input_CR_shared shared_links{1,i}.CR];
end

mean_GO=mean(input_Hit_shared,2,'omitnan');
std_GO=std(input_Hit_shared,[],2,'omitnan')./sqrt(mice);

mean_NOGO=mean(input_CR_shared,2,'omitnan');
std_NOGO=std(input_CR_shared,[],2,'omitnan')./sqrt(mice);

mean_performance=mean(cell2mat(performance_mice),2,'omitnan');
std_performance=std(cell2mat(performance_mice),[],2,'omitnan')./sqrt(mice);

days=1:41;
subplot(2,1,1)
shadedErrorBar(days,mean_performance,std_performance,'b',0)
%ylim([0 1])
hold on
subplot(2,1,2)
shadedErrorBar(days,mean_GO,std_GO,'b',0);
%ylim([0 20])
hold on
shadedErrorBar(days,mean_NOGO,std_NOGO,'r',0);
%ylim([0 20])
hold off

%% plot precent of Shared links in Total links shaded error bar
input_fraction_Hit= (100*input_Hit_shared)./input_Hit;
input_fraction_CR= (100*input_CR_shared)./input_CR;

mean_GO=mean(input_fraction_Hit,2,'omitnan');
std_GO=std(input_fraction_Hit,[],2,'omitnan')./sqrt(mice);

mean_NOGO=mean(input_fraction_CR,2,'omitnan');
std_NOGO=std(input_fraction_CR,[],2,'omitnan')./sqrt(mice);

mean_performance=mean(cell2mat(performance_mice),2,'omitnan');
std_performance=std(cell2mat(performance_mice),[],2,'omitnan')./sqrt(mice);

days=1:41;
subplot(2,1,1)
shadedErrorBar(days,mean_performance,std_performance,'b',0)
%ylim([0 1])
hold on
subplot(2,1,2)
shadedErrorBar(days,mean_GO,std_GO,'b',0);
%ylim([0 20])
hold on
shadedErrorBar(days,mean_NOGO,std_NOGO,'r',0);
%ylim([0 20])
hold off

%% plot assorativity/transitivity shaded error bar
% network_assortativity{1,i} or network_transitivity{1,i}
input_Hit=[];
for i=1:mice
input_Hit= [input_Hit network_assortativity{1,i}.Hit];
end
input_CR=[];
for i=1:mice
input_CR= [input_CR network_assortativity{1,i}.CR];
end

mean_GO=mean(input_Hit,2,'omitnan');
std_GO=std(input_Hit,[],2,'omitnan')./sqrt(mice);

mean_NOGO=mean(input_CR,2,'omitnan');
std_NOGO=std(input_CR,[],2,'omitnan')./sqrt(mice);

mean_performance=mean(cell2mat(performance_mice),2,'omitnan');
std_performance=std(cell2mat(performance_mice),[],2,'omitnan')./sqrt(mice);

days=1:41;
subplot(2,1,1)
shadedErrorBar(days,mean_performance,std_performance,'b',0)
%ylim([0 1])
hold on
subplot(2,1,2)
shadedErrorBar(days,mean_GO,std_GO,'b',0);
%ylim([0 20])
hold on
shadedErrorBar(days,mean_NOGO,std_NOGO,'r',0);
%ylim([0 20])
hold off

%% GO scatter plot of probability to generate shared links at random

[dim1,dim2]=size(cell2mat(performance_mice));
input_1 = reshape(cell2mat(performance_mice),[1 dim1*dim2]);

input_Hit=[];
for i=1:mice
input_Hit= [input_Hit network_links_shared_p_value{1,i}.Hit];
end
input_CR=[];
for i=1:mice
input_CR= [input_CR network_links_shared_p_value{1,i}.CR];
end

input_2 = reshape(input_Hit,[1 dim1*dim2]);

% plot size as scaled by the number of shared links normalized by
% total number of possible links N*(N-1)
N=12;
all_links_possible=N*(N-1);
sz = 100*reshape(input_Hit_shared,[1 dim1*dim2])/all_links_possible;

%plot the color as the total number links
c = 10*reshape(input_Hit_shared,[1 dim1*dim2])/all_links_possible;
% if to color plot by animal use linspace(1,10,dim1*dim2);

mdl = fitlm(table(input_1' , input_2'),'linear')
scatter(input_1,input_2,sz,c,'filled','MarkerEdgeColor',[0.1 .1 .1],'LineWidth',0.5)
ylim([0 0.2])
hold on
plot(mdl)
% hold on
% hline = refline([0 0.05]);
% hline.Color = 'b';
hold off

%% NOGO scatter plot of to generate shared links at random

input_2 = reshape(input_CR,[1 dim1*dim2]);

% plot size as scaled by the number of shared links normalized by
% total number of possible links N*(N-1)
N=12;
all_links_possible=N*(N-1);
sz = 200*((reshape(input_CR_shared,[1 dim1*dim2])+1)/all_links_possible);

%plot the color as the total number links
c = reshape(input_CR_shared,[1 dim1*dim2])/all_links_possible;
% if to color plot by animal use linspace(1,10,dim1*dim2);

mdl = fitlm(table(input_1' , input_2'),'linear')
scatter(input_1,input_2,sz,c,'filled','MarkerEdgeColor',[0.1 .1 .1],'LineWidth',0.5)
ylim([0 0.2])
hold on
plot(mdl)
% hold on
% hline = refline([0 0.05]);
% hline.Color = 'b';
hold off

%% save total number of links and LCC

cd(paths_mice{1,count_mice});
mydir  = cd(paths_mice{1,count_mice});
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1)

full_path=[newdir '\network_links_lcc_learning.mat'];
    %filename = 'network_links_lcc_learning.mat';
    save(full_path, 'network_links_GO_low', 'network_links_GO_high', ...
        'network_links_NOGO_low', 'network_links_NOGO_high', 'lcc_GO_size_low', 'lcc_GO_size_high', ...
        'lcc_NOGO_size_low', 'lcc_NOGO_size_high');
    
% for count_mice=1:mice
%     
%     cd(paths_mice{1,count_mice})
%     full_path=[paths_mice{1,count_mice} '\network_links_lcc_learning.mat'];
%     %filename = 'network_links_lcc_learning.mat';
%     save(full_path, 'network_links_GO_low', 'network_links_GO_high', ...
%         'network_links_NOGO_low', 'network_links_NOGO_high', 'lcc_GO_size_low', 'lcc_GO_size_high', ...
%         'lcc_NOGO_size_low', 'lcc_NOGO_size_high');
%     
% end

