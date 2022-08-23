performance_threshold=0.7;
count_predictors=12;

p_joined=[];

mice=length(mouse_IDs);
learning_ind=cell(1,mice);

% ------ preallocate cell arrays for correlation matrixes  ------

for count_mice = 1:mice
        
        sessions = length(sessions_IDs{1,count_mice});
        
        low_perf_sessions_indices{1,count_mice} = cell(1,sessions);
        high_perf_sessions_indices{1,count_mice} = cell(1,sessions);
        
        % Hit vs CR discrimination 
        ACC_stimulus{1,count_mice}= cell(1,sessions);
        
        % Hit vs MISS discrimination 
        ACC_choice{1,count_mice}= cell(1,sessions);
        
        % CR vs FA discrimination 
        ACC_action{1,count_mice}= cell(1,sessions);
        
end

% ------ loop across mice  ------
for count_mice=1:mice
    
    
    sessions = length(sessions_IDs{1,count_mice});

    high_perf_sessions_indices{1,count_mice} = find(cell2mat(performance{1,count_mice})...
    > performance_threshold);
    low_perf_sessions_indices{1,count_mice} = find(cell2mat(performance{1,count_mice})...
    <= performance_threshold);

    learning_ind{1,count_mice} =find(cell2mat(performance{1,count_mice})...
    > performance_threshold, 1, 'first');

    % ------ loop across sessions to collect whisk envelope ------
    
        % ------ loop across sessions to normalize calcium data ------

    for count_sessions=1:sessions

        if isempty(data{1,count_mice}{1,count_sessions}.data)==1
            
            % make ratio vectors of the same size (channels,1)
            % Hit vs CR discrimination 
            ACC_stimulus{1,count_mice}{1,count_sessions}.mean= NaN*ones(count_predictors,1);
            ACC_stimulus{1,count_mice}{1,count_sessions}.var= NaN*ones(count_predictors,1);

            % Hit vs MISS discrimination 
            ACC_choice{1,count_mice}{1,count_sessions}.mean= NaN*ones(count_predictors,1);
            ACC_choice{1,count_mice}{1,count_sessions}.var= NaN*ones(count_predictors,1);

            % CR vs FA discrimination 
            ACC_action{1,count_mice}{1,count_sessions}.mean= NaN*ones(count_predictors,1);
            ACC_action{1,count_mice}{1,count_sessions}.var= NaN*ones(count_predictors,1);

        else
            
            % ------ read all trials  ------
            Hit_trials = behaviorvar{1, count_mice}{1, count_sessions}.iGO;
            CR_trials = behaviorvar{1, count_mice}{1, count_sessions}.iNOGO;
            
            Miss_trials = behaviorvar{1, count_mice}{1, count_sessions}.iMISS;
            FA_trials = behaviorvar{1, count_mice}{1, count_sessions}.iFA;

            % check for the size of dataset 
            [fn,~,~]=size(data{1, count_mice}{1, count_sessions}.data);
            Hit_trials = Hit_trials(Hit_trials<fn);
            CR_trials = CR_trials(CR_trials<fn);
            Miss_trials = Miss_trials(Miss_trials<fn);
            FA_trials = FA_trials(FA_trials<fn);
            
            % remove CR trials if preceeded by 2 FA
%             [CR1,i_cr1,i_fa1] = intersect(CR_trials,FA_trials+1);
%             [CR2,i_cr2,i_fa2] = intersect(CR_trials,FA_trials+2);
%             CR_trials(intersect(i_cr1,i_cr2))=[];
            
            
            test_Hit = data{1,count_mice}{1,count_sessions}.data(Hit_trials,:,:);
            test_CR = data{1,count_mice}{1,count_sessions}.data(CR_trials,:,:);
            test_Miss = data{1,count_mice}{1,count_sessions}.data(Miss_trials,:,:);
            test_FA = data{1,count_mice}{1,count_sessions}.data(FA_trials,:,:);
            
            stimulus= 60:70; % from 2.5 to 3.5 seconds of a trial time
            %action= 116:136; % from 5 to 6.5 seconds of a trial time
            
            % --- Hit vs CR discrimination ---
            % use 12 channels as predictors and a pool of 12 channels as
            % 13th predictor
            count_predictors=size(test_Hit,3);
            limit_by_trials=40;
            ACC1=[];
            ACC1_var=[];  
            
            for count_channel=1:count_predictors 
                [ACC,ACC_var,f1]= calculate_ACC(test_Hit,test_CR,stimulus,count_channel,limit_by_trials);
                ACC1=[ACC1 ACC];
                ACC1_var=[ACC1_var ACC_var];
            end
            
            ACC_stimulus{1,count_mice}{1,count_sessions}.mean= ACC1;
            ACC_stimulus{1,count_mice}{1,count_sessions}.var= ACC1_var;
            ACC1=[];
            ACC1_var=[];
            
            % --- Hit vs FA discrimination ---

                for count_channel=1:count_predictors
                    [ACC,ACC_var,f1]= calculate_ACC(test_Hit,test_FA,stimulus,count_channel,limit_by_trials);
                    ACC1=[ACC1 ACC];
                    ACC1_var=[ACC1_var ACC_var];
                end
            
            ACC_choice{1,count_mice}{1,count_sessions}.mean= ACC1;
            ACC_choice{1,count_mice}{1,count_sessions}.var= ACC1_var;
            ACC1=[];
            ACC1_var=[];
            
            % --- CR vs FA discrimination ---
            
            for count_channel=1:count_predictors  
                [ACC,ACC_var,f1]= calculate_ACC(test_FA,test_CR,stimulus,count_channel,limit_by_trials);
                ACC1=[ACC1 ACC];
                ACC1_var=[ACC1_var ACC_var];
            end

            ACC_action{1,count_mice}{1,count_sessions}.mean= ACC1;
            ACC_action{1,count_mice}{1,count_sessions}.var= ACC1_var;
            ACC1=[];
            ACC1_var=[];
            
            end
            
        end
end

%% ------ sort ACC by channel ------
mice=length(mouse_IDs);

% add the "all" channel #13 to all mice/channel labels
% [~,d1]=size(channel_labels_all);
% for i=1:d1
%     channel_labels_all{1,i}(13)={'all'};
% end

%channel_label = {'Cpu','GP','DG','CA1_Py','VL','VM','LDVL','Rt','M1','S1_bf','all'};
channel_label = {'GP', 'Cpu','VL', 'VM', 'LDVL', 'Rt', 'M1', 'S1_bf', 'DG', 'CA1_Py', 'CA1_Mol','dCA1'};
%channel_label = flip(channel_label);

ch=length(channel_label);
label_ind=cell(1,ch);
k=1;

for count_channel_labels=1:length(channel_label)
    for count_mice=1:mice
        label=[];
        label = strncmp(channel_labels_all{count_mice},channel_label{1,count_channel_labels},5);
        channel_sub_network = find(label>0);
        if length(channel_sub_network)>1
            labels_mice(count_channel_labels,count_mice,1)=channel_sub_network(1,1);
            labels_mice(count_channel_labels,count_mice,2)=channel_sub_network(2,1);
        elseif isempty(channel_sub_network) 
            labels_mice(count_channel_labels,count_mice,1)=NaN;
            labels_mice(count_channel_labels,count_mice,2)=NaN;
        else
            labels_mice(count_channel_labels,count_mice,1)=channel_sub_network;
            labels_mice(count_channel_labels,count_mice,2)=NaN;
        end
    end
end

labels_mice_regional=reshape(labels_mice,12,[]);
% increase number of mice by 2 due to doubled recordings from some regions
% put NaNs otherwise
mice=size(labels_mice_regional,2);

% preallocate stat_label/learning
for count_mice=1:mice

        pool_ACC_choice{1,count_mice}=cell(1,ch);
        pool_ACC_stimulus{1,count_mice}=cell(1,ch);
        pool_ACC_action{1,count_mice}=cell(1,ch);
        
end

test=[];
count_mice_vector=[1:14 1:14];

for count_label=1:ch
    
    for count_i_mice=1:mice
    
    label=labels_mice_regional(count_label,count_i_mice);
    count_mice= count_mice_vector(count_i_mice);
    
    if isnan(label)
        
        sessions = length(sessions_IDs{1,count_mice});
        pool_ACC_choice{1,count_i_mice}{1,count_label}=NaN(1,sessions);
        pool_ACC_stimulus{1,count_i_mice}{1,count_label}=NaN(1,sessions);
        pool_ACC_action{1,count_i_mice}{1,count_label}=NaN(1,sessions);
            
    else 
        
        
        sessions = length(sessions_IDs{1,count_mice});
        for count_sessions=1:sessions
            
            temp_data=ACC_choice{1,count_mice}{1,count_sessions}.mean;
            temp_ACC_choice(count_sessions)=temp_data(label);
            
            temp_data=ACC_stimulus{1,count_mice}{1,count_sessions}.mean;
            temp_ACC_stimulus(count_sessions)=temp_data(label);
            
            temp_data=ACC_action{1,count_mice}{1,count_sessions}.mean;
            temp_ACC_action(count_sessions)=temp_data(label);
            
        end
            pool_ACC_choice{1,count_i_mice}{1,count_label}=temp_ACC_choice;
            pool_ACC_stimulus{1,count_i_mice}{1,count_label}=temp_ACC_stimulus;
            pool_ACC_action{1,count_i_mice}{1,count_label}=temp_ACC_action;
            
            temp_ACC_choice=[];
            temp_ACC_stimulus=[];
            temp_ACC_action=[];
        end
    
    end

    
end

% --- align ACC to learning session ---

% align all to learning session (first above 75% performance) 
% temp_loaded_data.performance

% ------ loop across mice  ------
learning_ind=[];

pool_mice=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 1 2 3 4 5 6 7 8 9 10 11 12 13 14];


for i=1:length(pool_mice)
    
    count_mice=pool_mice(i);
    
    sessions = length(sessions_IDs{1,count_mice});

    learning_ind(i) =find(double(cell2mat(performance{1,count_mice})> performance_threshold), 1, 'first');
    
    % create a vector +/-20 sessions around learning
    length_to_pad=20;
    
    for count_label=1:ch
    
        pool_ACC_choice{1,i}{1,count_label}= align_vectors_to_learning(pool_ACC_choice{1,count_mice}{1,count_label},...
        learning_ind(i), length_to_pad);
    
        pool_ACC_stimulus{1,i}{1,count_label}= align_vectors_to_learning(pool_ACC_stimulus{1,count_mice}{1,count_label},...
        learning_ind(i), length_to_pad);
    
        pool_ACC_action{1,i}{1,count_label}= align_vectors_to_learning(pool_ACC_action{1,count_mice}{1,count_label},...
        learning_ind(i), length_to_pad);
    
    end

end


% ---do a matrix plot for all regions ---

% new fig channel labels 
channel_label = {'GP', 'Cpu','VL', 'VM', 'LDVL', 'Rt', 'M1', 'S1_bf', 'DG', 'CA1_Py', 'CA1_Mol','dCA1'};
%channel_label = flip(channel_label);

mean_ACC_choice=[];
std_ACC_choice=[];
mean_ACC_stimulus=[];
std_ACC_stimulus=[];
mean_ACC_action=[];
std_ACC_action=[];

% rename pool_mice according to the counter of the pooled mice
pool_mice=1:28;

for count_channel=1:length(channel_label)
    
    all_ACC_choice= [];
    all_ACC_stimulus= [];
    all_ACC_action= [];

    for i=1:length(pool_mice)
        
        count_mice=pool_mice(i);

        all_ACC_choice= [all_ACC_choice pool_ACC_choice{1,count_mice}{1,count_channel}];
        all_ACC_stimulus= [all_ACC_stimulus pool_ACC_stimulus{1,count_mice}{1,count_channel}];
        all_ACC_action= [all_ACC_action pool_ACC_action{1,count_mice}{1,count_channel}];

    end

% all_ACC_choice(all_ACC_choice == 0) = NaN;
% all_ACC_stimulus(all_ACC_stimulus == 0) = NaN;
% all_ACC_action(all_ACC_action == 0) = NaN;

mean_ACC_choice(count_channel,:)=mean(all_ACC_choice,2,'omitnan')';
std_ACC_choice(count_channel,:)=std(all_ACC_choice,[],2,'omitnan')';

mean_ACC_stimulus(count_channel,:)=mean(all_ACC_stimulus,2,'omitnan')';
std_ACC_stimulus(count_channel,:)=std(all_ACC_stimulus,[],2,'omitnan')';

mean_ACC_action(count_channel,:)=mean(all_ACC_action,2,'omitnan')';
std_ACC_action(count_channel,:)=std(all_ACC_action,[],2,'omitnan')';

end

%
% input_to_plot=mean_ACC_action;
% th = 2*pi/length(channel_label):2*pi/length(channel_label):2*pi;
% 
% test_angles=repmat(th',1,length(days));
% 
% for every_channel=1:length(channel_label)
%     
%     angle = test_angles(every_channel,:);
%     r= [days];
%     sz = 2000*(input_to_plot(every_channel,:)-0.5);
%     c = 100*(input_to_plot(every_channel,:)-0.5);
%     polarscatter(angle,r,sz,c,'filled','MarkerFaceAlpha',.5)
%     hold on
%     
% end
% 
% hold off

test_channels=repmat([1:length(channel_label)]',1,length(days))/4;

input_days=1:41;
start_day=10;
end_day=30;
scale_size=100;
shift_input=0.5;
scale_input=0.4;
tick_step=0.1;
all_channels=(1:12)/4;
xchannel=3.1;

subplot(3,1,1)
input_to_plot=mean_ACC_choice;
input_std_to_plot=std_ACC_choice;
f1=[];

for every_channel=1:length(channel_label)
    
    sz = scale_size*(input_to_plot(every_channel,:)-scale_input);
    delta= (input_to_plot(every_channel,:)-shift_input);
    c = input_std_to_plot(every_channel,:);
    
    % fit polynomial
    input_to_fit=test_channels(every_channel,:)+ delta;
    input_to_fit(isnan(input_to_fit))=test_channels(every_channel);
    p = polyfit(input_days,input_to_fit,5);
    f1(every_channel,:) = polyval(p,input_days);
        
    plot(f1(every_channel,:),input_days,'k--')
    hold on
    
    scatter(test_channels(every_channel,:)+ delta,input_days,sz,c,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.7 .7 .7],...
              'LineWidth',0.5)
    hold on

end

ylim([start_day end_day])
xlim([0 xchannel])
    ax = gca;
    ax.XTick = sort([ all_channels ,all_channels+tick_step],'ascend');
    title('Choice')
hold off

subplot(3,1,2)
input_to_plot=mean_ACC_stimulus;
input_std_to_plot=std_ACC_stimulus;
f1=[];

for every_channel=1:length(channel_label)
    
    sz = scale_size*(input_to_plot(every_channel,:)-scale_input);
    delta= (input_to_plot(every_channel,:)-shift_input);
    c = input_std_to_plot(every_channel,:);
    
    % fit polynomial
    input_to_fit=test_channels(every_channel,:)+ delta;
    input_to_fit(isnan(input_to_fit))=test_channels(every_channel);
    p = polyfit(input_days,input_to_fit,5);
    f1(every_channel,:) = polyval(p,input_days);
    
    plot(f1(every_channel,:),input_days,'k--')
    hold on
    
    scatter(test_channels(every_channel,:)+ delta,input_days,sz,c,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.7 .7 .7],...
              'LineWidth',0.5)
    hold on
    

end

ylim([start_day end_day])
xlim([0 xchannel])
    ax = gca;
    ax.XTick = sort([ all_channels ,all_channels+tick_step],'ascend');
%     ax.XTickLabel = channel_label;
    title('Stimulus')
    
hold off

subplot(3,1,3)
input_to_plot=mean_ACC_action;
input_std_to_plot=std_ACC_action;
f1=[];


for every_channel=1:length(channel_label)
    
    sz = scale_size*(input_to_plot(every_channel,:)-scale_input);
    delta= (input_to_plot(every_channel,:)-shift_input);
    c = input_std_to_plot(every_channel,:);
    
    % fit polynomial
    input_to_fit=test_channels(every_channel,:)+ delta;
    input_to_fit(isnan(input_to_fit))=test_channels(every_channel);
    p = polyfit(input_days,input_to_fit,5);
    f1(every_channel,:) = polyval(p,input_days);
    
    plot(f1(every_channel,:),input_days,'k--')
    hold on
    
    scatter(test_channels(every_channel,:)+delta,input_days,sz,c,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.7 .7 .7],...
              'LineWidth',0.5)
    hold on
    

end

ylim([start_day end_day])
xlim([0 xchannel])
    ax = gca;
    ax.XTick = sort([ all_channels ,all_channels+tick_step],'ascend');
%     tick_labels = string(ax.XTickLabel);
%     tick_labels(1:2:end) = nan; % remove every other one
%     ax.XTickLabel = channel_label;
    title('Action')
    
hold off

%% Ranksum test all_ACC_stimulus and action for significance

for count_channel=1:length(channel_label)
    
    all_ACC_choice= [];
    all_ACC_stimulus= [];
    all_ACC_action= [];

    for i=1:length(pool_mice)
        
        count_mice=pool_mice(i);

        all_ACC_choice= [all_ACC_choice pool_ACC_choice{1,count_mice}{1,count_channel}];
        all_ACC_stimulus= [all_ACC_stimulus pool_ACC_stimulus{1,count_mice}{1,count_channel}];
        all_ACC_action= [all_ACC_action pool_ACC_action{1,count_mice}{1,count_channel}];

    end
    % test naive vs expert
    input_1= mean(all_ACC_stimulus(1:18,:),1,'omitnan');
    %input_1=reshape(all_ACC_stimulus(1:19,:),1,[]);
    input_2= mean(all_ACC_stimulus(23:end,:),1,'omitnan');
    %input_2= reshape(all_ACC_stimulus(21:end,:),1,[]);
    [p,h] = ranksum(input_1,input_2);
    p_val.stimulus(count_channel)=p;
    
    input_1= mean(all_ACC_action(20:21,:),1,'omitnan');
    %input_1=reshape(all_ACC_action(20:21,:),1,[]);
    input_2= mean(all_ACC_action([1:19, 22:end],:),1,'omitnan');
    %input_2= reshape(all_ACC_action([1:19, 22:end],:),1,[]);
    
%     input_1=mean(all_ACC_action(20:21,:),1,'omitnan');
%     input_2= mean(all_ACC_action([1:19, 22:end],:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.action(count_channel)=p;
    
    %input_1= mean(all_ACC_choice(20,:),1,'omitnan');
    input_1=reshape(all_ACC_choice(19:21,:),1,[]);
    %input_2= mean(all_ACC_choice([1:19, 21:end],:),1,'omitnan');
    input_2= reshape(all_ACC_choice([1:18, 22:end],:),1,[]);
    [p,h] = ranksum(input_1,input_2);
    p_val.choice(count_channel)=p;
    
    % test learning +/-3 day vs late expert (3 days from each side) 
    input_1= mean(all_ACC_stimulus(17:23,:),1,'omitnan');
    input_2= mean(all_ACC_stimulus([5:15, 25:end],:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.learning_stimulus(count_channel)=p;
    
    input_1= mean(all_ACC_action(17:23,:),1,'omitnan');
    input_2= mean(all_ACC_action([5:15, 25:end],:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.learning_action(count_channel)=p;
    
    input_1= mean(all_ACC_choice(17:23,:),1,'omitnan');
    input_2= mean(all_ACC_choice([5:15, 25:end],:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.learning_choice(count_channel)=p;
    
end


%% concatenate across mice

n_ch=4;
all_ACC_choice= [];
all_ACC_stimulus= [];
all_ACC_action= [];

for i=1:length(pool_mice)

    count_mice=pool_mice(i);
    
    all_ACC_choice= [all_ACC_choice pool_ACC_choice{1,count_mice}{1,n_ch}];
    all_ACC_stimulus= [all_ACC_stimulus pool_ACC_stimulus{1,count_mice}{1,n_ch}];
    all_ACC_action= [all_ACC_action pool_ACC_action{1,count_mice}{1,n_ch}];
 
end

all_ACC_choice(all_ACC_choice == 0) = NaN;
all_ACC_stimulus(all_ACC_stimulus == 0) = NaN;
all_ACC_action(all_ACC_action == 0) = NaN;

mean_ACC_choice=mean(all_ACC_choice,2,'omitnan');

N_choice=count_if_number(all_ACC_choice);
N_choice(N_choice==0)=NaN;
std_ACC_choice=std(all_ACC_choice,[],2,'omitnan')./sqrt(N_choice)';

mean_ACC_stimulus=mean(all_ACC_stimulus,2,'omitnan');

N_stimulus=count_if_number(all_ACC_stimulus);
N_stimulus(N_stimulus==0)=NaN;
std_ACC_stimulus=std(all_ACC_stimulus,[],2,'omitnan')./sqrt(N_stimulus)';

mean_ACC_action=mean(all_ACC_action,2,'omitnan');

N_action=count_if_number(all_ACC_action);
N_action(N_action==0)=NaN;
std_ACC_action=std(all_ACC_action,[],2,'omitnan')./sqrt(N_action)';

days=1:41;

shadedErrorBar(days,mean_ACC_choice,std_ACC_choice,'k',0);
ylim([0.5 1])
xlim([10 32])
hold on
shadedErrorBar(days,mean_ACC_stimulus,std_ACC_stimulus,'b',0)
ylim([0.5 1])
xlim([10 32])
hold on
shadedErrorBar(days,mean_ACC_action,std_ACC_action,'r',0);
ylim([0.5 1])
xlim([10 32])
hold off
