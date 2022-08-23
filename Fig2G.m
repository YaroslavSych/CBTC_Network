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
            [CR1,i_cr1,i_fa1] = intersect(CR_trials,FA_trials+1);
            [CR2,i_cr2,i_fa2] = intersect(CR_trials,FA_trials+2);
            CR_trials(intersect(i_cr1,i_cr2))=[];
            
            
            test_Hit = data{1,count_mice}{1,count_sessions}.data(Hit_trials,:,:);
            test_CR = data{1,count_mice}{1,count_sessions}.data(CR_trials,:,:);
            test_Miss = data{1,count_mice}{1,count_sessions}.data(Miss_trials,:,:);
            test_FA = data{1,count_mice}{1,count_sessions}.data(FA_trials,:,:);
            
            stimulus= 60:70; % from 2.5 to 3.5 seconds of a trial time
            %action= 100:130; % from 5 to 6.5 seconds of a trial time
            
            % --- Hit vs CR discrimination ---
            % use all 12 channels as predictors
            count_predictors=size(test_Hit,3);
            limit_by_trials=40;
            ACC1=[];
            ACC1_var=[];  
            
                [ACC,ACC_var,f1]= calculate_ACC_all_channels(test_Hit,test_CR,stimulus,limit_by_trials);
                ACC1=[ACC1 ACC];
                ACC1_var=[ACC1_var ACC_var];
            
            ACC_stimulus{1,count_mice}{1,count_sessions}.mean= ACC1;
            ACC_stimulus{1,count_mice}{1,count_sessions}.var= ACC1_var;
            ACC1=[];
            ACC1_var=[];
            
            % --- Hit vs FA discrimination ---

                 [ACC,ACC_var,f1]= calculate_ACC_all_channels(test_Hit,test_FA,stimulus,limit_by_trials);
                 ACC1=[ACC1 ACC];
                 ACC1_var=[ACC1_var ACC_var];
            
            ACC_choice{1,count_mice}{1,count_sessions}.mean= ACC1;
            ACC_choice{1,count_mice}{1,count_sessions}.var= ACC1_var;
            ACC1=[];
            ACC1_var=[];
            
            % --- CR vs FA discrimination ---
            
                [ACC,ACC_var,f1]= calculate_ACC_all_channels(test_FA,test_CR,stimulus,limit_by_trials);
                ACC1=[ACC1 ACC];
                ACC1_var=[ACC1_var ACC_var];

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

channel_label = {'all'};

ch=length(channel_label);
label_ind=cell(1,ch);

% preallocate stat_label/learning
for count_mice=1:mice

        pool_ACC_choice{1,count_mice}=cell(1,ch);
        pool_ACC_stimulus{1,count_mice}=cell(1,ch);
        pool_ACC_action{1,count_mice}=cell(1,ch);
        
end

test=[];

    
for count_mice=1:mice
        
        sessions = length(sessions_IDs{1,count_mice});
        for count_sessions=1:sessions
            
            temp_data=ACC_choice{1,count_mice}{1,count_sessions}.mean;
            temp_ACC_choice(count_sessions)=temp_data;
            
            temp_data=ACC_stimulus{1,count_mice}{1,count_sessions}.mean;
            temp_ACC_stimulus(count_sessions)=temp_data;
            
            temp_data=ACC_action{1,count_mice}{1,count_sessions}.mean;
            temp_ACC_action(count_sessions)=temp_data;
            
        end
            pool_ACC_choice{1,count_mice}=temp_ACC_choice;
            pool_ACC_stimulus{1,count_mice}=temp_ACC_stimulus;
            pool_ACC_action{1,count_mice}=temp_ACC_action;
            
            temp_ACC_choice=[];
            temp_ACC_stimulus=[];
            temp_ACC_action=[];
        end
    


% align ACC to learning session

% align all to learning session (first above 75% performance) 
% temp_loaded_data.performance

% ------ loop across mice  ------
learning_ind=[];
%pool_mice=[5 6 7 9 10 11 12 13 14];
pool_mice=[1 2 3 4 5 6 7 8 9 10 11 12 13 14];


for i=1:length(pool_mice)
    
    count_mice=pool_mice(i);
    
    sessions = length(sessions_IDs{1,count_mice});

    learning_ind(i) =find(double(cell2mat(performance{1,count_mice})> performance_threshold), 1, 'first');
    
    % create a vector +/-20 sessions around learning
    length_to_pad=20;
    
    
        pool_ACC_choice{1,count_mice}= align_vectors_to_learning(pool_ACC_choice{1,count_mice},...
        learning_ind(i), length_to_pad);
    
        pool_ACC_stimulus{1,count_mice}= align_vectors_to_learning(pool_ACC_stimulus{1,count_mice},...
        learning_ind(i), length_to_pad);
    
        pool_ACC_action{1,count_mice}= align_vectors_to_learning(pool_ACC_action{1,count_mice},...
        learning_ind(i), length_to_pad);
    
    end


% do a matrix plot for all regions

channel_label = {'all'};

mean_ACC_choice=[];
std_ACC_choice=[];
mean_ACC_stimulus=[];
std_ACC_stimulus=[];
mean_ACC_action=[];
std_ACC_action=[];

    
    all_ACC_choice= [];
    all_ACC_stimulus= [];
    all_ACC_action= [];

    for i=1:length(pool_mice)
        
        count_mice=pool_mice(i);

        all_ACC_choice= [all_ACC_choice pool_ACC_choice{1,count_mice}];
        all_ACC_stimulus= [all_ACC_stimulus pool_ACC_stimulus{1,count_mice}];
        all_ACC_action= [all_ACC_action pool_ACC_action{1,count_mice}];

    end

% all_ACC_choice(all_ACC_choice == 0) = NaN;
% all_ACC_stimulus(all_ACC_stimulus == 0) = NaN;
% all_ACC_action(all_ACC_action == 0) = NaN;

mean_ACC_choice=mean(all_ACC_choice,2,'omitnan')';
std_ACC_choice=std(all_ACC_choice,[],2,'omitnan')'./sqrt(length(pool_mice));

mean_ACC_stimulus=mean(all_ACC_stimulus,2,'omitnan')';
std_ACC_stimulus=std(all_ACC_stimulus,[],2,'omitnan')'./sqrt(length(pool_mice));

mean_ACC_action=mean(all_ACC_action,2,'omitnan')';
std_ACC_action=std(all_ACC_action,[],2,'omitnan')'./sqrt(length(pool_mice));

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

test_channels=repmat([1:length(channel_label)]',1,length(days));

input_days=1:41;
start_day=10;
end_day=30;
scale_size=100;
shift_input=0.5;
scale_input=0.4;
regions=2;

subplot(3,1,1)
input_to_plot=mean_ACC_choice;
input_std_to_plot=std_ACC_choice;
f1=[];

for every_channel=1:length(channel_label)
    
    sz = scale_size*(input_to_plot(every_channel,:)-scale_input);
    delta= 1*(input_to_plot(every_channel,:)-shift_input);
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
xlim([0 regions])
    ax = gca;
    all_channels=1:12;
    ax.XTick = sort([ all_channels ,all_channels+0.2, all_channels+0.4],'ascend');
    title('Choice')
hold off

subplot(3,1,2)
input_to_plot=mean_ACC_stimulus;
input_std_to_plot=std_ACC_stimulus;
f1=[];

for every_channel=1:length(channel_label)
    
    sz = scale_size*(input_to_plot(every_channel,:)-scale_input);
    delta= 1*(input_to_plot(every_channel,:)-shift_input);
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
xlim([0 regions])
    ax = gca;
    all_channels=1:12;
    ax.XTick = sort([ all_channels ,all_channels+0.2, all_channels+0.4],'ascend');
%     ax.XTickLabel = channel_label;
    title('Stimulus')
    
hold off

subplot(3,1,3)
input_to_plot=mean_ACC_action;
input_std_to_plot=std_ACC_action;
f1=[];


for every_channel=1:length(channel_label)
    
    sz = scale_size*(input_to_plot(every_channel,:)-scale_input);
    delta= 1*(input_to_plot(every_channel,:)-shift_input);
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
xlim([0 regions])
    ax = gca;
    all_channels=1:12;
    ax.XTick = sort([ all_channels ,all_channels+0.2, all_channels+0.4],'ascend');
%     tick_labels = string(ax.XTickLabel);
%     tick_labels(1:2:end) = nan; % remove every other one
%     ax.XTickLabel = channel_label;
    title('Action')
    
hold off

%% Ranksum test all_ACC_stimulus and action for significance

    
    all_ACC_choice= [];
    all_ACC_stimulus= [];
    all_ACC_action= [];

    for i=1:length(pool_mice)
        
        count_mice=pool_mice(i);

        all_ACC_choice= [all_ACC_choice pool_ACC_choice{1,count_mice}];
        all_ACC_stimulus= [all_ACC_stimulus pool_ACC_stimulus{1,count_mice}];
        all_ACC_action= [all_ACC_action pool_ACC_action{1,count_mice}];

    end
    % test naive vs expert
    input_1= mean(all_ACC_stimulus(1:20,:),1,'omitnan');
    input_2= mean(all_ACC_stimulus(21:end,:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.stimulus=p;
    
    input_1= mean(all_ACC_action(1:20,:),1,'omitnan');
    input_2= mean(all_ACC_action(21:end,:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.action=p;
    
    input_1= mean(all_ACC_choice(1:20,:),1,'omitnan');
    input_2= mean(all_ACC_choice(21:end,:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.choice=p;
    
    % test learning +/-3 day vs late expert (3 days from each side) 
    input_1= mean(all_ACC_stimulus(17:23,:),1,'omitnan');
    input_2= mean(all_ACC_stimulus([10:12, 28:30],:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.learning_stimulus=p;
    
    input_1= mean(all_ACC_action(18:22,:),1,'omitnan');
    input_2= mean(all_ACC_action([10:12, 28:30],:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.learning_action=p;
    
    input_1= mean(all_ACC_choice(18:22,:),1,'omitnan');
    input_2= mean(all_ACC_choice([10:12, 28:30],:),1,'omitnan');
    [p,h] = ranksum(input_1,input_2);
    p_val.learning_choice=p;
    

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
