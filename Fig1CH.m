% Fig.1 C and G pool selected mice and calculate behavior

% add folders containing required scripts to the current MATLAB path
% required scripts:
%   - uigetfile_n_dir.m
% path_raw = pwd;

path_raw = 'D:\Documents\data\texture_discrimination\mtp_12\mtp_12_TE';

% define as 1 for expert animals i.e. pool sessions above the performance threshold
% define as 2 for learning animals i.e. pool sessions between 60 and 80% of performance

% use 80 as a performance threshold
performance_threshold = 0.7;

% select mice folders
[paths_mice] = uigetfile_n_dir(path_raw,'Select mice folders');
mice = length(paths_mice);

% get mouse ID
mouse_IDs = cell([1,mice]);
for count_mice = 1:mice
    [~,mouse_IDs{count_mice},~] = fileparts(paths_mice{count_mice});
end

% import all channel labels for comparison
channel_labels_all = cell([1,mice]);
for count_mice = 1:mice
    temp_file_name = fullfile(paths_mice{count_mice},'channel_labels.mat');
    if exist(temp_file_name,'file')
        temp = load(temp_file_name);
        channel_labels_all{count_mice} = temp.channel_labels;
    else
        error('ERROR: channel_labels.mat file missing in folder: %s',paths_mice{count_mice})
    end
end

% preallocate behavioral variables

for count_mice = 1:mice
        
        % Get a list of all files and folders in the mouse folder
        all_files = dir(paths_mice{count_mice});
        % Remove . and .. folders
        all_files(1:2) = [];
        % Get a logical vector that tells which is a directory
        dirFlags = [all_files.isdir];
        % Extract only those that are directories.
        sessions_IDs = {all_files(dirFlags).name};
        sessions = length(sessions_IDs);
        
        mouse_ID= mouse_IDs{count_mice};
        
        % performance
        performance{1,count_mice} = cell(1,sessions);
        % calcium dynamics
        data{1,count_mice} = cell(1,sessions);
        % reaction times
        behavior{1,count_mice} = cell(1,sessions);


end

for count_mice = 1:mice
        
        % Get a list of all files and folders in the mouse folder
        all_files = dir(paths_mice{count_mice});
        % Remove . and .. folders
        all_files(1:2) = [];
        % Get a logical vector that tells which is a directory
        dirFlags = [all_files.isdir];
        % Extract only those that are directories.
        sessions_IDs = {all_files(dirFlags).name};
        sessions = length(sessions_IDs);
        
        mouse_ID= mouse_IDs{count_mice};
        
        % import performance files and create performance vector

        for count_sessions = 1:sessions

            % -----------------------
            % import performance file
            % -----------------------
            
            path_session_id = fullfile(paths_mice{count_mice},sessions_IDs{1,count_sessions});
            temp_file_name = fullfile(path_session_id,'Transfer Entropy','performance.mat');
            if exist(temp_file_name,'file')
                temp_loaded_data = load(temp_file_name);
                performance{1,count_mice}{1,count_sessions} = temp_loaded_data.performance;
            else
                disp(['ERROR: performance.mat file missing in folder: %s', fullfile(path_session_id,'Transfer Entropy')])
            end
            
            % -----------------------
            % import calcium dynamics
            % -----------------------
            
            path_session_id = fullfile(paths_mice{count_mice},sessions_IDs{1,count_sessions});
            temp_file_name = fullfile(path_session_id,'data.mat');
            if exist(temp_file_name,'file')
                temp_loaded_data = load(temp_file_name);
                data{1,count_mice}{1,count_sessions} = temp_loaded_data;
            else
                disp(['ERROR: data.mat file missing in folder: %s',fullfile(path_session_id)])
            end
            
            
            % --------------------
            % import reaction time
            % --------------------
            
            path_reaction_session_id = strrep(path_session_id,'mtp_12_TE','mtp_12_lick_reaction');
            temp_file_name = fullfile(path_reaction_session_id,'RT_264.mat');
            if exist(temp_file_name,'file')==2
                temp_loaded_data = load(temp_file_name);
                behavior{1,count_mice}{1,count_sessions}.reaction_time = temp_loaded_data;
                
            else
%                 error('ERROR: RT_264.mat file missing in folder: %s',fullfile(path_reaction_session_id))
                disp(['ERROR: RT_264.mat file missing in folder: %s',fullfile(path_reaction_session_id)])
            end
            
            % ----------------
            % import lick rate
            % ----------------

            temp_file_name = fullfile(path_reaction_session_id,'lick_traces.mat');
            if exist(temp_file_name,'file')==2
                temp_loaded_data = load(temp_file_name);
                lick_rate_all = temp_loaded_data;
                
            lick_rate= [];
            input_to_decimate=[];
            % Decimate 100Hz to 20Hz 
            max_lick_thres=2.63;
            min_lick_thres=2.61;
            
            % ------ prepare all Hit trials as test trials  ------
            test_trials = sort([behavior{1,count_mice}{1,count_sessions}.reaction_time.go_trials ...
                behavior{1,count_mice}{1,count_sessions}.reaction_time.early_go_trials],'ascend');
            
            input_to_decimate= double(or(lick_rate_all.licks_go>max_lick_thres ...
                ,lick_rate_all.licks_go<min_lick_thres)');
           
            [n,m] = size(input_to_decimate);
            % sample rate 100Hz
            delta = 10; % samples (window of 100 ms)
            rate_to_decimate=[];
            for ii = 1:length(test_trials)
                for j = 1:1:m-delta
                    rate_to_decimate(ii,j)= sum(input_to_decimate(ii,j:j+delta))*(delta/4);
                end
            end
           
            lick_rate=[];
            
            for ii = 1:length(test_trials)
                lick_rate(ii,:) = resample(rate_to_decimate(ii,:),20,100);
            end
            % imagesc(lick_rate)          
            % take first 8 seconds of a trial
            Fs=20;
            trial_time=1:8*Fs-4; % 4 is a time window delta
            behavior{1,count_mice}{1,count_sessions}.lick_rate=lick_rate(:,trial_time);

            % ------ correlate calcium signal to a lick rate------
            % take first 8 seconds of a trial
            trial_time=1:Fs*8;
            input1=data{1, count_mice}{1, count_sessions}.data(test_trials,trial_time,:);
            input2=lick_rate(:,trial_time);
                        
            cc= crosscorr_behavior_channels(input1,input2);
            
%             corr_to_lick_mean(:,count_sessions)= cc.mean_coeff;
%             corr_to_lick_var(:,count_sessions)= cc.var_coeff;
            behavior{1,count_mice}{1,count_sessions}.corr_to_lick=cc;
            
%             subplot(2,1,1)
%             imagesc(input2)
%             subplot(2,1,2)
%             imagesc(input1(:,:,10))

            cc=[];
                
            else
%                 error('ERROR: lick_traces.mat file missing in folder: %s',fullfile(path_reaction_session_id))
                disp(['ERROR: lick_traces.mat file missing in folder: %s',fullfile(path_reaction_session_id)])
                behavior{1,count_mice}{1,count_sessions}.corr_to_lick=[];
                behavior{1,count_mice}{1,count_sessions}.lick_rate=[];
            end
            
            % -----------------------
            % import DLC results
            % -----------------------
            
            path_DLC_id = fullfile(path_raw,'dlc_results');
            
            % for Hit trials
            temp_file_name = [sessions_IDs{1,count_sessions} '.h5'];
            temp_path_name = fullfile(path_DLC_id,temp_file_name);
            
            if exist(temp_path_name,'file')==2
                
                % 1'Texture' 2'Nose' 3'Whisker1_bottom' 4'Whisker1_top'
                % 5'Whisker2_bottom' 6'Whisker2_top' 7'Whisker3_bottom'
                % 8'Whisker3_top' 9'Touch'
                node_names_table=h5read(temp_path_name,'/NODE_NAMES');
                % Dimentions: trials x points x frames
                p_table=h5read(temp_path_name,'/P');
                % create a mask if p<0.05
                p_mask=1-p_table;
                p_mask(p_mask>0.05)=NaN;
                p_mask(p_mask<=0.05)=1;
                
                x_table=h5read(temp_path_name,'/X');
                y_table=h5read(temp_path_name,'/Y');

                x_table=x_table.*p_mask;
                y_table=y_table.*p_mask;
                
                [tracked_trials,tracked_points,tracked_frames]=size(p_table);

                % --- calculate from coordinates to behavioral variables ---
                % define zero = the nose coordinates
                nose_xcoordinates= x_table(:,2,:);
                nose_xcoordinates=0.5 * (fillmissing(nose_xcoordinates, 'previous') ...
                    + fillmissing(nose_xcoordinates, 'next'));
                nose_ycoordinates= y_table(:,2,:);
                nose_ycoordinates= 0.5 * (fillmissing(nose_ycoordinates, 'previous') ...
                    + fillmissing(nose_ycoordinates, 'next'));
                nose_movement= squeeze(sqrt(nose_xcoordinates.^2+nose_ycoordinates.^2));
                
                % angle to the first whisker
                whisker_xbase=x_table(:,3,:);
                whisker_xbase=0.5 * (fillmissing(whisker_xbase, 'previous') ...
                    + fillmissing(whisker_xbase, 'next'));
                
                whisker_ybase=y_table(:,3,:);
                whisker_ybase=0.5 * (fillmissing(whisker_ybase, 'previous') ...
                    + fillmissing(whisker_ybase, 'next'));
                
                whisker_xtop=x_table(:,4,:);
                whisker_xtop=0.5 * (fillmissing(whisker_xtop, 'previous') ...
                    + fillmissing(whisker_xtop, 'next'));
                whisker_ytop=y_table(:,4,:);
                whisker_ytop=0.5 * (fillmissing(whisker_ytop, 'previous') ...
                    + fillmissing(whisker_ytop, 'next'));
                
                x1=whisker_xtop-whisker_xbase;
                y1=whisker_ytop-whisker_ybase;
                
                x0=nose_xcoordinates-whisker_xbase;
                y0=nose_ycoordinates-whisker_ybase;
                
                %v0=[x0 y0];
                %v1=[x1 y1];
                %whisker_angle_1=acosd(squeeze(dot(v1,v0,2)./(sqrt(x1.^2+y1.^2).*sqrt(x0.^2+y0.^2))));
                whisker_angle_1=acosd(squeeze((x1.*x0+y1.*y0)./(sqrt(x1.^2+y1.^2).*sqrt(x0.^2+y0.^2))));
                
                % angle to the second whisker
                whisker_xbase=x_table(:,5,:);
                whisker_xbase=0.5 * (fillmissing(whisker_xbase, 'previous') ...
                    + fillmissing(whisker_xbase, 'next'));
                
                whisker_ybase=y_table(:,5,:);
                whisker_ybase=0.5 * (fillmissing(whisker_ybase, 'previous') ...
                    + fillmissing(whisker_ybase, 'next'));
                
                whisker_xtop=x_table(:,6,:);
                whisker_xtop=0.5 * (fillmissing(whisker_xtop, 'previous') ...
                    + fillmissing(whisker_xtop, 'next'));
                whisker_ytop=y_table(:,6,:);
                whisker_ytop=0.5 * (fillmissing(whisker_ytop, 'previous') ...
                    + fillmissing(whisker_ytop, 'next'));
                
                x1=whisker_xtop-whisker_xbase;
                y1=whisker_ytop-whisker_ybase;
                
                x0=nose_xcoordinates-whisker_xbase;
                y0=nose_ycoordinates-whisker_ybase;
                whisker_angle_2=acosd(squeeze((x1.*x0+y1.*y0)./(sqrt(x1.^2+y1.^2).*sqrt(x0.^2+y0.^2))));

                % texture movement
                %texture_xcoordinates= x_table(:,1,:).*p_mask(:,1,:);
                %texture_ycoordinates= y_table(:,1,:).*p_mask(:,1,:);
                
                %texture_movement= squeeze(sqrt(texture_xcoordinates.^2+texture_ycoordinates.^2));
                
                % touch moments
                touch= squeeze(p_mask(:,9,:));
                touch(isnan(touch))=0;
                
    behavior{1,count_mice}{1,count_sessions}.nose=nose_movement;
    behavior{1,count_mice}{1,count_sessions}.whisker_angle_1=whisker_angle_1;
    behavior{1,count_mice}{1,count_sessions}.whisker_angle_2=whisker_angle_2;
    behavior{1,count_mice}{1,count_sessions}.touch=touch;
                

    % calculate envelope and downsample to 40 Hz
    %input_trials=behaviorvar{1,count_mice}{1,count_sessions};

    % --- specify trial type ---
    input_trial_type=test_trials; % taken from reaction time (from lick above)
    [trials, frames]=size(behavior{1,count_mice}{1,count_sessions}.nose);
    
    % truncate trial number by the max of behavior
    if max(input_trial_type)>trials
        input_trial_type=input_trial_type(1:find(input_trial_type<trials,1,'last'));
    end
    % check for the frame rate
    if frames>1500
        Fs=200;
        resample_to_Fs=20;
        trial_time= 1:8*Fs;

        % nose
        test_input=behavior{1,count_mice}{1,count_sessions}.nose(input_trial_type,trial_time);
        nose=resample(calculate_envelope(test_input)',resample_to_Fs,Fs)';
        [count_i, ~]=find(isnan(nose));
        nose_NaNtrials=unique(count_i);

        % whisker1
        test_input=behavior{1,count_mice}{1,count_sessions}.whisker_angle_1(input_trial_type,trial_time);
        whisker1=resample(calculate_envelope(test_input)',resample_to_Fs,Fs)';
        [count_i, ~]=find(isnan(whisker1));
        whisker1_NaNtrials=unique(count_i);

        % whisker2
        test_input=behavior{1,count_mice}{1,count_sessions}.whisker_angle_2(input_trial_type,trial_time);
        whisker2=resample(calculate_envelope(test_input)',resample_to_Fs,Fs)';
        [count_i, ~]=find(isnan(whisker2));
        whisker2_NaNtrials=unique(count_i);

        %touch=mean(behavior{1,count_mice}{1,count_sessions}.touch(input_trial_type,trial_time),2);
        touch=resample(behavior{1,count_mice}{1,count_sessions}.touch(input_trial_type,trial_time)',resample_to_Fs,Fs)';
        [count_i, ~]=find(isnan(touch));
        touch_NaNtrials=unique(count_i);
        
        trials_to_remove=unique([whisker1_NaNtrials; whisker2_NaNtrials; nose_NaNtrials; touch_NaNtrials]);


    elseif frames<1500
        
        Fs=40;
        resample_to_Fs=20;
        trial_time= 1:8*Fs;
        
        % nose
        test_input=behavior{1,count_mice}{1,count_sessions}.nose(input_trial_type,trial_time);
        nose=resample(calculate_envelope(test_input)',resample_to_Fs,Fs)';
        [count_i, ~]=find(isnan(nose));
        nose_NaNtrials=unique(count_i);

        % whisker1
        test_input=behavior{1,count_mice}{1,count_sessions}.whisker_angle_1(input_trial_type,trial_time);
        whisker1=resample(calculate_envelope(test_input)',resample_to_Fs,Fs)';
        [count_i, ~]=find(isnan(whisker1));
        whisker1_NaNtrials=unique(count_i);
        
        % whisker2
        test_input=behavior{1,count_mice}{1,count_sessions}.whisker_angle_2(input_trial_type,trial_time);
        whisker2=resample(calculate_envelope(test_input)',resample_to_Fs,Fs)';
        [count_i, ~]=find(isnan(whisker2));
        whisker2_NaNtrials=unique(count_i);

        %touch
        touch=resample(behavior{1,count_mice}{1,count_sessions}.touch(input_trial_type,trial_time)',resample_to_Fs,Fs)';
        [count_i, ~]=find(isnan(touch));
        touch_NaNtrials=unique(count_i);
        
        trials_to_remove=unique([whisker1_NaNtrials; whisker2_NaNtrials; nose_NaNtrials; touch_NaNtrials]);
        
    end
    
    Fs=resample_to_Fs;
    
    % --- zscore to first half a second ---
    zscore_time_interval= 1:0.5*Fs;
    % remove NaN trials (usually first and last)
    trials_in=setdiff(1:size(nose,1),trials_to_remove);
    nose=(nose(trials_in,:)-mean(nose(trials_in,zscore_time_interval),2,'omitnan'))./std(nose(trials_in,zscore_time_interval),[],2,'omitnan');
    whisker1=(whisker1(trials_in,:)-mean(whisker1(trials_in,zscore_time_interval),2,'omitnan'))./std(whisker1(trials_in,zscore_time_interval),[],2,'omitnan');
    whisker2=(whisker2(trials_in,:)-mean(whisker2(trials_in,zscore_time_interval),2,'omitnan'))./std(whisker2(trials_in,zscore_time_interval),[],2,'omitnan');
    touch=touch(trials_in,:);
    
    % ------ correlate calcium signal to whisking envelope------
    % take first 8 seconds of a trial
    trial_time=1:Fs*8-2;
    input1=data{1, count_mice}{1, count_sessions}.data(trials_in,trial_time,:);
    % first whisker
    input2=whisker1;
                        
            cc= crosscorr_behavior_channels(input1,input2);
            
%             corr_to_whisk1_mean(:,count_sessions)= cc.mean_coeff;
%             corr_to_whisk1_var(:,count_sessions)= cc.var_coeff;
            behavior{1,count_mice}{1,count_sessions}.corr_to_whisk1=cc;
            behavior{1,count_mice}{1,count_sessions}.whisker_envelope1=whisker1;

            
            cc=[];
    % second whisker
    input2=whisker2;
                        
            cc= crosscorr_behavior_channels(input1,input2);
            
%             corr_to_whisk2_mean(:,count_sessions)= cc.mean_coeff;
%             corr_to_whisk2_var(:,count_sessions)= cc.var_coeff;
            behavior{1,count_mice}{1,count_sessions}.corr_to_whisk2=cc;
            behavior{1,count_mice}{1,count_sessions}.whisker_envelope2=whisker2;
            cc=[];
            
    % touch
    input2=[];
    h = exp(-trial_time/10);
    for count_trials=1:size(touch,1)
        test_conv=conv(touch(count_trials,trial_time),h);
        input2(count_trials,:)=test_conv(1:end-max(trial_time)+1);
    end
%            subplot(2,1,1)
%            plot(touch(count_trials,trial_time))
%            subplot(2,1,2)
%            plot(test_conv(1:end-max(trial_time)))
           
            cc= crosscorr_behavior_channels(input1,input2);
            
%             corr_to_whisk2_mean(:,count_sessions)= cc.mean_coeff;
%             corr_to_whisk2_var(:,count_sessions)= cc.var_coeff;
            behavior{1,count_mice}{1,count_sessions}.corr_to_touch=cc;
            behavior{1,count_mice}{1,count_sessions}.whisker_touch=touch;
            cc=[];
            
                            
   else
        disp(['ERROR: % DLC file missing: ' temp_path_name])
        
        behavior{1,count_mice}{1,count_sessions}.corr_to_whisk1=[];
        behavior{1,count_mice}{1,count_sessions}.corr_to_whisk2=[];
        behavior{1,count_mice}{1,count_sessions}.corr_to_touch=[];


        behavior{1,count_mice}{1,count_sessions}.whisker_envelope1=[];
        behavior{1,count_mice}{1,count_sessions}.whisker_envelope2=[];
        behavior{1,count_mice}{1,count_sessions}.whisker_touch=[];

    end


        end % count_sessions
end % count_mice
    



%% sort by the performance threshold
performance_threshold=0.7;
lick_rate_expert=[];
whisk_envelope1_expert=[];
whisk_envelope2_expert=[];
touch_expert=[];

lick_rate_naive=[];
whisk_envelope1_naive=[];
whisk_envelope2_naive=[];
touch_naive=[];

lick_rate=[];
reaction_time_naive=[];
reaction_time_expert=[];
whisk_envelope1=[];
whisk_envelope2=[];
whisk_touch=[];
mice=size(behavior,2);

for count_mice=1:mice
    
    high_perf_sessions_indices = find(cell2mat(performance{1,count_mice})...
    > performance_threshold);
    low_perf_sessions_indices = find(cell2mat(performance{1,count_mice})...
    <= performance_threshold);

%     if isempty(high_perf_sessions_indices)
%         lick_rate_expert=nan(1,156);
%         whisk_envelope1_expert=nan(1,158);
%         whisk_envelope2_expert=nan(1,158);
%         touch_expert=nan(1,160)
%     end
for count_sessions=1:length(high_perf_sessions_indices)
    
    lick_rate_expert=cat(1,lick_rate_expert, behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.lick_rate);
    % for Hit trials only as lick rate
    hit=behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.reaction_time.go_trials;
    reaction_time_expert= cat(1, reaction_time_expert, behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.reaction_time.reaction_time(hit));

    if isempty(behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.whisker_envelope1)

    else        
        whisk_envelope1_expert=cat(1,whisk_envelope1_expert, ...
            behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.whisker_envelope1);
        whisk_envelope2_expert=cat(1,whisk_envelope2_expert, ...
            behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.whisker_envelope2);
        touch_expert=cat(1,touch_expert, ...
            behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.whisker_touch);
    end
end

%     if isempty(low_perf_sessions_indices)
%         lick_rate_naive=nan(1,156);
%         whisk_envelope1_naive=nan(1,158);
%         whisk_envelope2_naive=nan(1,158);
%         touch_naive=nan(1,160);
%     end
for count_sessions=1:length(low_perf_sessions_indices)
    
    lick_rate_naive=cat(1,lick_rate_naive, behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.lick_rate);
    % for hit trials nly as lick rate
    hit=behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.reaction_time.go_trials;
    reaction_time_naive= cat(1, reaction_time_naive, behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.reaction_time.reaction_time(hit));

        if isempty(behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.whisker_envelope1)

        else        
            whisk_envelope1_naive=cat(1,whisk_envelope1_naive, ...
                behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.whisker_envelope1);
            whisk_envelope2_naive=cat(1,whisk_envelope2_naive, ...
                behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.whisker_envelope2);
            touch_naive=cat(1,touch_naive, ...
                behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.whisker_touch);
        end
end

lick_rate(count_mice).naive=mean(lick_rate_naive,1,'omitnan');
lick_rate(count_mice).expert=mean(lick_rate_expert,1,'omitnan');
lick_rate_naive=[];
lick_rate_expert=[];

% reaction_time_naive(reaction_time_naive<0)=NaN;
% reaction_time_expert(reaction_time_expert<0)=NaN;

reaction_time(count_mice).naive=mean(reaction_time_naive,1,'omitnan');
reaction_time(count_mice).expert=mean(reaction_time_expert,1,'omitnan');
reaction_time_naive=[];
reaction_time_expert=[];

whisk_envelope1(count_mice).naive=mean(whisk_envelope1_naive,1,'omitnan');
whisk_envelope1(count_mice).expert=mean(whisk_envelope1_expert,1,'omitnan');
whisk_envelope1_naive=[];
whisk_envelope1_expert=[];

whisk_envelope2(count_mice).naive=mean(whisk_envelope2_naive,1,'omitnan');
whisk_envelope2(count_mice).expert=mean(whisk_envelope2_expert,1,'omitnan');
whisk_envelope2_naive=[];
whisk_envelope2_expert=[];

% physical constrains on touch- not possible
% touch_naive(:,[1:50, 120:160])=0;
% touch_expert(:,[1:50, 120:160])=0;

whisk_touch(count_mice).naive=mean(touch_naive,1,'omitnan');
whisk_touch(count_mice).expert=mean(touch_expert,1,'omitnan');
touch_naive=[];
touch_expert=[];
 
end

mean_lick_expert=[];
mean_lick_naive=[];
mean_whisk_envelope1_naive=[];
mean_whisk_envelope1_expert=[];
mean_whisk_envelope2_naive=[];
mean_whisk_envelope2_expert=[];
mean_touch_expert=[];
mean_touch_naive=[];
mean_rt_naive=[];
mean_rt_expert=[];

for count_mice=1:mice

% ---lick rate
    if isempty(lick_rate(count_mice).naive)
        mean_lick_naive=[mean_lick_naive; nan(1,156)];
    else
        mean_lick_naive=[mean_lick_naive; lick_rate(count_mice).naive];
    end
    
    if isempty(lick_rate(count_mice).expert)
        mean_lick_expert=[mean_lick_expert; nan(1,156)];
    else
        mean_lick_expert=[mean_lick_expert; lick_rate(count_mice).expert];
    end
    
% ---reaction time
    if isempty(reaction_time(count_mice).naive)
        mean_rt_naive=[mean_rt_naive NaN];
    else
        mean_rt_naive=[mean_rt_naive reaction_time(count_mice).naive];
    end
    
    mean_rt_expert=[mean_rt_expert reaction_time(count_mice).expert];

% ---whisk envelope
    if isempty(whisk_envelope1(count_mice).naive)
        mean_whisk_envelope1_naive=[mean_whisk_envelope1_naive; nan(1,158)];
    else
        mean_whisk_envelope1_naive=[mean_whisk_envelope1_naive; whisk_envelope1(count_mice).naive];
    end
    
    if isempty(whisk_envelope1(count_mice).expert)
        mean_whisk_envelope1_expert=[mean_whisk_envelope1_expert; nan(1,158)];
    else
        mean_whisk_envelope1_expert=[mean_whisk_envelope1_expert; whisk_envelope1(count_mice).expert];
    end
    
    if isempty(whisk_envelope2(count_mice).naive)
        mean_whisk_envelope2_naive=[mean_whisk_envelope2_naive; nan(1,158)];
    else
        mean_whisk_envelope2_naive=[mean_whisk_envelope2_naive; whisk_envelope2(count_mice).naive];
    end
    
    if isempty(whisk_envelope2(count_mice).expert)
        mean_whisk_envelope2_expert=[mean_whisk_envelope2_expert; nan(1,158)];
    else
        mean_whisk_envelope2_expert=[mean_whisk_envelope2_expert; whisk_envelope2(count_mice).expert];
    end

% ---touch
    if isempty(whisk_touch(count_mice).expert)
        mean_touch_expert=[mean_touch_expert; nan(1,160)];
    else
        mean_touch_expert=[mean_touch_expert; whisk_touch(count_mice).expert];
    end
    
    if isempty(whisk_touch(count_mice).naive)
        mean_touch_naive=[mean_touch_naive; nan(1,160)];
    else
        mean_touch_naive=[mean_touch_naive; whisk_touch(count_mice).naive];
    end
end

%% plot behavioral variables 

p_threshold=0.05;
subplot(8,1,1)
t=1/Fs:1/Fs:size(mean_lick_naive,2)/Fs;
shadedErrorBar(t,mean(mean_lick_naive,1,'omitnan'),std(mean_lick_naive,'omitnan')/sqrt(mice),'k');
hold on
shadedErrorBar(t,mean(mean_lick_expert,1,'omitnan'),std(mean_lick_expert,'omitnan')/sqrt(mice),'b');
hold off

subplot(8,1,2)
for count_trial_time=1:size(mean_lick_naive,2)
    p=ranksum(mean_lick_naive(:,count_trial_time),mean_lick_expert(:,count_trial_time));
    p_val(count_trial_time)=p;
end
p_val(p_val>p_threshold)=1;
p_val(p_val<p_threshold)=0;
map=[0, 0, 0;
    1, 1, 1];
imagesc(p_val);
colormap(map)

%xlim([3 4])
subplot(8,1,3)
t=1/Fs:1/Fs:size(mean_whisk_envelope1_naive,2)/Fs;
shadedErrorBar(t,mean(mean_whisk_envelope1_naive,1,'omitnan'),std(mean_whisk_envelope1_naive,'omitnan')/sqrt(mice),'k');
hold on
shadedErrorBar(t,mean(mean_whisk_envelope1_expert,1,'omitnan'),std(mean_whisk_envelope1_expert,'omitnan')/sqrt(mice),'b');
hold off

subplot(8,1,4)
for count_trial_time=1:size(mean_whisk_envelope1_naive,2)
    p=ranksum(mean_whisk_envelope1_naive(:,count_trial_time),mean_whisk_envelope1_expert(:,count_trial_time));
    p_val(count_trial_time)=p;
end
p_val(p_val>p_threshold)=1;
p_val(p_val<p_threshold)=0;

imagesc(p_val);
colormap(map)

subplot(8,1,5)
t=1/Fs:1/Fs:size(mean_whisk_envelope2_naive,2)/Fs;
shadedErrorBar(t,mean(mean_whisk_envelope2_naive,1,'omitnan'),std(mean_whisk_envelope2_naive,'omitnan')/sqrt(mice),'k');
hold on
shadedErrorBar(t,mean(mean_whisk_envelope2_expert,1,'omitnan'),std(mean_whisk_envelope2_expert,'omitnan')/sqrt(mice),'b');
hold off

subplot(8,1,6)
for count_trial_time=1:size(mean_whisk_envelope2_naive,2)
    p=ranksum(mean_whisk_envelope2_naive(:,count_trial_time),mean_whisk_envelope2_expert(:,count_trial_time));
    p_val(count_trial_time)=p;
end
p_val(p_val>p_threshold)=1;
p_val(p_val<p_threshold)=0;

imagesc(p_val);
colormap(map)

subplot(8,1,7)
t=1/Fs:1/Fs:size(mean_touch_naive,2)/Fs;
shadedErrorBar(t,mean(mean_touch_naive,1,'omitnan'),std(mean_touch_naive,'omitnan')/sqrt(mice),'k');
hold on
shadedErrorBar(t,mean(mean_touch_expert,1,'omitnan'),std(mean_touch_expert,'omitnan')/sqrt(mice),'b');
hold off

subplot(8,1,8)
for count_trial_time=1:size(mean_touch_naive,2)
    p=ranksum(mean_touch_naive(:,count_trial_time),mean_touch_expert(:,count_trial_time));
    p_val(count_trial_time)=p;
end
p_val(p_val>p_threshold)=1;
p_val(p_val<p_threshold)=0;

imagesc(p_val);
colormap(map)
        
%% calculate mean envelope angle in the time window [3 3.5]
% size mice x trial time 0-8 sec , 20Hz
texture=(60:70);
m_inp=mean(mean(mean_whisk_envelope2_naive(:,texture),2),1,'omitnan');
sd_inp=std(mean(mean_whisk_envelope2_naive(:,texture),2),[],1,'omitnan');
disp(['naive ' num2str(m_inp) '+/-' num2str(sd_inp)])

m_inp=mean(mean(mean_whisk_envelope2_expert(:,texture),2),1,'omitnan');
sd_inp=std(mean(mean_whisk_envelope2_expert(:,texture),2),[],1,'omitnan');
disp(['expert ' num2str(m_inp) '+/-' num2str(sd_inp)])

p=ranksum(mean(mean_whisk_envelope2_naive(:,texture),2),mean(mean_whisk_envelope2_expert(:,texture),2))
disp(['Mann Whitney U test p value ' num2str(p)])

% calculate change in a reaction time
m_inp=mean(mean_rt_naive,2,'omitnan');
sd_inp=std(mean_rt_naive,[],2,'omitnan');
disp(['naive ' num2str(m_inp) '+/-' num2str(sd_inp)])

m_inp=mean(mean_rt_expert,2);
sd_inp=std(mean_rt_expert,[],2);
disp(['expert ' num2str(m_inp) '+/-' num2str(sd_inp)])

p=ranksum(mean_rt_naive',mean_rt_expert')
disp(['Mann Whitney U test p value ' num2str(p)])

%% plot touch

% define time window to plot the touch
texture_approach=50:70;
sz=20;
% define time window to check correlations
% usual texture time window is [3 3,5]
% here it is nesessary to check for the dependancy on early touch
texture=60:70;
c=1:mice;
subplot(4,2,1)
imagesc(mean_touch_naive(:,texture_approach), [0 1])
title('touch naive')
subplot(4,2,2)
imagesc(mean_touch_expert(1:end,texture_approach), [0 1])
title('touch expert')

subplot(4,2,3)
scatter(mean(mean_touch_naive(:,texture),2,'omitnan'), mean(mean_whisk_envelope2_naive(:,texture),2,'omitnan'),sz,c,'filled');
% xlim([0.5 1]);
% ylim([0 40]);
legend
lsline
title('touch and whisk naive')

subplot(4,2,4)
scatter(mean(mean_touch_expert(:,texture),2), mean(mean_whisk_envelope2_expert(:,texture),2,'omitnan'),sz,c,'filled');
lsline
% xlim([0.5 1]);
% ylim([0 40]);
title('touch and whisk expert')

subplot(4,2,5)
scatter(mean(mean_touch_naive(:,texture),2), mean_rt_naive',sz,c,'filled');
lsline

y =  mean_rt_naive';
x = mean(mean_touch_naive(:,texture),2);
nan_ind=find(isnan(x));
x(isnan(x))=[];
y(nan_ind)=[];
[p,S] = polyfit(x,y,1);
% Evaluate the first-degree polynomial fit in p at the points in x.
[y_fit,delta] = polyval(p,x,S);
yresid = y - y_fit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
% xlim([0.5 1]);
% ylim([0 3]);
title(['touch and reaction time naive r sq' num2str(rsq)])

subplot(4,2,6)
scatter(mean(mean_touch_expert(:,texture),2), mean_rt_expert',sz,c,'filled');
% xlim([0.5 1]);
% ylim([0 3]);
lsline
y =  mean_rt_expert';
x = mean(mean_touch_expert(:,texture),2);
nan_ind=find(isnan(x));
x(isnan(x))=[];
y(nan_ind)=[];
[p,S] = polyfit(x,y,1);
% Evaluate the first-degree polynomial fit in p at the points in x.
[y_fit,delta] = polyval(p,x,S);
yresid = y - y_fit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal
% xlim([0.5 1]);
% ylim([0 3]);
title(['touch and reaction time expert r sq' num2str(rsq)])

subplot(4,2,7)
scatter(mean(mean_whisk_envelope2_naive(:,texture),2,'omitnan'), mean_rt_naive',sz,c,'filled');
% xlim([0 40]);
% ylim([0 3]);
lsline
title('whisk and reaction time naive')

subplot(4,2,8)
scatter(mean(mean_whisk_envelope2_expert(:,texture),2,'omitnan'), mean_rt_expert',sz,c,'filled');
xlim([0 40]);
ylim([0 3]);
lsline
title('whisk and reaction time expert')

%% plot behavior acorss mice naive vs expert

texture=50:60;

subplot(3,1,1)
scatter(mean(mean_touch_naive(:,texture),2,'omitnan'), mean(mean_touch_expert(:,texture),2,'omitnan'),sz,c,'filled');
% xlim([0.5 1]);
% ylim([0 40]);
legend
lsline
title('touch ')

subplot(3,1,2)
scatter(mean(mean_whisk_envelope2_naive(:,texture),2,'omitnan'), mean(mean_whisk_envelope2_expert(:,texture),2,'omitnan'),sz,c,'filled');
lsline
% xlim([0.5 1]);
% ylim([0 40]);
title('whisk')

subplot(3,1,3)
scatter(mean_rt_naive', mean_rt_expert',sz,c,'filled');
lsline
% xlim([0.5 1]);
% ylim([0 3]);
title('reaction time ')

%% sort correlations 
    
label=[];
channel_label = {'GP', 'Cpu', 'VL','VM', 'LDVL', 'Rt','M1','S1_bf','DG','CA1_Py','CA1_Mol','dCA1'};

comb=length(channel_label);


for count_mice=1:mice

        high_perf_sessions_indices = find(cell2mat(performance{1,count_mice})...
        > performance_threshold);
        low_perf_sessions_indices = find(cell2mat(performance{1,count_mice})...
        <= performance_threshold);
    
    for count_channel_labels=1:comb
                    
    label = strncmp(channel_labels_all{count_mice},channel_label{1,count_channel_labels},5);
    channel_sub_network = find(label>0);
    if length(channel_sub_network)>1
        channel_sub_network=channel_sub_network(1,1);
    end
    
    if isempty(channel_sub_network)  
        %channel_sub_network=NaN;
    else

    cc_lick_rate_expert=[];
    cc_whisk_envelope1_expert=[];
    cc_whisk_envelope2_expert=[];
    cc_touch_expert=[];

    cc_lick_rate_naive=[];
    cc_whisk_envelope1_naive=[];
    cc_whisk_envelope2_naive=[];
    cc_touch_naive=[];


    for count_sessions=1:length(high_perf_sessions_indices)
        
        cc_lick_rate_expert=cat(1,cc_lick_rate_expert, behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.corr_to_lick(1,channel_sub_network));
        
        if isempty(behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.corr_to_whisk1)
            cc_whisk_envelope1_expert=NaN;
        else
            cc_whisk_envelope1_expert=cat(1,cc_whisk_envelope1_expert, ...
                behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.corr_to_whisk1(1,channel_sub_network));
        end
        if isempty(behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.corr_to_whisk2)
             cc_whisk_envelope2_expert=NaN;
        else
            cc_whisk_envelope2_expert=cat(1,cc_whisk_envelope2_expert, ...
                behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.corr_to_whisk2(1,channel_sub_network));
        end
        
        if isempty(behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.corr_to_touch)
             cc_touch_expert=NaN;
        else
            cc_touch_expert=cat(1,cc_touch_expert, ...
                behavior{1,count_mice}{1,high_perf_sessions_indices(count_sessions)}.corr_to_touch(1,channel_sub_network));
        end
        
    end
    
    if isempty(low_perf_sessions_indices)
            cc_lick_rate_naive=NaN;
            cc_whisk_envelope1_naive=NaN;
            cc_whisk_envelope2_naive=NaN;
            cc_touch_naive=NaN;
    else

    for count_sessions=1:length(low_perf_sessions_indices)
        
        cc_lick_rate_naive=cat(1,cc_lick_rate_naive, behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.corr_to_lick(1,channel_sub_network));
        if isempty(behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.corr_to_whisk1)
            cc_whisk_envelope1_naive=NaN;
        else
            cc_whisk_envelope1_naive=cat(1,cc_whisk_envelope1_naive, ...
                behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.corr_to_whisk2(1,channel_sub_network));
        end
        if isempty(behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.corr_to_whisk2)
             cc_whisk_envelope2_naive=NaN;
        else
            cc_whisk_envelope2_naive=cat(1,cc_whisk_envelope2_naive, ...
                behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.corr_to_whisk2(1,channel_sub_network));
        end
        
        if isempty(behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.corr_to_touch)
             cc_touch_naive=NaN;
        else
            cc_touch_naive=cat(1,cc_touch_naive, ...
                behavior{1,count_mice}{1,low_perf_sessions_indices(count_sessions)}.corr_to_touch(1,channel_sub_network));
        end
        
        end
    end

    cc_lr_naive(count_channel_labels,count_mice)=mean(cc_lick_rate_naive,1,'omitnan');
    cc_lr_expert(count_channel_labels,count_mice)=mean(cc_lick_rate_expert,1,'omitnan');

    cc_we1_naive(count_channel_labels,count_mice)=mean(cc_whisk_envelope1_naive,1,'omitnan');
    cc_we1_expert(count_channel_labels,count_mice)=mean(cc_whisk_envelope1_expert,1,'omitnan');

    cc_we2_naive(count_channel_labels,count_mice)=mean(cc_whisk_envelope2_naive,1,'omitnan');
    cc_we2_expert(count_channel_labels,count_mice)=mean(cc_whisk_envelope2_expert,1,'omitnan');
    
    cc_tm_naive(count_channel_labels,count_mice)=mean(cc_touch_naive,1,'omitnan');
    cc_tm_expert(count_channel_labels,count_mice)=mean(cc_touch_expert,1,'omitnan');
    
    end % isempty channel sub network

    end % count_channel labels
    
end % count_mice


%% --- boxplot plus lines ---

input_field_1=cc_tm_naive;
input_field_2=cc_tm_expert;

% [p,tbl,stats]=anova1([input_field_1(1,:); input_field_2(1,:)]');
% multcompare(stats)
p_val=[];
for i=1:12
    [h,p] =ttest2(input_field_1(i,:), input_field_2(i,:));
    p_val(i)=p<0.05
end

input_1= std(input_field_1,0,2,'omitnan')./sqrt(mice);
input_2= mean(input_field_1,2,'omitnan');

input_3= std(input_field_2,0,2,'omitnan')./sqrt(mice);
input_4= mean(input_field_2,2,'omitnan');

% [p,tbl,stats]=anova1([input_field_1 input_field_2]');
% multcompare(stats)

h1=barwitherr([input_1 input_3],[input_2 input_4]);
set(gca,'XTickLabel',channel_label)
set(h1(1),'FaceColor',[0.0 0.0 0.6]);
%ylim([0 1]);
hold on

for i=1:12

plot([i-0.1:0.2:0.1+i],[input_field_1(i,:); input_field_2(i,:)],'-o',...
    'LineWidth',1,'Color',[0.8,0.8,0.8], 'MarkerSize',5,...
    'MarkerEdgeColor','b','MarkerFaceColor',[1,1,1]);

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

%% plot scatter plot

% Touch 
input_field_1=cc_tm_naive;
input_field_2=cc_tm_expert;

subplot(1,4,1)

p_val=[];
for count_channel=1:12
    [h,p] =ttest2(input_field_1(count_channel,:), input_field_2(count_channel,:));
    p_val(count_channel)=p;
end

input_1= std(input_field_1,0,2,'omitnan')./sqrt(mice);
input_2= mean(input_field_1,2,'omitnan');

input_3= std(input_field_2,0,2,'omitnan')./sqrt(mice);
input_4= mean(input_field_2,2,'omitnan');

tick_step=1;
all_channels=1:12;
scale_size=300;
% sz = [scale_size*(p_val+0.01) scale_size*(p_val+0.01)];
% c = [input_2' input_4'];
sz = [scale_size*(input_2') scale_size*(input_4')];
c = [p_val<0.05 p_val<0.05];
regions=[(1:12) (1:12)];

scatter([ones(1,12) 2*ones(1,12)],regions,sz,c,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.7 .7 .7],...
              'LineWidth',0.5)
xlim([0 3])     

ay = gca;
ay.YTickLabel = channel_label;
ay.YTick = 1:12;
title('Correlation to Touch')


% --- Whisk 1 ---
input_field_1=cc_we1_naive;
input_field_2=cc_we1_expert;

subplot(1,4,2)

p_val=[];
for count_channel=1:12
    [h,p] =ttest2(input_field_1(count_channel,:), input_field_2(count_channel,:));
    p_val(count_channel)=p;
end

input_1= std(input_field_1,0,2,'omitnan')./sqrt(mice);
input_2= mean(input_field_1,2,'omitnan');

input_3= std(input_field_2,0,2,'omitnan')./sqrt(mice);
input_4= mean(input_field_2,2,'omitnan');

tick_step=1;
all_channels=1:12;
scale_size=300;
% sz = [scale_size*(p_val+0.01) scale_size*(p_val+0.01)];
% c = [input_2' input_4'];
sz = [scale_size*(input_2') scale_size*(input_4')];
c = [p_val<0.05 p_val<0.05];
regions=[(1:12) (1:12)];

scatter([ones(1,12) 2*ones(1,12)],regions,sz,c,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.7 .7 .7],...
              'LineWidth',0.5)
xlim([0 3])     

ay = gca;
ay.YTickLabel = channel_label;
ay.YTick = 1:12;
title('Correlation to Whisk Envelope1')

% --- Whisk 2 ---
input_field_1=cc_we2_naive;
input_field_2=cc_we2_expert;

subplot(1,4,3)

p_val=[];
for count_channel=1:12
    [h,p] =ttest2(input_field_1(count_channel,:), input_field_2(count_channel,:));
    p_val(count_channel)=p;
end

input_1= std(input_field_1,0,2,'omitnan')./sqrt(mice);
input_2= mean(input_field_1,2,'omitnan');

input_3= std(input_field_2,0,2,'omitnan')./sqrt(mice);
input_4= mean(input_field_2,2,'omitnan');

tick_step=1;
all_channels=1:12;
scale_size=300;
% sz = [scale_size*(p_val+0.01) scale_size*(p_val+0.01)];
% c = [input_2' input_4'];
sz = [scale_size*(input_2') scale_size*(input_4')];
c = [p_val<0.05 p_val<0.05];
regions=[(1:12) (1:12)];

scatter([ones(1,12) 2*ones(1,12)],regions,sz,c,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.7 .7 .7],...
              'LineWidth',0.5)
xlim([0 3])     

ay = gca;
ay.YTickLabel = channel_label;
ay.YTick = 1:12;
title('Correlation to Whisk Envelope2')

% --- Lick Rate ---
input_field_1=cc_lr_naive;
input_field_2=cc_lr_expert;

subplot(1,4,4)

p_val=[];
for count_channel=1:12
    [h,p] =ttest2(input_field_1(count_channel,:), input_field_2(count_channel,:));
    p_val(count_channel)=p;
end

input_1= std(input_field_1,0,2,'omitnan')./sqrt(mice);
input_2= mean(input_field_1,2,'omitnan');

input_3= std(input_field_2,0,2,'omitnan')./sqrt(mice);
input_4= mean(input_field_2,2,'omitnan');

tick_step=1;
all_channels=1:12;
scale_size=300;
% sz = [scale_size*(p_val+0.01) scale_size*(p_val+0.01)];
% c = [input_2' input_4'];
sz = [scale_size*(input_2') scale_size*(input_4')];
c = [p_val<0.05 p_val<0.05];
regions=[(1:12) (1:12)];

scatter([ones(1,12) 2*ones(1,12)],regions,sz,c,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.7 .7 .7],...
              'LineWidth',0.5)
xlim([0 3])     

ay = gca;
ay.YTickLabel = channel_label;
ay.YTick = 1:12;
title('Correlation to Lick Rate')

% FileName=['lick_rate_pooled_',datestr(now, 'dd-mmm-yyyy'),'.mat'];
% save(FileName,...
%         'mouse_IDs',...
%         'sessions_IDs',...
%         'channel_labels_all',...
%         'performance',...
%         'behavior','-v7.3');
    