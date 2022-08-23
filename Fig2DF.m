count_mice=7;
count_sessions=13;

for count_sessions=[6 17]
% ------ read all trials  ------
Hit_trials = behaviorvar{1, count_mice}{1, count_sessions}.iGO;
CR_trials = behaviorvar{1, count_mice}{1, count_sessions}.iNOGO;
            
Miss_trials = behaviorvar{1, count_mice}{1, count_sessions}.iMISS;
FA_trials = behaviorvar{1, count_mice}{1, count_sessions}.iFA;

test_Hit = data{1,count_mice}{1,count_sessions}.data(Hit_trials,:,:);
test_CR = data{1,count_mice}{1,count_sessions}.data(CR_trials,:,:);
test_Miss = data{1,count_mice}{1,count_sessions}.data(Miss_trials,:,:);
test_FA = data{1,count_mice}{1,count_sessions}.data(FA_trials,:,:);

input_data.Hit=test_Hit;
input_data.CR=test_CR;

% find min sample to decode
end_sample=min([size(input_data.Hit,2) size(input_data.CR,2)]);

number_of_regions=min([size(input_data.Hit,3) size(input_data.CR,3)]);

% --- Predict Strategy with Linear/Quadratic Discriminant Analysis ---

%clear Y cp X predictions Accuracy Accuracy_shuffle f1_score_shuffle

% compute the confusion matrix using stratified 10-fold cross validation
for count_time_interval = 1:end_sample-40

    input_1=[];
    input_2=[];
    
    for count_source=1:number_of_regions
   
        % create predictor
        input_1=input_data.Hit(:,count_time_interval,count_source);
        input_2=input_data.CR(:,count_time_interval,count_source);
        
        input_1(isnan(input_1))=0;
        input_2(isnan(input_2))=0;

        predictor= cat(1,input_1,input_2);
        %sessions=size(predictor,1);

        label=[ones(1,size(input_1,1)) zeros(1,size(input_2,1))];
        Y = label';
        cp = classperf(Y);
        X = predictor;

        % split dataset into train and test and run 5 times 
        number_of_estimations=5;

        for count_runs=1:number_of_estimations

            [train,test] = crossvalind('LeaveMOut',Y,200);
            %mdl = fitcdiscr(X,Y,'DiscrimType','pseudolinear');
            mdl = fitcdiscr(X(train),Y(train),'DiscrimType','pseudolinear');
            predictions = predict(mdl,X(test));
            estimated_performance=classperf(cp,predictions,test);

            % compute the confusion matrix using stratified 10-fold cross validation
            [C,order] = confusionmat(Y(test),predictions);

            Accuracy(count_time_interval,count_source,count_runs)=sum(diag(C))/sum(sum(C));
            f1_score(count_time_interval,count_source,count_runs)=2*C(1,1)/(2*C(1,1)+C(1,2)+C(2,1));

            % do the same for shuffled labels
            Y_shuffle= label(randperm(length(Y)))';
            mdl = fitcdiscr(X(train),Y_shuffle(train),'DiscrimType','pseudolinear');
            predictions = predict(mdl,X(test));
            %estimated_performance=classperf(cp,predictions,test);

            [C,order] = confusionmat(Y_shuffle(test),predictions);            
            Accuracy_shuffle(count_time_interval,count_source,count_runs)=sum(diag(C))/sum(sum(C));
            f1_score_shuffle(count_time_interval,count_source,count_runs)=2*C(1,1)/(2*C(1,1)+C(1,2)+C(2,1));

        end
       
    end

end


% plot accuracy

Fs=20; % sampling rate
time= 0:1/Fs:(end_sample-1)/Fs;

for source=1:12

subplot(12,1,source)

input_to_plot=[];
input_to_plot.mean=mean(Accuracy(:,source,:),3,'omitnan');
input_to_plot.std=std(Accuracy(:,source,:),[],3,'omitnan')./sqrt(number_of_estimations);

if count_sessions<10
     color='m';
else
    color='b';
end

%shadedErrorBar(time,input_to_plot.mean,input_to_plot.std,color,0)
plot(time,input_to_plot.mean,color)


hold on

input_to_plot=[];
input_to_plot.mean=mean(Accuracy_shuffle(:,source,:),3,'omitnan');
input_to_plot.std=std(Accuracy_shuffle(:,source,:),[],3,'omitnan')./sqrt(number_of_estimations);
% input with shuffled labels
if count_sessions<10
     color='w';
else
    color='k';
end

%shadedErrorBar(time,input_to_plot.mean,input_to_plot.std,color,0);
plot(time,input_to_plot.mean,color)
ylim([0.4 0.9])
xlim([0 8])
xticks(1:8)
%title(channel_labels_all{1, count_mice})

end

end % count 1 naive and 1 expert sessions
%% predict from groups of regions

%
joint_sources=cell(1,5);
% 1 motor regions
joint_sources{1,1}=26:28;
% 2 somatosensory regions
joint_sources{1,2}=18:25;
% 3 visual regions
joint_sources{1,3}=1:4;
% 4 auditory regions
joint_sources{1,4}=14:17;
% 5 all regions
joint_sources{1,5}=1:length(ROIs_names);

number_of_groups=size(joint_sources,2);

% compute the confusion matrix using stratified 10-fold cross validation
for count_time_interval = 1:end_sample

    input_1=[];
    input_2=[];
    label=[];
    
    for count_source=1:number_of_groups
   
        % create predictor
        input_1=input_data.passive.tdt(:,count_time_interval,joint_sources{1,count_source});
        input_2=input_data.active.tdt(:,count_time_interval,joint_sources{1,count_source});
        
%         input_1(isnan(input_1))=0;
%         input_2(isnan(input_2))=0;
        
        predictor= cat(1,squeeze(input_1),squeeze(input_2));
        sessions=size(predictor,1);

        label=[ones(1,size(input_1,1)) zeros(1,size(input_2,1))];
        Y = label';
        cp = classperf(Y);
        X = predictor;

    %     % to reduce dimentionality of the predictor you may use PCA transformed predictor
    %     [PCALoadings,PCAScores,PCAVar] = pca(predictor);
    %    
    %     Y = label';
    %     cp = classperf(Y);
    %     % and use 3 first components
    %     X = PCAScores(:,1:3);

        % split dataset into train and test and run 5 times 
        number_of_estimations=10;

        for count_runs=1:number_of_estimations
            
            if sessions<700
                [train,test] = crossvalind('LeaveMOut',Y,150);
            else
                [train,test] = crossvalind('LeaveMOut',Y,200);
            end
            
            %mdl = fitcdiscr(X,Y,'DiscrimType','pseudolinear');
            mdl = fitcdiscr(X(train,:),Y(train),'DiscrimType','pseudoquadratic');
            predictions = predict(mdl,X(test,:));
            estimated_performance=classperf(cp,predictions,test);

            % compute the confusion matrix using stratified 10-fold cross validation
            [C,order] = confusionmat(Y(test),predictions);

            Accuracy(count_time_interval,count_source,count_runs)=sum(diag(C))/sum(sum(C));
            f1_score(count_time_interval,count_source,count_runs)=2*C(1,1)/(2*C(1,1)+C(1,2)+C(2,1));

            % do the same for shuffled labels
            Y_shuffle= label(randperm(length(Y)))';
            X_shuffle= X(randperm(length(X)),:);
            mdl_shuffle = fitcdiscr(X_shuffle(train,:),Y_shuffle(train),'DiscrimType','pseudoquadratic');
            predictions = predict(mdl_shuffle,X_shuffle(test,:));
            %estimated_performance=classperf(cp,predictions,test);

            [C,order] = confusionmat(Y_shuffle(test),predictions);            
            Accuracy_shuffle(count_time_interval,count_source,count_runs)=sum(diag(C))/sum(sum(C));
            f1_score_shuffle(count_time_interval,count_source,count_runs)=2*C(1,1)/(2*C(1,1)+C(1,2)+C(2,1));

        end
       
    end

end

%% plot accuracy

Fs=20; % sampling rate
time= 0:1/Fs:(end_sample-1)/Fs;
source=2;

input_to_plot=[];
input_to_plot.mean=mean(Accuracy(:,source,:),3,'omitnan');
input_to_plot.std=std(Accuracy(:,source,:),[],3,'omitnan')./sqrt(number_of_estimations);

shadedErrorBar(time,input_to_plot.mean,input_to_plot.std,'b',0)
ylim([0.5 0.9])

hold on

input_to_plot=[];
input_to_plot.mean=mean(Accuracy_shuffle(:,source,:),3,'omitnan');
input_to_plot.std=std(Accuracy_shuffle(:,source,:),[],3,'omitnan')./sqrt(number_of_estimations);
% input with shuffled labels
shadedErrorBar(time,input_to_plot.mean,input_to_plot.std,'k',0)

ylim([0.5 0.9])
title([ROIs_names{source,1}])
%end
hold off


%%  classification using the k-nearest neighbor classifier
% cross-validate the model 10 times

for ii=1:length(channel_label)
    input_Hit=[];
    input_CR=[];
   
        for jj=1:mice
            input_Hit=cat(1, input_Hit, outdegree_vector.Hit.multivar{1,ii}{1,jj});
            input_CR=cat(1, input_CR, outdegree_vector.CR.multivar{1,ii}{1,jj});
        end

    input_Hit(isnan(input_Hit))=0;
    input_CR(isnan(input_CR))=0;
   
    sessions=size(cat(1,input_Hit,input_CR),1);
    predictor= cat(1,input_Hit,input_CR);
    % PCA transformed predictor
    [PCALoadings,PCAScores,PCAVar] = pca(predictor);
    label=[ones(1,size(input_Hit,1)) zeros(1,size(input_CR,1))];
   
    Y = label';
    cp = classperf(Y);
    X = PCAScores(:,1:6);
   
    mdl = fitcdiscr(X,Y);
    predictions = predict(mdl,X);
    classperf(cp,predictions);
%     for i = 1:10
%        
%     [train,test] = crossvalind('LeaveMOut',Y,5);
%     mdl = fitcknn(X(train,:),Y(train),'NumNeighbors',3);
%     %mdl = fitglm(X(train,:),Y(train),'Distribution','binomial','Link','logit');
%     predictions = predict(mdl,X(test,:));
%     classperf(cp,predictions,test);
%    
%     end
   
end
cp.ErrorRate

subplot(1,4,1)
imagesc(predictor)
subplot(1,4,2)
imagesc(X)
subplot(1,4,3)
imagesc(Y)
subplot(1,4,4)
imagesc(predictions)
% subplot(1,5,5)
% test_ind=find(test>0);
% Y(test_ind)=predictions;
% imagesc(Y)