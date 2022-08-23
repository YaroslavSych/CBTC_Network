function [ACC,ACC_var,f1]= calculate_ACC_all_channels(test_Hit,test_CR,stimulus,limit_by_trials)

% input dimentions [number_of_trials,trial_time,fiber_channel]
[n1,~,number_of_channels]=size(test_Hit);
[n2,~,~]=size(test_CR);
min_of_trials=min(n1, n2);
    
% mean as a predictor compute the confusion matrix with cross validation
if or(isempty(test_Hit),n1<limit_by_trials)
    ACC = NaN;
    ACC_var = NaN;
    f1= NaN;
elseif or(isempty(test_CR),n2<limit_by_trials)
    ACC = NaN;
    ACC_var = NaN;
    f1= NaN;
else
    predictor=[];
    label=[];


    % if number_of_channels==12
    % 
    %     tmp_Hit= reshape(test_Hit(:,stimulus,1:12),n1,[],1);
    %     tmp_CR= reshape(test_CR(:,stimulus,1:12),n2,[],1);
    % 
    %     label=[ones(n1,1); zeros(n2,1)];
    %     %label_shuffle= label(randperm(length(label)));
    % 
    %     predictor=[tmp_Hit; tmp_CR];
    %     
    % elseif number_of_channels==48
    % end

    number_of_estimations=5;

    for count_runs=1:number_of_estimations
        
        % average df/f across stimulus time interval 
        if n1>n2
            tmp_Hit= datasample(test_Hit(:,stimulus,:),n2,1,'Replace',false);
            tmp_CR= test_CR(:,stimulus,:);
        else
            tmp_Hit= test_Hit(:,stimulus,:);
            tmp_CR= datasample(test_CR(:,stimulus,:),n1,1,'Replace',false);
        end
            
        % concatenate all channels in the predictor
        if number_of_channels==48
            channel_sub_network=[26 31 30 25 34 46 44 47 15 48 33 35];
            tmp_Hit= tmp_Hit(:,:,channel_sub_network);
            tmp_CR= tmp_CR(:,:,channel_sub_network);
        end
        
        % reshape matrix of predictors to use temporal information during
        % stimulus time
        rtmp_Hit=reshape(tmp_Hit,min_of_trials,[],1);
        rtmp_CR=reshape(tmp_CR,min_of_trials,[],1);
        
        label=[ones(min_of_trials,1); zeros(min_of_trials,1)];
        %predictor=[tmp_Hit;tmp_CR]; 
        predictor=[rtmp_Hit; rtmp_CR]; 

        % PCA transformed predictor
        % [PCALoadings,PCAScores,PCAVar] = pca(predictor);
        % X = PCAScores(:,1:3);

        Y = label;
        cp = classperf(Y);
        X = predictor;

        [train,test] = crossvalind('LeaveMOut',Y,limit_by_trials-10);
%         mdl = fitcdiscr(X,Y,'DiscrimType','pseudolinear');
%         predictions = predict(mdl,X);
%         [C,order] = confusionmat(Y,predictions);
%         Accuracy(count_runs)=sum(diag(C))/sum(sum(C));
%         f1_score(count_runs)=2*C(1,1)/(2*C(1,1)+C(1,2)+C(2,1));

        mdl = fitcdiscr(X(train,:),Y(train),'DiscrimType','pseudolinear');
        predictions = predict(mdl,X(test,:));
        estimated_performance=classperf(cp,predictions,test);

        [C,order] = confusionmat(Y(test),predictions);
        Accuracy(count_runs)=sum(diag(C))/sum(sum(C));
        f1_score(count_runs)=2*C(1,1)/(2*C(1,1)+C(1,2)+C(2,1));

        % shuffle labels
%         mdl_shuffle = fitcdiscr(X,Y(randperm(length(Y))),'DiscrimType','pseudolinear');
%         predictions_shuffle = predict(mdl_shuffle,X);
%         [C_shuffle,order] = confusionmat(Y,predictions_shuffle);
%         
%         Accuracy_shuffle(count_runs)=sum(diag(C_shuffle))/sum(sum(C_shuffle));
%         f1_score_shuffle(count_runs)=2*C_shuffle(1,1)/(2*C_shuffle(1,1)+C_shuffle(1,2)+C_shuffle(2,1));

    end

    ACC=mean(Accuracy,2,'omitnan');
    ACC_var=std(Accuracy,[],2,'omitnan')/sqrt(number_of_estimations);
    f1=mean(f1_score,2,'omitnan');
    
%     ACC_shuffle=mean(Accuracy_shuffle,2,'omitnan');
%     ACC_var_shuffle=std(Accuracy_shuffle,[],2,'omitnan')/sqrt(number_of_estimations);
%     f1_shuffle=mean(f1_score_shuffle,2,'omitnan');
%     
%     
%     if ACC<ACC_shuffle
%       ACC=0.5;  
%       f1=f1_shuffle;
%     end
end
end