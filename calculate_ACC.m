function [ACC,ACC_var,f1]= calculate_ACC(test_Hit,test_CR,stimulus,count_channel,limit_by_trials)

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

    number_of_estimations=30;

    for count_runs=1:number_of_estimations
        
            if n1>n2
                tmp_Hit= datasample(test_Hit(:,stimulus,count_channel),n2,1,'Replace',false);
                tmp_CR= test_CR(:,stimulus,count_channel);
            else
                tmp_Hit= test_Hit(:,stimulus,count_channel);
                tmp_CR= datasample(test_CR(:,stimulus,count_channel),n1,1,'Replace',false);
            end

        label=[ones(min_of_trials,1); zeros(min_of_trials,1)];
        predictor=[tmp_Hit;tmp_CR]; 
        %predictor=[mean(tmp_Hit,2); mean(tmp_CR,2)]; 

        % PCA transformed predictor
        % [PCALoadings,PCAScores,PCAVar] = pca(predictor);
        % X = PCAScores(:,1:3);

        Y = label;
        cp = classperf(Y);
        X = predictor;

        % --- use Linear Discriminant Analysis for Classification ---
%         mdl = fitcdiscr(X,Y,'DiscrimType','pseudolinear');
%         predictions = predict(mdl,X);
%         [C,order] = confusionmat(Y,predictions);
%         Accuracy(count_runs)=sum(diag(C))/sum(sum(C));
%         f1_score(count_runs)=2*C(1,1)/(2*C(1,1)+C(1,2)+C(2,1));
        
        % --- use Random Forest for Classification instead of LDA ---
%         [train,test] = crossvalind('LeaveMOut',Y,limit_by_trials-10);
%         rng(1); % For reproducibility
%         mdl = TreeBagger(100,X(train,:),Y(train),'OOBPrediction','On',...
%             'Method','classification');
%         [C,order] = confusionmat(Y(test),str2double(predictions));
%         Accuracy(count_runs)=sum(diag(C))/sum(sum(C));
%         f1_score(count_runs)=2*C(1,1)/(2*C(1,1)+C(1,2)+C(2,1));

        % --- use Linear Discriminant Analysis with cross-validation ---
        [train,test] = crossvalind('LeaveMOut',Y,limit_by_trials-10);
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