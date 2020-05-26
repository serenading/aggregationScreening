clear
close all

%% Script uses supervised machine learning algorithms to train classifiers to discriminate between 40 vs. 5 worm videos
% based on extracted Tierpsy features. It also has the option to apply
% sequantial feature selection to identify top features to use for
% classification

% Todo: use the same holdout for sequential feature selection as for hold
% out validation

%% Specify analysis parameters

% set which feature extraction timestamp to use
featExtractTimestamp = '20200511_162714'; %'20200511_162714' or '20191024_122847' or '20181203_141111'
featExtractWindow = '1'; % 'none','0','1','2'
n_nonFeatVar = 17;

% select which features to drop or keep
dropPathBlobFeatures = true;
applyKeyFeaturesForClassifierTraining = false;

% select which tasks to perform
trainClassifier = false;
performSequentialFeatureSelection = true;

% get date time window stamp
if strcmp(featExtractWindow,'none')
    extractStamp = featExtractTimestamp;
else
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
end

if performSequentialFeatureSelection
    if applyKeyFeaturesForClassifierTraining
        error('applyKeyFeaturesForClassifierTraining must be set to false for sequentialFeatureSelection.')
    end
end

%% Load and process featureTable
% load
featureTable = readtable(['/Users/sding/Dropbox/aggScreening/results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
% trim table down to retain necessary info
wormNums = featureTable.wormNum; % get classification labels
featureTable = featureTable(:,n_nonFeatVar+1:end); % get full features matrix
% shortens the names of features
featureTable = shortenFeatNamesInFeatTable(featureTable);
% Drop path-related and blob features as these are already expected to be good predictors for this task due to tracking artefacts
if dropPathBlobFeatures
    featuresTable = dropPathBlobFeats(featuresTable);
end

%% Pre-process features matrix and turn back into table
% drop NaN's, z-normalise, etc. from feature table
featureMat = table2array(featureTable);
featNames = featureTable.Properties.VariableNames;
[featureMat,droppedCols] = preprocessFeatMat(featureMat);
% remove corresponding feature names for dropped features
featNamesLogInd = true(size(featNames));
featNamesLogInd(droppedCols) = false;
featNames = featNames(featNamesLogInd);
% put the table back together
featureTable = array2table(featureMat,'VariableNames',featNames);

%% Select pre-determined key features

if applyKeyFeaturesForClassifierTraining
    
    % use key features from sequential feature selection (without path and blob features)
    keyFeats = {'turn_intra_frequency','length_norm_IQR','length_norm_10th','orientation_food_edge_w_forward_IQR',...
        'motion_mode_paused_fraction','orientation_food_edge_in_edge_IQR','orientation_food_edge_w_backward_IQR',...
        'd_rel_to_head_base_radial_vel_head_tip_w_backward_IQR','motion_mode_forward_fraction','dist_from_food_edge_w_backward_90th'};
    
    % % use key features from sequential feature selection (from all features)
    % keyFeats = {'path_coverage_head','turn_intra_frequency','path_coverage_body','path_coverage_midbody'...
    %     'path_coverage_tail','path_coverage_head_norm','path_coverage_body_norm'...
    %     'path_coverage_midbody_norm','path_coverage_tail_norm','blob_box_width_w_paused_90th'};
    
    % % use key features from compareSwFeatDensityEffect.m
    % keyFeats = {'motion_mode_forward_duration_50th','motion_mode_forward_fraction','motion_mode_paused_fraction',...
    %     'food_region_inside_duration_50th','food_region_edge_duration_50th','turn_intra_duration_50th',...
    %     'turn_intra_frequency','rel_to_head_base_radial_vel_head_tip_w_backward_50th',...
    %     'd_speed_w_backward_50th','orientation_food_edge_w_forward_IQR'};
    
    % slice out key features from the full featureTable
    featCols = [];
    for keyFeatCtr = 1:numel(keyFeats)
        keyFeat = keyFeats{keyFeatCtr};
        featCols = [featCols, find(strcmp(featureTable.Properties.VariableNames,keyFeat))];
    end
    featureTable = featureTable(:,featCols);
    disp([num2str(numel(featCols)) ' key features retained for classification training'])
end

%% Append wormNum classification labels back into the table for cvpartition
featureTable.wormNum = wormNums;

if trainClassifier
    %% Partition dataset into training set (75% of files) and holdout test set (25% of files)
    c = cvpartition(featureTable.wormNum,'HoldOut',0.25);
    trainIdx = gather(c.training);
    testIdx = gather(c.test);
    
    %% Extract features matrix for training and test sets
    trainFeatureTable = featureTable(trainIdx,:);
    testFeatureTable = featureTable(testIdx,:);
    
    %% Model selection using classificationLearner GUI
    classificationLearner
    % use wormnum as response labels and all features as predictors.
    % use five-fold cross validation with training dataset.
    % select model and export
    
    %% Test hold-out validation score
    %     yfit = trainedModel.predictFcn(testFeatureTable);
    %     accuracy = nnz(yfit == testFeatureTable.wormNum)/c.TestSize;
end

%% Sequential feature selection
% following instructions from https://uk.mathworks.com/help/stats/examples/selecting-features-for-classifying-high-dimensional-data.html

if performSequentialFeatureSelection & ~applyKeyFeaturesForClassifierTraining
    % name variables to match that of matlab page
    grp = featureTable.wormNum;
    obs = table2array(featureTable(:,1:end-1));
    
    % hold out partition
    holdoutCVP = cvpartition(grp,'holdout',501);
    dataTrain = obs(holdoutCVP.training,:);
    grpTrain = grp(holdoutCVP.training,:);
    
    % ecdf plot of t test values
    dataTrainG1 = dataTrain(grp2idx(grpTrain)==1,:);
    dataTrainG2 = dataTrain(grp2idx(grpTrain)==2,:);
    [h,p,ci,stat] = ttest2(dataTrainG1,dataTrainG2,'Vartype','unequal');
    ecdf(p);
    xlabel('P value');
    ylabel('CDF value')
    
    % sort features
    [~,featureIdxSortbyP] = sort(p,2); % sort the features
    classf = @(xtrain,ytrain,xtest,ytest) ...
        sum(~strcmp(ytest,classify(xtest,xtrain,ytrain,'quadratic')));
    
    % sequential feature selection with k-fold cross-validation
    fivefoldCVP = cvpartition(grpTrain,'kfold',5);
    fs1 = featureIdxSortbyP(1:500);
    fsLocal = sequentialfs(classf,dataTrain(:,fs1),grpTrain,'cv',fivefoldCVP,'Nf',10);
    featCols = fs1(fsLocal);
    keyFeats = featureTable.Properties.VariableNames(featCols)'
    
    % evaluate the performance of the selected model with the features
    testMCELocal = crossval(classf,obs(:,fs1(fsLocal)),grp,'partition',...
        holdoutCVP)/holdoutCVP.TestSize;
    
    % plot of the cross-validation MCE as a function of the number of features for up to 50 features
    [fsCVfor10,historyCV] = sequentialfs(classf,dataTrain(:,fs1),grpTrain,...
        'cv',fivefoldCVP,'Nf',10);
    plot(historyCV.Crit,'o');
    xlabel('Number of Features');
    ylabel('CV MCE');
    title('Forward Sequential Feature Selection with cross-validation');
    
    % plot of resubstitution MCE values on the training set
    % (i.e., without performing cross-validation during the feature selection procedure)
    % as a function of the number of features:
    [fsResubfor50,historyResub] = sequentialfs(classf,dataTrain(:,fs1),...
        grpTrain,'cv','resubstitution','Nf',10);
    plot(1:10, historyCV.Crit,'bo',1:10, historyResub.Crit,'r^');
    xlabel('Number of Features');
    ylabel('MCE');
    legend({'10-fold CV MCE' 'Resubstitution MCE'},'location','NE');
end