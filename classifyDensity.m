clear
close all

%% Script uses supervised machine learning algorithms to train classifiers to discriminate between 40 vs. 5 worm videos
% based on extracted Tierpsy features. It also has the option to apply
% sequantial feature selection to identify top features to use for
% classification

%% Specify analysis parameters

% set which feature extraction timestamp to use
featExtractTimestamp = '20200519_153722'; %'20200519_153722' (feat 3006),'20200511_162714' (feat 3006 three windows) or '20191024_122847' (feat 4548) or '20181203_141111' (feat 4548)
if strcmp(featExtractTimestamp,'20200511_162714')
    featExtractWindow = '1'; %'0','1','2'
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
else
    extractStamp = featExtractTimestamp;
end
n_nonFeatVar = 17;

% select which features to drop or keep
feats2drop = {'path','blob'};
applyKeyFeaturesForClassifierTraining = false; % apply pre-determined key features

% select which tasks to perform
performSequentialFeatureSelection = true;
trainClassifier = false;

if performSequentialFeatureSelection
    % Note: Currently using quadratic discriminant analysis for sfs.
    % Otherwise redefine classf function.
    crossVal_k = 5;
    n_sortedFeatsInput = 500; % 'NaN' or use top '500' features with lowest t-test p-values for sfs to make it run faster.
    n_topFeatsOutput = 10;
    plotSFS = false; % generte additional diagnostic plots for SFS
    if applyKeyFeaturesForClassifierTraining
        error('applyKeyFeaturesForClassifierTraining must be set to false for sequentialFeatureSelection results to be used for classification.')
    end
end

%% Load and process featureTable
% load
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
% trim table down to retain necessary info
wormNums = featureTable.wormNum; % get classification labels
strainNames = featureTable.strain_name; % get classification labels
featureTable = featureTable(:,n_nonFeatVar+1:end); % get full features matrix
% shortens the names of features
featureTable = shortenFeatNamesInFeatTable(featureTable);
% drop features as specified
if ~isempty(feats2drop)
    [featureTable, ~] = dropFeats(featureTable,feats2drop);
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

%% Holdout partition to separate training (75%) and test (25%) set
holdoutCVP = cvpartition(wormNums,'holdout',0.25);
dataTrain = featureMat(holdoutCVP.training,:);
grpTrain = wormNums(holdoutCVP.training,:);
dataTest = featureMat(holdoutCVP.test,:);
grpTest = wormNums(holdoutCVP.test,:);

%% Sequential feature selection
% following instructions from https://uk.mathworks.com/help/stats/examples/selecting-features-for-classifying-high-dimensional-data.html

if performSequentialFeatureSelection
    
    % ecdf plot of t-test values
    dataTrainG1 = dataTrain(grp2idx(grpTrain)==1,:);
    dataTrainG2 = dataTrain(grp2idx(grpTrain)==2,:);
    [h,p,ci,stat] = ttest2(dataTrainG1,dataTrainG2,'Vartype','unequal');
    if plotSFS
        ecdf(p);
        xlabel('P value');
        ylabel('CDF value')
    end
    
    % sort features
    [~,featureIdxSortbyP] = sort(p,2); 
    
    % define function for sfs. Quadratic discriminant analysis currently but should probably change for favoured model. 
    classf = @(xtrain,ytrain,xtest,ytest) ...
        sum(~strcmp(ytest,classify(xtest,xtrain,ytrain,'quadratic'))); 
    
    % sequential feature selection with k-fold cross-validation
    fivefoldCVP = cvpartition(grpTrain,'kfold',crossVal_k);
    if ~isnan(n_sortedFeatsInput)
        featureIdxSortbyP = featureIdxSortbyP(1:n_sortedFeatsInput);
    end
    fsLocal = sequentialfs(classf,dataTrain(:,featureIdxSortbyP),grpTrain,'cv',fivefoldCVP,'nfeatures',n_topFeatsOutput);
    featCols = featureIdxSortbyP(fsLocal);
    keyFeats = featureTable.Properties.VariableNames(featCols)'
    
    if plotSFS
        % evaluate the performance of the selected model with the features
        testMCELocal = crossval(classf,featureMat(:,featureIdxSortbyP(fsLocal)),wormNums,'partition',...
            holdoutCVP)/holdoutCVP.TestSize;
        
        % plot of the cross-validation MCE as a function of the number of features for up to 10 features
        [fsCVfor10,historyCV] = sequentialfs(classf,dataTrain(:,featureIdxSortbyP),grpTrain,...
            'cv',fivefoldCVP,'nfeatures',n_topFeatsOutput);
        plot(historyCV.Crit,'o');
        xlabel('Number of Features');
        ylabel('CV MCE');
        title('Forward Sequential Feature Selection with cross-validation');
        
        % plot of resubstitution MCE values on the training set
        % (i.e., without performing cross-validation during the feature selection procedure)
        % as a function of the number of features:
        [fsResubfor10,historyResub] = sequentialfs(classf,dataTrain(:,featureIdxSortbyP),...
            grpTrain,'cv','resubstitution','nfeatures',n_topFeatsOutput);
        plot(1:n_topFeatsOutput, historyCV.Crit,'bo',1:n_topFeatsOutput, historyResub.Crit,'r^');
        xlabel('Number of Features');
        ylabel('MCE');
        legend({'10-fold CV MCE' 'Resubstitution MCE'},'location','NE'); %%%% 10 fold? 5 fold? ?? don't understand this plot
    end
end

%% Select pre-determined key features (this overrides keyFeats that may have been generated during sequential feature selection).

if applyKeyFeaturesForClassifierTraining
    
% use key features from sequential feature selection (without path and blob features)
keyFeats = {'turn_intra_frequency','length_norm_IQR','length_norm_10th',...
    'rel_to_head_base_radial_vel_head_tip_w_backward_50th','orientation_food_edge_w_forward_IQR',...
    'turn_intra_duration_50th','motion_mode_paused_fraction','orientation_food_edge_in_edge_IQR',...
    'orientation_food_edge_in_edge_10th','orientation_food_edge_w_backward_IQR'};
     
% % use key features from compareSwFeatDensityEffect.m
% keyFeats = {'motion_mode_forward_duration_50th','motion_mode_forward_fraction','motion_mode_paused_fraction',...
%         'food_region_inside_duration_50th','food_region_edge_duration_50th','turn_intra_duration_50th',...
%         'turn_intra_frequency','rel_to_head_base_radial_vel_head_tip_w_backward_50th',...
%         'd_speed_w_backward_50th','orientation_food_edge_w_forward_IQR'};

end

%% Slice out key features from the full featureTable
if performSequentialFeatureSelection | applyKeyFeaturesForClassifierTraining
    
    featCols = [];
    for keyFeatCtr = 1:numel(keyFeats)
        keyFeat = keyFeats{keyFeatCtr};
        keyFeat = string(keyFeat);
        featCols = [featCols, find(strcmp(featureTable.Properties.VariableNames,keyFeat))];
    end
    featureTable = featureTable(:,featCols);
    disp([num2str(numel(featCols)) ' key features retained.'])
end

%% Train classifier

if trainClassifier
    
    % Add wormNum to featureTable for classification
    featureTable.wormNum = wormNums;
    
    % Extract features matrix for training and test sets
    trainFeatureTable = featureTable(holdoutCVP.training,:);
    testFeatureTable = featureTable(holdoutCVP.test,:);
    
    % Model selection using classificationLearner GUI
    disp('Please use GUI to train model with trainFeatureTable. Use five-fold cross validation, select the best model, and export to workspace.')
    classificationLearner
    
    % Wait until trained model has been selected and saved
    while exist('trainedModel')~=1
        pause
    end
    
    % Classification accuracy on unseen hold-out test set
    yfit = trainedModel.predictFcn(testFeatureTable);
    accuracy = nnz(yfit == testFeatureTable.wormNum)/holdoutCVP.TestSize;
    disp(['Accuracy from trained model is ' num2str(accuracy) ' on unseen test data.'])
end