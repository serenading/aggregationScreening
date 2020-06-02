clear
close all

%% Script uses supervised machine learning algorithms to train classifiers for a specified variable
% based on extracted Tierpsy features. It also has the option to apply
% sequantial feature selection to identify top features to use for
% classification

% author: serenading. June 2020.

%% Specify analysis parameters

% set which variable to classify for
classVar = 'strain_name'; % 'wormNum','strain_name'. Must be a variable field of the featureTable. 
n_nonFeatVar = 17;

% set which feature extraction timestamp to use
featExtractTimestamp = '20200519_153722'; %'20200519_153722' (feat 3016),'20200511_162714' (feat 3016 three windows) or '20191024_122847' (feat 4548) or '20181203_141111' (feat 4548)
if strcmp(featExtractTimestamp,'20200511_162714')
    featExtractWindow = '1'; %'0','1','2'
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
else
    extractStamp = featExtractTimestamp;
end

% select which features and features to drop or keep
strains2keep = {'all'}; % Use all strains if cell left empty. {'all'} or {'divergent'} or {'controls'} or {'strain1', 'strain2'}. Cell array containing strains to keep for analysis. 
strains2drop = {}; % {'N2','CB4856','DA609','ECA252','LSJ1'}; Cell array containing strains to drop from analysis.
feats2keep = {'Tierpsy_256'}; % {'Tierpsy_256'} or {'feat1','feat2'}. Cell array containing features to use for analysis. Use all features if left empty.
feats2drop = {}; % {'path'}; % Cell array containing features to drop from analysis. Partial name of feature allowed. 
applyKeyFeaturesForClassifierTraining = false; % apply pre-determined key features

% select which tasks to perform
performSequentialFeatureSelection = false;
trainClassifier = true;

if performSequentialFeatureSelection
    % Note: Currently using linear discriminant analysis for sfs.
    % Otherwise redefine classf function.
    crossVal_k = 5;
    if isempty(feats2keep)
        n_sortedFeatsInput = 500; % 'NaN' or use top '500' features with lowest t-test p-values for sfs to make it run faster.
    end
    n_topFeatsOutput = 10;
    plotSFS = false; % generte additional diagnostic plots for SFS
    if applyKeyFeaturesForClassifierTraining
        error('applyKeyFeaturesForClassifierTraining must be set to false for sequentialFeatureSelection results to be used for classification.')
    end
end

%% Load and process featureTable
% load
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);

% check that the specified variable name exists
if nnz(contains(featureTable.Properties.VariableNames,classVar))~=1
    error('Invalid classVar name. Variable name must exist inside featureTable.')
end

% use only 5 worm data to classify strains
if strcmp(classVar,'strain_name')
    fivewormLogInd = featureTable.wormNum == 5;
    featureTable = featureTable(fivewormLogInd,:);
end

% filter featureTable based on specified strain and features
[featureTable, classLabels] = filterFeatureTable(featureTable,classVar,n_nonFeatVar,strains2keep,strains2drop,feats2keep,feats2drop);

%% Pre-process features matrix and turn back into table 
% split table into matrix and featNames
featureMat = table2array(featureTable);
featNames = featureTable.Properties.VariableNames;
% preprocess feature matrix: drop zero standard deviation, NaN's, z-normalise, etc. 
[featureMat,dropLogInd] = preprocessFeatMat(featureMat);
featNames = featNames(~dropLogInd);
% put the table back together
featureTable = array2table(featureMat,'VariableNames',featNames);

%% Holdout partition to separate training (75%) and test (25%) set
holdoutCVP = cvpartition(classLabels,'holdout',0.25);
dataTrain = featureMat(holdoutCVP.training,:);
grpTrain = classLabels(holdoutCVP.training,:);
dataTest = featureMat(holdoutCVP.test,:);
grpTest = classLabels(holdoutCVP.test,:);

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
    
    % define function for sfs. Linear discriminant analysis currently but should probably change for favoured model. 
    classf = @(xtrain,ytrain,xtest,ytest) ...
        sum(~strcmp(ytest,classify(xtest,xtrain,ytrain,'linear'))); 
    
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
        testMCELocal = crossval(classf,featureMat(:,featureIdxSortbyP(fsLocal)),classLabels,'partition',...
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
    
% % use key features from sequential feature selection (without path and blob features)
% keyFeats = {'turn_intra_frequency','length_norm_IQR','length_norm_10th',...
%     'rel_to_head_base_radial_vel_head_tip_w_backward_50th','orientation_food_edge_w_forward_IQR',...
%     'turn_intra_duration_50th','motion_mode_paused_fraction','orientation_food_edge_in_edge_IQR',...
%     'orientation_food_edge_in_edge_10th','orientation_food_edge_w_backward_IQR'};

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
    
    % Add classLabels to featureTable for classification
    featureTable.(classVar) = classLabels;
    
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
    if strcmp(classVar,'strain_name')
        accuracy = nnz(strcmp(yfit,testFeatureTable.(classVar)))/holdoutCVP.TestSize;
    elseif strcmp(classVar,'wormNum')
        accuracy = nnz(yfit == testFeatureTable.(classVar))/holdoutCVP.TestSize;
    else
        error('Please specify how accuracy should be assessed for this classVar.')
    end
    disp(['Test accuracy from trained model is ' num2str(accuracy) ' on unseen data.'])
end