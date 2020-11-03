clear
close all

%% Script uses supervised machine learning algorithms to train classifiers for a specified variable
% based on extracted Tierpsy features. It also has the option to apply
% sequantial feature selection to identify top features to use for
% classification

% author: @serenading. June 2020.

%% Specify analysis parameters

% Set which variable to classify for
classVar = 'strain_name'; % 'wormNum' or 'strain_name'. Must be a variable field of the featureTable.
wormNum = 5; % only taken into account if classVar is 'strain_name'.

% Set which feature extraction timestamp to use
featExtractTimestamp = '20200630_171156'; %'20201021_114135' (augmented feat 3016),'20200630_171156' (augmented feat 3016 windows),'20200519_153722' (feat 3016),'20200511_162714' (feat 3016 three windows) or '20191024_122847' (feat 4548)
if strcmp(featExtractTimestamp,'20200511_162714')
    featExtractWindow = '1'; %'0','1','2'. Windows: window0: 0-15 min, window1: 15-30 min, window2: 30-45 min.
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
    n_nonFeatVar = 17;
elseif strcmp(featExtractTimestamp,'20201021_114135') % Augmented with 5 fold at 0.8 trajectory ratio, up to 10 min.
    extractStamp = ['augmented_' featExtractTimestamp];
    n_nonFeatVar = 19;
elseif strcmp(featExtractTimestamp,'20200630_171156') % Augmented with 5 fold at 0.8 trajectory ratio, up to 30 min.
    featExtractWindow = '4'; %'0','1','2','3','4'. Windows: window0: 0-15 min window1: 15-30 min, window2: 30-45 min, window3: 15-45 min, window 4: 0-45 min.
    extractStamp = ['augmented_' featExtractTimestamp '_window_' featExtractWindow];
    n_nonFeatVar = 19;
else
    extractStamp = featExtractTimestamp;
    n_nonFeatVar = 17;
end

% Select which features and features to drop or keep
strains2keep = {}; % Use all strains if cell left empty. {'all'} or {'divergent'} or {'controfeatureTablels'} or {'strain1', 'strain2'} or {'swept_liberal'} or {'nonSwept_liberal'} or {'swept_conservative'} or {'nonSwept_conservative'}. Cell array containing strains to keep for analysis.
strains2drop = {'DA609'}; % {'N2','CB4856','DA609','ECA252','LSJ1'}; Cell array containing strains to drop from analysis.
feats2keep = {}; % Use all features if left empty. {'Tierpsy_256'} or {'feat1','feat2'}. Cell array containing features to use for analysis.
feats2drop = {}; % {'path'}; % Cell array containing features to drop from analysis. Partial name of feature allowed.

% Select which tasks to perform
applyKeyFeaturesForClassifierTraining = true; % apply pre-determined key features. Note for classification on truly unseen data you possibly want to choose features using SFS each time to train a new model as opposed to using pre-determined feats?
if applyKeyFeaturesForClassifierTraining
    keyFeats2load = 'feats_all'; % Options:'feats_all','feats_swept_liberal','feats_nonSwept_liberal','feats_swept_conservative','feats_nonSwept_conservative','feats_top_H2','feats_top_h2'
end
performSequentialFeatureSelection = false;
trainClassifier = true;

%% Optional parameters

% Randomly sub-select a number of strains from the full panel
if isempty(strains2keep) || strcmp(strains2keep,'all')
    useRandomStrains = false;
    n_RandomStrains = 8;
end

% Randomly subsample replicates so that all strains have the same number of replicates
subsampleReps = true;

% SFS parameters. Note: Currently using linear discriminant analysis for sfs. Otherwise redefine classf function.
if performSequentialFeatureSelection
    % Generate additional diagnostic plots for SFS
    plotSFS = true;
    if isempty(feats2keep)
        % Use top n sorted features with lowest p-values to make sfs run faster.
        % Also if n_feat too high for n_obs then some classifiers do not work.
        % Likely needs tweaking based on input observations.
        n_sortedFeatsInput = 500;
    else
        n_sortedFeatsInput = NaN;
    end
    if strcmp(classVar,'wormNum')
        n_topFeatsOutput = 10;
    elseif strcmp(classVar,'strain_name')
        n_topFeatsOutput = 30;
    else
        n_topFeatsOutput = 200;
    end
end

% Classification parameters. How many replicates to use for hold out and cross validation
if strcmp(classVar,'strain_name') && wormNum ==40
    holdoutRatio = 1/5; % 1 out of (typically) 5 replicates held out
    crossVal_k = 4;
else
    holdoutRatio = 0.25;
    crossVal_k = 5;
    if strcmp(classVar,'strain_name') && wormNum ==5
        holdoutRatio = 1/3; % 1 out of (typically) 3 replicates held out
        crossVal_k = 2;
        if strcmp(featExtractTimestamp,'20201021_114135') & strcmp(featExtractTimestamp,'20200630_171156') % if using augmented data
            crossVal_k = 5;
        end
    end
end

%% Load and process featureTable

% Load featuresTable and set classification variables according to the task
% 40 worm manual features for strain classification
if strcmp(classVar,'strain_name') && wormNum ==40
    % load 40 worm manually extracted features
    featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_20200519_153722_new_20200620.csv'],'Delimiter',',','preserveVariableNames',true);
else
    %  5 and 40 worm Tierpsy features for density classification
    featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
    if strcmp(classVar,'strain_name') && wormNum ==5
        % 5 worm Tierpsy features for strain classification
        fivewormLogInd = featureTable.wormNum == 5;
        featureTable = featureTable(fivewormLogInd,:);
    end
end

% Load the full list of strains
load('strainsList/all.mat');

% Check that the specified variable name exists
if nnz(contains(featureTable.Properties.VariableNames,classVar))~=1
    error('Invalid classVar name. Variable name must exist inside featureTable.')
end

% If specified, subsample replicates to ensure all strains have the same number of reps
if strcmp(classVar,'strain_name') && subsampleReps
    featureTable = balanceStrainReps(featureTable,featExtractTimestamp,crossVal_k,strains);
end

% If specified, use a randomly sub-selected panel of strains from the full list
if (isempty(strains2keep) || strcmp(strains2keep,'all')) && useRandomStrains
    randInd = randperm(numel(strains),n_RandomStrains);
    strains2keep = strains(randInd);
end

% Filter featureTable based on specified strain and features
[featureTable, classLabels] = filterFeatureTable(featureTable,classVar,n_nonFeatVar,strains2keep,strains2drop,feats2keep,feats2drop);

% Apply pre-determined keyFeatures if selected
if applyKeyFeaturesForClassifierTraining
    load('keyFeats.mat',keyFeats2load)
    keyFeats = eval(genvarname(keyFeats2load)); % renames the loaded features as "keyFeats"
    % trim down to keyFeats
    featureTable = featureTable(:,keyFeats);
end

%% Pre-process features matrix and turn back into table

% Split table into matrix and featNames
featureMat = table2array(featureTable);
featNames = featureTable.Properties.VariableNames;
% Preprocess feature matrix: drop zero standard deviation, NaN's, z-normalise, etc.
[featureMat,dropLogInd] = preprocessFeatMat(featureMat);
featNames = featNames(~dropLogInd);
% Put the table back together
featureTable = array2table(featureMat,'VariableNames',featNames);

%% Holdout partition to separate training and test set

% Generate test/training partition
if strcmp(featExtractTimestamp,'20201021_114135') | strcmp(featExtractTimestamp,'20200630_171156')
    % If using augmented data, make sure one entire replicate is held out so it is not at all seen in training. 
    holdoutCVP = cvpartition_augmented(classLabels);
else 
    % If not using augmented data, use cvpartition to divide up data into training and test sets
    holdoutCVP = cvpartition(classLabels,'holdout',holdoutRatio);
end

% Extract training and test datasets
dataTrain = featureMat(holdoutCVP.training,:);
grpTrain = classLabels(holdoutCVP.training,:);
dataTest = featureMat(holdoutCVP.test,:);
grpTest = classLabels(holdoutCVP.test,:);

%% Sequential feature selection
% Following instructions from https://uk.mathworks.com/help/stats/examples/selecting-features-for-classifying-high-dimensional-data.html

if performSequentialFeatureSelection
    
    % ecdf plot of t-test or anova test values
    grpTrainInd = grp2idx(grpTrain);
    if strcmp(classVar,'wormNum')
        % t-test for two density class labels: 5 vs. 40
        dataTrainG1 = dataTrain(grpTrainInd==1,:);
        dataTrainG2 = dataTrain(grpTrainInd==2,:);
        [~,p] = ttest2(dataTrainG1,dataTrainG2,'Vartype','unequal');
    elseif strcmp(classVar,'strain_name')
        % One-way anova for multiple class labels: strain names
        p = NaN(1,size(dataTrain,2));
        for featCtr = 1:size(dataTrain,2)
            p(featCtr) = anova1(dataTrain(:,featCtr),grpTrainInd,'off');
        end
    else
        error('Appropriate statistical test not specified to rank features for this classVar')
    end
    if plotSFS
        figure; ecdf(p);
        xlabel('P value');
        ylabel('CDF value')
    end
    
    % Sort features
    [~,featureIdxSortbyP] = sort(p,2);
    
    % Define function for sfs
    
    % subspace discriminant model - TODO: needs to write this as a proper function to enable multiple lines
    %     Mdl = fitcensemble(xtrain,ytrain,'Method','subspace'); % trained object using subspace discriminant model.
    %     classf = @(xtrain,ytrain,xtest,ytest) ...
    %         sum(~strcmp(ytest,fitcensemble(xtrain,ytrain,'Method','subspace').predictFcn(xtest)));
    
    % Linear discriminant model
    classf = @(xtrain,ytrain,xtest,ytest) ...
        sum(~strcmp(ytest,classify(xtest,xtrain,ytrain,'linear')));
    
    % Sequential feature selection with k-fold cross-validation
    kfoldCVP = cvpartition(grpTrain,'kfold',crossVal_k);
    if ~isnan(n_sortedFeatsInput)
        featureIdxSortbyP = featureIdxSortbyP(1:n_sortedFeatsInput);
    end
    [fsLocal,historyCV] = sequentialfs(classf,dataTrain(:,featureIdxSortbyP),grpTrain,'cv',kfoldCVP,'nfeatures',n_topFeatsOutput);
    featCols = featureIdxSortbyP(fsLocal);
    keyFeats = featureTable.Properties.VariableNames(featCols)'
    
    if plotSFS
        % Evaluate the performance of the selected model with the features
%         testMCELocal = crossval(classf,featureMat(:,featureIdxSortbyP(fsLocal)),classLabels,'partition',...
%             holdoutCVP)/holdoutCVP.TestSize;
        figure;
        
        % Plot of the cross-validation MCE as a function of the number of top features
        subplot(1,2,1)
        plot(historyCV.Crit,'o');
        xlabel('Number of Features');
        ylabel('Cross-validation Misclassification Error');
        ylim([0 1])
        
        % Plot the derivative of the cross-validation MCE as a function of the number of top features
        dCrit = [historyCV.Crit 0]-[1 historyCV.Crit];
        subplot(1,2,2)
        plot(dCrit,'o');
        xlim([0 numel(historyCV.Crit)])
        xlabel('Number of Features');
        ylabel('Change in CV MSE');
        
        %         % plot of resubstitution MCE values on the training set
        %         % (i.e., without performing cross-validation during the feature selection procedure)
        %         % as a function of the number of top features:
        %         [fsResub,historyResub] = sequentialfs(classf,dataTrain(:,featureIdxSortbyP),grpTrain,...
        %             'cv','resubstitution','nfeatures',n_topFeatsOutput);
        %         figure; plot(1:n_topFeatsOutput, historyCV.Crit,'bo',1:n_topFeatsOutput, historyResub.Crit,'r^');
        %         xlabel('Number of Features');
        %         ylabel('MCE');
        %         legend({[num2str(crossVal_k) '-fold CV MCE'], 'Resubstitution MCE'},'location','NE');
    end
    
    % Slice out key features from the full featureTable
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
    
    % Divide featureTable into training and test sets
    trainFeatureTable = featureTable(holdoutCVP.training,:);
    testFeatureTable = featureTable(holdoutCVP.test,:);
    
    % Perform model selection using classificationLearner GUI
    disp(['Please use GUI to train model with trainFeatureTable. Use ' num2str(crossVal_k) '-fold cross validation, select the best model, and export to workspace.'])
    classificationLearner
    
    % Wait until trained model has been selected and saved
    while exist('trainedModel')~=1
        pause
    end
    
    % For classifying strains,  Ensemble subspace discrimnant and linear discriminant models tend to work the best.
    % For classifying density, coarse tree model works the best but many models perform very well.
    
    %% Classification accuracy on unseen hold-out test set
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