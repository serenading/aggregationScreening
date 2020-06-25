clear
close all

%% Script uses supervised machine learning algorithms to train classifiers for a specified variable
% based on extracted Tierpsy features. It also has the option to apply
% sequantial feature selection to identify top features to use for
% classification

% author: @serenading. June 2020.

%% Specify analysis parameters

% set which variable to classify for
classVar = 'strain_name'; % 'wormNum' or 'strain_name'. Must be a variable field of the featureTable.
wormNum = 5; % only taken into account if classVar is 'strain_name'.
n_nonFeatVar = 17;

% set which feature extraction timestamp to use
featExtractTimestamp = '20200519_153722'; %'20200519_153722' (feat 3016),'20200511_162714' (feat 3016 three windows) or '20191024_122847' (feat 4548)
if strcmp(featExtractTimestamp,'20200511_162714')
    featExtractWindow = '1'; %'0','1','2'
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
else
    extractStamp = featExtractTimestamp;
end

% select which features and features to drop or keep
strains2keep = {}; % Use all strains if cell left empty. {'all'} or {'divergent'} or {'controls'} or {'strain1', 'strain2'}. Cell array containing strains to keep for analysis.
strains2drop = {}; % {'N2','CB4856','DA609','ECA252','LSJ1'}; Cell array containing strains to drop from analysis.
feats2keep = {}; % Use all features if left empty. {'Tierpsy_256'} or {'feat1','feat2'}. Cell array containing features to use for analysis.
feats2drop = {}; % {'path'}; % Cell array containing features to drop from analysis. Partial name of feature allowed.

% select which tasks to perform
applyKeyFeaturesForClassifierTraining = false; % apply pre-determined key features
performSequentialFeatureSelection = true;
trainClassifier = true;

%% Optional parameters
% randomly sub-select a number of strains from the full panel
if isempty(strains2keep) || strcmp(strains2keep,'all')
    useRandomStrains = false;
    n_RandomStrains = 13;
end

% SFS parameters. Note: Currently using linear discriminant analysis for sfs. Otherwise redefine classf function.
if performSequentialFeatureSelection
    % generate additional diagnostic plots for SFS
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

%% Load and process featureTable
% load featuresTable and set classification variables according to the task
% 40 worm manual features for strain classification
if strcmp(classVar,'strain_name') && wormNum ==40
    % load 40 worm manually extracted features
    featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_20200519_153722_new_20200620.csv'],'Delimiter',',','preserveVariableNames',true);
    holdoutRatio = 1/5; % 1 out of (typically) 5 replicates held out
    crossVal_k = 4;
else
    %  5 and 40 worm Tierpsy features for density classification
    featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
    holdoutRatio = 0.25;
    crossVal_k = 5;
    if strcmp(classVar,'strain_name') && wormNum ==5
        % 5 worm Tierpsy features for strain classification
        fivewormLogInd = featureTable.wormNum == 5;
        featureTable = featureTable(fivewormLogInd,:);
        holdoutRatio = 1/3; % 1 out of (typically) 3 replicates held out
        crossVal_k = 2;
    end
end

% check that the specified variable name exists
if nnz(contains(featureTable.Properties.VariableNames,classVar))~=1
    error('Invalid classVar name. Variable name must exist inside featureTable.')
end

% if specified, use a randomly sub-selected panel of strains from the full list
if (isempty(strains2keep) || strcmp(strains2keep,'all')) && useRandomStrains
    load('strainsList/all.mat');
    randInd = randperm(numel(strains),n_RandomStrains);
    strains2keep = strains(randInd);
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

%% Holdout partition to separate training and test set
holdoutCVP = cvpartition(classLabels,'holdout',holdoutRatio);
dataTrain = featureMat(holdoutCVP.training,:);
grpTrain = classLabels(holdoutCVP.training,:);
dataTest = featureMat(holdoutCVP.test,:);
grpTest = classLabels(holdoutCVP.test,:);

%% Sequential feature selection
% following instructions from https://uk.mathworks.com/help/stats/examples/selecting-features-for-classifying-high-dimensional-data.html

if performSequentialFeatureSelection
    
    % ecdf plot of t-test or anova test values
    grpTrainInd = grp2idx(grpTrain);
    if strcmp(classVar,'wormNum')
        % t-test for two density class labels: 5 vs. 40
        dataTrainG1 = dataTrain(grpTrainInd==1,:);
        dataTrainG2 = dataTrain(grpTrainInd==2,:);
        [~,p] = ttest2(dataTrainG1,dataTrainG2,'Vartype','unequal');
    elseif strcmp(classVar,'strain_name')
        % one-way anova for multiple class labels: strain names
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
    
    % sort features
    [~,featureIdxSortbyP] = sort(p,2);
    
    % define function for sfs. Linear discriminant analysis currently but should probably change for favoured model.
    % TODO:  % check out fitcensemble for subspace discriminant model.
    classf = @(xtrain,ytrain,xtest,ytest) ...
        sum(~strcmp(ytest,classify(xtest,xtrain,ytrain,'linear')));
    
    % sequential feature selection with k-fold cross-validation
    kfoldCVP = cvpartition(grpTrain,'kfold',crossVal_k);
    if ~isnan(n_sortedFeatsInput)
        featureIdxSortbyP = featureIdxSortbyP(1:n_sortedFeatsInput);
    end
    [fsLocal,historyCV] = sequentialfs(classf,dataTrain(:,featureIdxSortbyP),grpTrain,'cv',kfoldCVP,'nfeatures',n_topFeatsOutput);
    featCols = featureIdxSortbyP(fsLocal);
    keyFeats = featureTable.Properties.VariableNames(featCols)'
    
    if plotSFS
        % evaluate the performance of the selected model with the features
        testMCELocal = crossval(classf,featureMat(:,featureIdxSortbyP(fsLocal)),classLabels,'partition',...
            holdoutCVP)/holdoutCVP.TestSize;
        
        % plot of the cross-validation MCE as a function of the number of top features
        figure; plot(historyCV.Crit,'o');
        xlabel('Number of Features');
        ylabel('CV MCE');
        title('Forward Sequential Feature Selection with cross-validation');
        
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
    disp('Please use GUI to train model with trainFeatureTable. Use five-fold cross validation, select the best model, and export to workspace.')
    classificationLearner
    
    % Wait until trained model has been selected and saved
    while exist('trainedModel')~=1
        pause
    end
    
    % For classifying strains, linear discriminant model and sometimes Ensemble subspace discrimnant model work the best.
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