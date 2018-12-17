% Script for multiworm tracking data of worms treated with antipsychotics


% what is the minimum trajectory length in frames?
minLength = 3000;

% import the experiment description file to get drug names and doses
% 
%first read as text data
fileID = fopen('/Volumes/behavgenom_archive$/Adam/screening/antipsychotics/ExtraFiles/exp_list.csv');
textData = textscan(fileID, '%s %s',  'Delimiter', ',');
fclose(fileID);
% 
% % split the text data into parts
fileNames = textData{1};
drugNames = textData{2};

% set the root directory
directory = '/Volumes/behavgenom_archive$/Adam/screening/antipsychotics/Results/Antipsychotics_Agar_Screening_070617/';

% get a list of trajectory files
[fileList, ~] = dirSearch(directory, '_features.hdf5');

% get the drug names from the filelist
drugNames = cell(size(fileList));
for ii = 1:numel(fileList)
    [~, bn, ~] = fileparts(fileList{ii});
    dd = strsplit(bn, '_');
    drugNames(ii) = strcat(dd(3), '_', dd(4));
end

% get unique class names
uniqueNames = unique(drugNames);

% impport the mean features data into a trajectory-level feature matrix
featFlag = true(726, 1);

% get featNames (which don't vary) from first file
[~, ~, featNames] = featStruct2Mat(fileList{1}, 0, featFlag);

% loop over files to grow total feature matrix and mean feature matrix
% (averaged across all long trajectories on a given plate)
featMatTotal = [];
featMatMean = NaN(numel(fileList), numel(featNames));
namesTotal = {};
count = 1;
for ii = 1:numel(fileList)
    disp(ii/numel(fileList))
    
    % import feature data
    [featMat, ~, ~] = ...
        featStruct2Mat(fileList{ii}, minLength, featFlag);
    
    % some feature files are empty, skip these
    if isempty(featMat)
        continue
    end
    
    % add to featMatTotal
    featMatTotal = [featMatTotal; featMat];
    
    % expand list of total names, including the doses
    for jj = 1:size(featMat, 1)
        namesTotal{count, 1} = drugNames{ii};
        count = count + 1;
    end
    
    % average and add to featMatMean
    featMatMean(ii, :) = nanmean(featMat, 1);
end

% drop features with too many NaNs
nanCols = false(size(featMat, 2), 1);
for ii = 1:size(featMat, 2)
    if sum(isnan(featMat(:, ii))) > 0.5 * size(featMat, 1)
        nanCols(ii) = true;
    end
end
featMatMean(:, nanCols) = [];
featMatTotal(:, nanCols) = [];
featNames(nanCols) = [];

% impute NaN values to the column means
meanVals = nanmean(featMatTotal);
featMatTotalNoNaN = featMatTotal;
for ii = 1:size(featMatTotal, 1)
    % get the NaN values in the current row
    nanVals = isnan(featMatTotal(ii, :));
    
    % replace the nanVals
    featMatTotalNoNaN(ii, nanVals) = meanVals(nanVals);
end

meanVals = nanmean(featMatTotal);
featMatMeanNoNaN = featMatMean;
for ii = 1:size(featMatMean, 1)
    % get the NaN values in the current row
    nanVals = isnan(featMatMean(ii, :));
    
    % replace the nanVals
    featMatMeanNoNaN(ii, nanVals) = meanVals(nanVals);
end

% z-normalise
featMatTotalNorm = bsxfun(@minus, featMatTotalNoNaN, mean(featMatTotalNoNaN));
featMatTotalNorm = bsxfun(@rdivide, featMatTotalNorm, std(featMatTotalNoNaN));

% do PCA on the total feature matrix
[pc, score, ~, ~, explained] = pca(featMatTotalNorm);

% plot data on first two PCs
figure
for ii = 1:numel(uniqueNames)
    rowInds = strcmp(namesTotal, uniqueNames{ii});
    if strcmp(uniqueNames{ii}, 'DMSO')
        plot(score(rowInds, 1), score(rowInds, 2), ...
            '.', 'MarkerSize', 12, 'Color', 'r')
        hold on
    else
        plot(score(rowInds, 1), score(rowInds, 2), ...
            '.', 'MarkerSize', 12, 'Color', rand(1, 3))
        hold on
    end
end
hold off


% cluster data using first few PCs as features
clg = clustergram(score(:, 1:15), 'RowLabels', namesTotal, ...
    'Linkage', 'complete');
cmap = cbrewer('div', 'RdBu', 50); % RdYlBu
colormap(cmap(end:-1:1, :))
set(clg, 'Colormap', cmap(end:-1:1, :), 'DisplayRange', 2)



% % also make a biplot of the PCA
figure;
biplot(pc(:, 1:2), 'Scores', score(:, 1:2), 'VarLabels', featNames)


% z-normalise featMatMean
featMatNorm = bsxfun(@minus, featMatMeanNoNaN, mean(featMatMeanNoNaN));
featMatNorm = bsxfun(@rdivide, featMatNorm, std(featMatMeanNoNaN));


% compare each compound to the control data in each
% feature using rank sum tests or t-tests
controlMeans = featMatMean(strcmp(drugNames, 'DMSO_14.08'), :);
pVals = NaN(numel(uniqueNames), size(featMatMean, 2));
for ii = 1:numel(uniqueNames)
    % skip iteration for the control data
    if strcmp(uniqueNames{ii}, 'DMSO_14.08')
        continue
    end
    
    % get the indices for the current comparison
    currentInds = strcmp(drugNames, uniqueNames{ii});
    
    for jj = 1:size(featMatMean, 2)
        % also do rank sum tests
%         p = ranksum(featMatMean(currentInds, jj), controlMeans(:, jj));
        [~, p] = ttest2(featMatMean(currentInds, jj), controlMeans(:, jj));

        % get the p-value
        pVals(ii, jj) = p;
    end
end

% Correct for multiple comparisons using Benjamini procedure.
[~, pCrit, bhFDR] = fdr_bh(pVals, 0.05 , 'dep');



% do PCA on the mean feature matrix
[~, score] = pca(featMatNorm);

% plot data on first two PCs
figure
for ii = 1:numel(uniqueNames)
    rowInds = strcmp(drugNames, uniqueNames{ii});
    if strcmp(uniqueNames{ii}, 'DMSO_14.08')
        plot(score(rowInds, 1), score(rowInds, 2), '.', 'MarkerSize', 25, 'color', [0 0 0])
    elseif ~isempty(strfind(uniqueNames{ii}, 'Clozapine_100'))
        plot(score(rowInds, 1), score(rowInds, 2), '.', 'MarkerSize', 25, 'color', [1 0 0])
    else
        plot(score(rowInds, 1), score(rowInds, 2), ...
            '.', 'MarkerSize', 12, 'Color', rand(1, 3))
    end
    hold on
end
hold off

% get a preliminary ranking of the features using F-statistic
[rankF, fStats] = rankFeaturesFStat(featMatTotalNorm, namesTotal);


% cluster videos and features
clg = clustergram(featMatNorm(:, rankF(1:100)), 'RowLabels', drugNames, ...
    'ColumnLabels', featNames(rankF(1:100)), 'Linkage', 'complete');
cmap = cbrewer('div', 'RdBu', 50); % RdYlBu
colormap(cmap(end:-1:1, :))
set(clg, 'Colormap', cmap(end:-1:1, :), 'DisplayRange', 2)

% % cluster data using first few PCs as features
% clg = clustergram(score(:, 1:45), 'RowLabels', drugNames, ...
%     'Linkage', 'complete');
% cmap = cbrewer('div', 'RdBu', 50); % RdYlBu
% colormap(cmap(end:-1:1, :))
% set(clg, 'Colormap', cmap(end:-1:1, :), 'DisplayRange', 2)



% also make a version of feature matrix averaged across all videos in each
% condition
featMatMean2 = NaN(numel(uniqueNames), size(featMatTotal, 2));
for ii = 1:numel(uniqueNames)
    % get the indices of the current video in featMatTotal
    currentInds = strcmp(namesTotal, uniqueNames{ii});
    
    % average features across worms in the current video
    featMatMean2(ii, :) = nanmean(featMatTotalNoNaN(currentInds, :));
end

% z-normalise featMatTotal
featMatNorm2 = bsxfun(@minus, featMatMean2, nanmean(featMatMean2));
featMatNorm2 = bsxfun(@rdivide, featMatNorm2, nanstd(featMatMean2));

% some conditions have no good data and so can still be nan
nanRows = isnan(featMatNorm2(:, 1));

% do PCA on the mean feature matrix
[~, score] = pca(featMatNorm2);

% plot data on first two PCs
figure
for ii = 1:numel(uniqueNames)
    plot(score(ii, 1), score(ii, 2), ...
        '.', 'MarkerSize', 12, 'Color', rand(1, 3))
    hold on
end
hold off

% cluster conditions and features
clg = clustergram(featMatNorm2(~nanRows, :), 'RowLabels', uniqueNames(~nanRows), ...
    'ColumnLabels', featNames, 'Linkage', 'complete');
set(clg, 'Colormap', cmap(end:-1:1, :), 'DisplayRange', 2)



% t-sne
no_dims = 2;
initial_dims = 12;
[~, ~, numLabels] = unique(namesTotal);
ydata = tsne(featMatTotalNorm, numLabels, no_dims, initial_dims);

% plot the points from each compound
figure
for ii = 1:numel(uniqueNames)
    currentInds = strcmp(namesTotal, uniqueNames{ii});
    if strcmp(uniqueNames{ii}, 'DMSO_14.08')
        plot(ydata(currentInds, 1), ydata(currentInds, 2), '.', 'MarkerSize', 25, 'color', [0 0 0])
    elseif ~isempty(strfind(uniqueNames{ii}, 'Clozapine_100'))
        plot(ydata(currentInds, 1), ydata(currentInds, 2), '.', 'MarkerSize', 25, 'color', [1 0 0])
    else
        plot(ydata(currentInds, 1), ydata(currentInds, 2), 'o', 'MarkerSize', 5)
    end
    %     plot3(ydata(currentInds, 1), ydata(currentInds, 2), ydata(currentInds, 3), '.', 'MarkerSize', 50)
    hold on
end
hold off



% % partition the data into a training set for feature selection and a test
% % set for subsequent testing
% clozapineNames = namesTotal;
% notClozapine100Inds = find(~strcmp(namesTotal, 'Clozapine_100'));
% for ii = 1:numel(notClozapine100Inds)
%     clozapineNames{notClozapine100Inds(ii)} = 'Not Clozapine';
% end
% 
% 
% 
% holdoutCVP = cvpartition(clozapineNames, 'holdout', 0.2);
% 
% % make a training and test feature matrix
% featMatTrain = featMatTotalNorm(holdoutCVP.training, :);
% classesTrain = namesTotal(holdoutCVP.training);
% featMatTest = featMatTotalNorm(holdoutCVP.test, :);
% classesTest = namesTotal(holdoutCVP.test);
% 
% 
% % because featMat is relatively small, we can use sequential feature
% % selection on the whole feature matrix
% tenfoldCVP = cvpartition(classesTrain, 'kfold', 10);
% options = statset('Display', 'iter');
% [~, historyCV] = ...
%     sequentialfs(@svmError, featMatTrain(:, rankF(1:20)), classesTrain, ...
%     'cv', tenfoldCVP, 'Nf', 10, 'options', options);
% 
% 
% % plot the cross-validation error rate as a function of the number of
% % features included
% figure
% plot(historyCV.Crit)
% xlim([0, 11])
% 
% % use plot to choose feature number
% selectedFeatures = historyCV.In(4, :);
% 
% % now, use these selected features to train an SVM classifier and to
% % predict the labels of the held-out test set
% t = templateSVM('Standardize', 1, 'KernelFunction', 'linear');
% model = fitcecoc(featMatTrain(:, selectedFeatures), classesTrain, 'Learners', t);
% predictedLabels = predict(model, featMatTest(:, selectedFeatures));
% 
% % calculate the total error
% error = sum(~strcmp(predictedLabels, classesTest)) / size(featMatTest, 1);
% disp(['error on test set: ' num2str(error) '.'])
% 
% 
% % make a scatter plot of the data
% firstFeature = logical(historyCV.In(1,:));
% secondFeature = logical(historyCV.In(2,:) - firstFeature);
% figure
% gscatter(featMatTrain(:, firstFeature), featMatTrain(:, secondFeature), ...
%     classesTrain);
% figure
% gscatter(featMatTest(:, firstFeature), featMatTest(:, secondFeature), ...
%     classesTest);
% 
% 
% 
% 
% 
% 
