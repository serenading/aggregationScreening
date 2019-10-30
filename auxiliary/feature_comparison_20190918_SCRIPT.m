% script to load a feature matrix generated by Tierpsy Tracker summarising
% each plate in an experiment and to do some basic comparisons (clustering
% and feature-by-feature t-tests to compare the samples to a control.


% for this analysis, select the cultivation temperature to analyse
cultTempStr = 'Grown at 20';

% load the feature matrix, corresponding filenames, and metadata
tierpsyFeatureTable = readtable('./Results/features_summary_tierpsy_plate_20190911_214248.csv');
tierpsyFileTable = readtable('./Results/filenames_summary_tierpsy_plate_20190911_214248.csv', 'Format', '%d%s%s');
metadataTable = readtable('./AuxiliaryFiles/metadata.xlsx');

% join the Tierpsy tables to match filenames with file_id. Required in case 
% features were not extracted for any files.
combinedTierpsyTable = outerjoin(tierpsyFileTable, tierpsyFeatureTable, ...
    'MergeKeys', true);

% get just the filenames from the full path in the tables
[~, fileNamesMetadata] = cellfun(@fileparts, metadataTable.filename, 'UniformOutput', false);
metadataTable.filename = fileNamesMetadata;
[~, fileNamesTierpsy] = cellfun(@fileparts, combinedTierpsyTable.file_name, 'UniformOutput', false);
combinedTierpsyTable.file_name = strrep(fileNamesTierpsy, '_featuresN', '');

% rename Tierpsy output to match metadata output
combinedTierpsyTable.Properties.VariableNames{'file_name'} = 'filename';

% finally, join tables to get strain names for each set of features
featureTable = outerjoin(metadataTable, combinedTierpsyTable, ...
    'MergeKeys', true);

% we will analyse the worms cultivated at different temperatures separately
% Select the temperature data here
rowsToKeep = contains(featureTable.comments, cultTempStr);
featureTable = featureTable(rowsToKeep, :);

% get worm strain names and combine with cultivation temperatures
wormNames = featureTable.strain;
uniqueNames = unique(wormNames);

% load the set of 256 features selected using based on classification
% accuracy on a set of mutant strains
filename256_features = ['/Users/abrown/Andre/projects/' ...
    '2018-tierpsy-features/' ...
    'top256_tierpsy_no_blob_no_eigen_only_abs_no_norm.csv'];
top256_all = readtable(filename256_features);
top256 = top256_all.Var2(2:end); % take just one set of 256, drop header

% shorten variable names (done in Excel for imported feature table)
featNames = strrep(top256, 'relative', 'rel');
featNames = strrep(featNames, 'velocity', 'vel');
featNames = strrep(featNames, 'angular', 'ang');

% advantage of indexing by feature name in tables is outweighed by
% complication of most function calls so proceed with feature matrix,
% wormNames, and feature indices
featMat = featureTable{:, featNames};

% impute nan values to global mean

% introduce check for number of imputed values

featMeans = nanmean(featMat);
for ii = 1:size(featMat, 2)
    nanInds = isnan(featMat(:, ii));
    featMat(nanInds, ii) = featMeans(ii);
end

% z-normalise feature matrix
featMatNorm = normalize(featMat);

% make a matrix of mean values for each strain for all features
groupMeanMat = NaN(numel(uniqueNames), size(featMatNorm, 2));
groupNames = grpstats(featMatNorm(:, 1), wormNames, 'gname');
for ii = 1:size(featMatNorm, 2)
    groupMeanMat(:, ii) = grpstats(featMatNorm(:, ii), wormNames, 'mean');
end


% % cluster at the level of strains to check for reproducibility
% clg = clustergram(groupMeanMat, 'RowLabels', groupNames, ...
%     'ColumnLabels', featNames, 'Linkage', 'complete');
% cmap = cbrewer('div', 'RdBu', 50); % RdYlBu
% colormap(cmap(end:-1:1, :))
% set(clg, 'Colormap', cmap(end:-1:1, :)) % 'DisplayRange', 2




% compare all of the MOTT strains to eachother and to each of the OW
% strains

% get the indices of the two groups
mottNames = unique(wormNames( contains(wormNames, 'MOTT') ));
owNames = unique(wormNames( contains(wormNames, 'OW') ));

% initialise
pValMatMott = NaN(numel(mottNames), numel(mottNames), size(featMat, 2));

% loop over unique strains
for ii = 1:numel(mottNames)-1
    disp(ii/(numel(mottNames)-1))
    
    % get the indices of the current strain
    iiNameInds = contains(wormNames, mottNames{ii});
    
    % loop over other strains
    for kk = ii+1:numel(mottNames)
        
        % get the indices of the current strain
        kkNameInds = contains(wormNames, mottNames{kk});
        
        % loop over features
        for jj = 1:size(featMat, 2)
            [~, pVal] = ttest2(featMat(iiNameInds, jj), ...
                featMat(kkNameInds, jj));
            pValMatMott(ii, kk, jj) = pVal;
        end
    end
end


% initialise
pValMatOw = NaN(numel(mottNames), numel(owNames), size(featMat, 2));

% loop over unique strains
for ii = 1:numel(mottNames)
    disp(ii/(numel(mottNames)))
    
    % get the indices of the current strain
    iiNameInds = contains(wormNames, mottNames{ii});
    
    % loop over other strains
    for kk = 1:numel(owNames)
        
        % get the indices of the current strain
        kkNameInds = contains(wormNames, owNames{kk});
        
        % loop over features
        for jj = 1:size(featMat, 2)
            [~, pVal] = ttest2(featMat(iiNameInds, jj), ...
                featMat(kkNameInds, jj));
            pValMatOw(ii, kk, jj) = pVal;
        end
    end
end

% correct for multiple comparisons across all comparisons
% compare fdr_bh output to scikit-learn

combinedPVals = [pValMatMott(:); pValMatOw(:)];
[~, pCrit, bhFDR2] = fdr_bh(combinedPVals(~isnan(combinedPVals)), 0.05 , 'dep');

% get the features with significant differences for each comparison set
hitIndsMott = squeeze(any(any(pValMatMott < pCrit, 1), 2));
hitIndsOw = squeeze(any(any(pValMatOw < pCrit, 1), 2));

% list the group names to set the plot order
plotGroupOrder = [owNames; mottNames];

% make some boxplots
for ii = 1:size(featMat, 2)
    disp(ii/size(featMat, 2))
    
    if ~hitIndsMott(ii) && ~hitIndsOw(ii)
        continue
    end
    
    % check for setting related to visibility to increase speed
    figure('visible', 'off')
    CategoricalScatterplot(featMat(:, ii), wormNames, ...
        'GroupOrder', plotGroupOrder)
    
    % add bars to indicate any significant differences
    
    % find OW strain hits for current feature
    [a, b] = ...
        ind2sub(size(pValMatOw(:, :, ii)), find(pValMatOw(:, :, ii) < pCrit));
    hitNames1 = mottNames(a);
    hitNames2 = owNames(b);
    
    % use group order to find start and end coordinates for each bar
    for nn = 1:size(hitNames1, 1)
        % get coordinates for current bar
        startCoord = find(strcmp(plotGroupOrder, hitNames1(nn)));
        endCoord = find(strcmp(plotGroupOrder, hitNames2(nn)));

        % plot bar
        plot([startCoord, endCoord], ...
            [max(featMat(:, ii)), max(featMat(:, ii))] ...
            + nn * 0.05 * range(featMat(:, ii)), ...
            'LineWidth', 2)
    end
    
    ylabel(strrep(featNames(ii), '_', ' '))
    
%     xlim([0.5, 4.5])
    
    set(gca, 'FontSize', 14, 'XTickLabelRotation', 45)
%     set(gca, 'PlotBoxAspectRatio', [0.6, 1, 1])
    saveas(gcf, ['./plots/hit-plots-' cultTempStr '/feature-plot-' num2str(ii) '.png'])
    close
end







