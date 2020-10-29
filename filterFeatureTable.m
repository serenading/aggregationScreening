%% Function filters features table by specified strain and feature requirements before classification tasks.

%% Inputs:
% featureTable: file-by-feature featureTable
% classVar: String or cell array of strings containing the name of the variable to classify for. Each string must be a variable field of the featureTable.
% n_nonFeatVar: Scalar specifying the first number of columns that do not contain features to use for classification.
% strains2keep: Cell array containing strains to keep for analysis. Use all strains if left empty. {'all'} or {'divergent'} or {'controls'} or {'strain1', 'strain2'}.
% strains2drop: Cell array containing strains to drop from analysis.
% feats2keep: Cell array containing features to use for analysis. Use all features if left empty. % {'Tierpsy_256'} or  {'feat1','feat2'}
% feats2drop: Cell array containing features to drop from analysis. Partial name of feature allowed. 

%% Outputs: 
% featureTable: file-by-feature featureTable with non-feature metadata variables removed.
% classLabels: file x 1 class labels to feed into the classification task. If multiple multiple variables are specied for classVar then classLabels is a struct where each field contains the labels for that variable.

function [featureTable, classLabels,filenames] = filterFeatureTable(featureTable,classVar,n_nonFeatVar,strains2keep,strains2drop,feats2keep,feats2drop)

%% process rows (observations) first

% retain strains as specified
if isempty(strains2keep) | strcmp(strains2keep,'all')
    load('strainsList/all.mat','strains');
    strains2keep = strains;
elseif strcmp(strains2keep,'divergent')
    load('strainsList/divergent.mat','strains');
    strains2keep = strains;
elseif strcmp(strains2keep,'controls')
    load('strainsList/controls.mat','strains');
    strains2keep = strains;
elseif strcmp(strains2keep,'swept_liberal')
    load('strainsList/swept_strains_liberal.mat','sweptStrains')
    strains2keep = sweptStrains;
elseif strcmp(strains2keep,'nonSwept_liberal')
    load('strainsList/swept_strains_liberal.mat','nonSweptStrains')
    strains2keep = nonSweptStrains;
elseif strcmp(strains2keep,'swept_conservative')
    load('strainsList/swept_strains_conservative.mat','sweptStrains')
    strains2keep = sweptStrains;
elseif strcmp(strains2keep,'nonSwept_conservative')
    load('strainsList/swept_strains_conservative.mat','nonSweptStrains')
    strains2keep = nonSweptStrains;
end
strainLogInd = ismember(featureTable.strain_name,strains2keep);
featureTable = featureTable(strainLogInd,:);
% drop strains as specified
[featureTable,~] = dropStrains(featureTable,strains2drop);
% get strain classification labels
if ischar(classVar)
    classLabels = featureTable.(classVar); 
else
    for varCtr = 1:numel(classVar)
        var = classVar{varCtr};
        classLabels.(var) = featureTable.(var);
    end
end

%% process columns (features) second

% trim table down to retain necessary info
featureTable = featureTable(:,n_nonFeatVar+1:end); % get full features matrix
% retain features as specified
if ~isempty(feats2keep)
    if strcmp(feats2keep,'Tierpsy_256')
        feats2keep = table2cell(readtable('strainsList/Tierpsy_256_short.csv','PreserveVariableNames',true,'ReadVariableNames',false));
    end
    featureTable = featureTable(:,feats2keep);
end
% drop features as specified
[featureTable,~] = dropFeats(featureTable,feats2drop);