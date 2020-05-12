clear
close all

% choose which feature to plot
feat2use = {'path_curvature_midbody_norm_abs_50th'};

%% Extract value for strain names and the desired feature
T = readtable('/Users/sding/Desktop/AggregationScreening/fiveWorm/fiveWormFeaturesTable_20191024_122847.csv');
strainNameColIdx = find(strcmp(T.Properties.VariableNames,'strain_name'));
featColIdx = find(strcmp(T.Properties.VariableNames,feat2use));
assert(numel(featColIdx)==1,'More than one feature with the specified name is found');
boxplotTable = T(:,[strainNameColIdx,featColIdx]);
boxplotVal = table2array(T(:,featColIdx));
boxplotGrp = table2cell(T(:,strainNameColIdx));

%% Sort grouping variable by mean
% extract group stats
[means,strains,n] = grpstats(boxplotVal,boxplotGrp,{'mean','gname','numel'});
% concatenate extracted stats
stats = table(means,strains,n);
% sort stats by means
stats_sorted = sortrows(stats);

% go through each strain generate new grouping variables 
boxplotGrp_sorted = NaN(size(boxplotGrp));
colour_sorted = zeros(numel(boxplotGrp_sorted),3); % [0 0 0] gives black default colour

for grpCtr = 1:numel(stats_sorted.strains)
    strain = stats_sorted.strains(grpCtr);
    strainLogInd = strcmp(boxplotGrp,strain);
    assert(nnz(strainLogInd) == stats_sorted.n(grpCtr));
    boxplotGrp_sorted(strainLogInd) = grpCtr;
    % generate colour labels of interest
    if strcmp(strain,'N2')
        colour_sorted(grpCtr,:) = [1 0 0]; % red
    elseif strcmp(strain,'CB4856')
        colour_sorted(grpCtr,:) = [0 1 1]; % cyan
    elseif strcmp(strain,'DA609')
        colour_sorted(grpCtr,:) = [0 0 1]; % blue
    end
end
assert(nnz(isnan(boxplotGrp_sorted))==0,'Some boxplotGrp_sorted values are NaN')

%% Boxplot
H = boxplot(boxplotVal,boxplotGrp_sorted,...
    'PlotStyle','compact','BoxStyle','filled',...
    'Labels',stats_sorted.strains,'LabelOrientation','inline','Colors',colour_sorted);
title(feat2use,'Interpreter','none')
savefig(['/Users/sding/Desktop/AggregationScreening/fiveWorm/' char(feat2use) '.fig'])