clear
close all

%% script takes saved features data (i.e. area, speed),
% 1) calculates metrics of the features data for cluster, mw, and sw blob categories (i.e. median of each replicate), and
% 2) generates box plots using user-defined metric (i.e. median), sorted by median or 90th percentile values (i.e. median of median, or 90th percentile of median), and
% 3) performs one-way ANOVA of the means and include p-values in the figure title, and
% 4) exports median and 90th percentile feature median values (and removes DA609) for CeNDR GWAS mapping

%% set parameter
feature = 'perdurance'; % specify feature as string. 'area','compactness','perimeter','quirkiness','solidity','speed','perdurance'
strainSet = 'all'; %'all'
metric = 4; % 4 by default to use median from each replicate. %1 = min, 2 = 10th prc, 3 = 25th prc, 4 = median, 5 = 75th prc, 6 = 90th prc, 7 = max
saveResults = true;
export4Mapping = true;

%% set feature-specific parameter
if strcmp(feature,'area')
    unit = ' (mm^2)'; % micron squared
    applySwNormalisation = true;
elseif strcmp(feature,'perdurance')
    unit = ' (seconds elapsed)';
    applySwNormalisation = false; % this feature should not be normalised against single worm
elseif strcmp(feature,'speed')
    unit = ' (\mum/s)';
    applySwNormalisation = true;
elseif strcmp(feature,'compactness')
    unit = '';
    applySwNormalisation = true;
elseif strcmp(feature,'quirkiness')
    unit = '';
    applySwNormalisation = true;
elseif strcmp(feature,'solidity')
    unit = '';
    applySwNormalisation = true;
elseif strcmp(feature,'perimeter')
    unit = ' (mm)';
    applySwNormalisation = true;
end

%% prep work
metricText = [{'min'};{'10th_percentile'};{'25th_percentile'};{'median'};{'75th_percentile'};{'90th_percentile'};{'max'}];
addpath('auxiliary/')
load(['results/' feature '_all.mat'])

% create empty figures
cluster_boxplotFig_sortMedian = figure;
mw_boxplotFig_sortMedian = figure;
sw_boxplotFig_sortMedian = figure;
cluster_boxplotFig_sort90prc = figure;
mw_boxplotFig_sort90prc = figure;
sw_boxplotFig_sort90prc = figure;
if applySwNormalisation
    clusterNorm_boxplotFig_sortMedian = figure;
    mwNorm_boxplotFig_sortMedian = figure;
    clusterNorm_boxplotFig_sort90prc = figure;
    mwNorm_boxplotFig_sort90prc = figure;
end

% set eps export options
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',80,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',3);

%% calculate metric values if necessary (i.e. median, 90th percentile, etc.) for each strain and replicate
try % see if metric value has already been calculated and saved for loading
    load(['results/featMetricVals/' feature '_all.mat'])
    % calculate total number of reps for all strains
    numStrains = numel(fieldnames(cluster_feature));
    totalRepCtr = 0; % initialise total rep counter
    for strainCtr = 1:numStrains
        strainNames = fieldnames(cluster_feature);
        numReps = numel(cluster_feature.(strainNames{strainCtr}));
        totalRepCtr = totalRepCtr + numReps; % update total rep counter
    end
catch % calculate metric value only if it doesn't already exist
    % go through each strain
    numStrains = numel(fieldnames(cluster_feature));
    totalRepCtr = 0; % initialise total rep counter
    for strainCtr = 1:numStrains
        strainNames = fieldnames(cluster_feature);
        strain = strainNames{strainCtr}; % get strain name
        numReps = numel(cluster_feature.(strain));
        totalRepCtr = totalRepCtr + numReps; % update total rep counter
        % preallocate variable to hold calculated feature metric values for that strain/rep
        cluster_featVals.(strain) = NaN(numReps,7); % 7 columns to contain [min,10th,25th,50th,75th,90th percentiles,max] values
        mw_featVals.(strain) = NaN(numReps,7);
        sw_featVals.(strain) = NaN(numReps,7);
        if applySwNormalisation
            clusterNorm_featVals.(strain) = NaN(numReps,7); % 7 columns to contain [min,10th,25th,50th,75th,90th percentiles,max] values
            mwNorm_featVals.(strain) = NaN(numReps,7);
        end
        % go through each replicate
        for repCtr = 1:numReps
            cluster_repFeat = cluster_feature.(strain){repCtr};
            mw_repFeat = mw_feature.(strain){repCtr};
            sw_repFeat = sw_feature.(strain){repCtr};
            % apply swNormalisation as specified
            if applySwNormalisation
                swMedian = nanmedian(sw_repFeat);
                clusterNorm_repFeat = cluster_repFeat/swMedian;
                mwNorm_repFeat = mw_repFeat/swMedian;
            end
            % populate variable with values
            if ~isempty(cluster_repFeat)
                cluster_featVals.(strain)(repCtr,1) = min(cluster_repFeat);
                cluster_featVals.(strain)(repCtr,2) = prctile(cluster_repFeat,10);
                cluster_featVals.(strain)(repCtr,3) = prctile(cluster_repFeat,25);
                cluster_featVals.(strain)(repCtr,4) = nanmedian(cluster_repFeat);
                cluster_featVals.(strain)(repCtr,5) = prctile(cluster_repFeat,75);
                cluster_featVals.(strain)(repCtr,6) = prctile(cluster_repFeat,90);
                cluster_featVals.(strain)(repCtr,7) = max(cluster_repFeat);
            end
            if ~isempty(mw_repFeat)
                mw_featVals.(strain)(repCtr,1) = min(mw_repFeat);
                mw_featVals.(strain)(repCtr,2) = prctile(mw_repFeat,10);
                mw_featVals.(strain)(repCtr,3) = prctile(mw_repFeat,25);
                mw_featVals.(strain)(repCtr,4) = nanmedian(mw_repFeat);
                mw_featVals.(strain)(repCtr,5) = prctile(mw_repFeat,75);
                mw_featVals.(strain)(repCtr,6) = prctile(mw_repFeat,90);
                mw_featVals.(strain)(repCtr,7) = max(mw_repFeat);
            end
            if ~isempty(sw_repFeat)
                sw_featVals.(strain)(repCtr,1) = min(sw_repFeat);
                sw_featVals.(strain)(repCtr,2) = prctile(sw_repFeat,10);
                sw_featVals.(strain)(repCtr,3) = prctile(sw_repFeat,25);
                sw_featVals.(strain)(repCtr,4) = nanmedian(sw_repFeat);
                sw_featVals.(strain)(repCtr,5) = prctile(sw_repFeat,75);
                sw_featVals.(strain)(repCtr,6) = prctile(sw_repFeat,90);
                sw_featVals.(strain)(repCtr,7) = max(sw_repFeat);
            end
            if applySwNormalisation
                if ~isempty(clusterNorm_repFeat)
                    clusterNorm_featVals.(strain)(repCtr,1) = min(clusterNorm_repFeat);
                    clusterNorm_featVals.(strain)(repCtr,2) = prctile(clusterNorm_repFeat,10);
                    clusterNorm_featVals.(strain)(repCtr,3) = prctile(clusterNorm_repFeat,25);
                    clusterNorm_featVals.(strain)(repCtr,4) = nanmedian(clusterNorm_repFeat);
                    clusterNorm_featVals.(strain)(repCtr,5) = prctile(clusterNorm_repFeat,75);
                    clusterNorm_featVals.(strain)(repCtr,6) = prctile(clusterNorm_repFeat,90);
                    clusterNorm_featVals.(strain)(repCtr,7) = max(clusterNorm_repFeat);
                end
                if ~isempty(mwNorm_repFeat)
                    mwNorm_featVals.(strain)(repCtr,1) = min(mwNorm_repFeat);
                    mwNorm_featVals.(strain)(repCtr,2) = prctile(mwNorm_repFeat,10);
                    mwNorm_featVals.(strain)(repCtr,3) = prctile(mwNorm_repFeat,25);
                    mwNorm_featVals.(strain)(repCtr,4) = nanmedian(mwNorm_repFeat);
                    mwNorm_featVals.(strain)(repCtr,5) = prctile(mwNorm_repFeat,75);
                    mwNorm_featVals.(strain)(repCtr,6) = prctile(mwNorm_repFeat,90);
                    mwNorm_featVals.(strain)(repCtr,7) = max(mwNorm_repFeat);
                end
            end
        end
    end
    % save calculated rep metrics
    if saveResults
        if applySwNormalisation
            save(['results/featMetricVals/' feature '_all.mat'],'cluster_featVals','mw_featVals','sw_featVals','clusterNorm_featVals','mwNorm_featVals');
        else
            save(['results/featMetricVals/' feature '_all.mat'],'cluster_featVals','mw_featVals','sw_featVals');
        end
    end
end

%% generate box plot inputs for user-defined feature metric
% initialise
rowCtr = 1;
boxplotStrainNames = cell(totalRepCtr,1); % grouping variable
cluster_boxplotValues = NaN(totalRepCtr,1); % values to plot
cluster_boxPlotSort = cell(numStrains,3); % variable to sort grouping order by
mw_boxplotValues = NaN(totalRepCtr,1);
mw_boxPlotSort = cell(numStrains,3);
sw_boxplotValues = NaN(totalRepCtr,1);
sw_boxPlotSort = cell(numStrains,3);
if applySwNormalisation
    clusterNorm_boxplotValues = NaN(totalRepCtr,1);
    clusterNorm_boxPlotSort = cell(numStrains,3);
    mwNorm_boxplotValues = NaN(totalRepCtr,1);
    mwNorm_boxPlotSort = cell(numStrains,3);
end

% fill boxplot values and labels
for strainCtr = 1:numel(strainNames)
    strain = strainNames{strainCtr};
    numReps = numel(cluster_feature.(strain));
    % get appropriate user-defined metric values (i.e. median value of each replicate)
    cluster_featMetricVals = cluster_featVals.(strain)(:,metric);
    mw_featMetricVals = mw_featVals.(strain)(:,metric);
    sw_featMetricVals = sw_featVals.(strain)(:,metric);
    if applySwNormalisation
        clusterNorm_featMetricVals = clusterNorm_featVals.(strain)(:,metric);
        mwNorm_featMetricVals = mwNorm_featVals.(strain)(:,metric);
    end
    % populate boxplot variables
    boxplotStrainNames(rowCtr:rowCtr+numReps-1) = {strain};
    cluster_boxplotValues(rowCtr:rowCtr+numReps-1) = cluster_featMetricVals;
    mw_boxplotValues(rowCtr:rowCtr+numReps-1) = mw_featMetricVals;
    sw_boxplotValues(rowCtr:rowCtr+numReps-1) = sw_featMetricVals;
    if applySwNormalisation
        clusterNorm_boxplotValues(rowCtr:rowCtr+numReps-1) = clusterNorm_featMetricVals;
        mwNorm_boxplotValues(rowCtr:rowCtr+numReps-1) = mwNorm_featMetricVals;
    end
    % get boxplot median and 90th percentile for sorting
    cluster_boxPlotSort{strainCtr,1} = strain;
    cluster_boxPlotSort{strainCtr,2} = nanmedian(cluster_featMetricVals); % get the median of the feature
    cluster_boxPlotSort{strainCtr,3} = prctile(cluster_featMetricVals,90); % get the 90th percentile of the feature
    mw_boxPlotSort{strainCtr,1} = strain;
    mw_boxPlotSort{strainCtr,2} = nanmedian(mw_featMetricVals); % get the median of the feature
    mw_boxPlotSort{strainCtr,3} = prctile(mw_featMetricVals,90); % get the 90th percentile of the feature
    sw_boxPlotSort{strainCtr,1} = strain;
    sw_boxPlotSort{strainCtr,2} = nanmedian(sw_featMetricVals); % get the median of the feature
    sw_boxPlotSort{strainCtr,3} = prctile(sw_featMetricVals,90); % get the 90th percentile of the feature
    if applySwNormalisation
        clusterNorm_boxPlotSort{strainCtr,1} = strain;
        clusterNorm_boxPlotSort{strainCtr,2} = nanmedian(clusterNorm_featMetricVals); % get the median of the feature
        clusterNorm_boxPlotSort{strainCtr,3} = prctile(clusterNorm_featMetricVals,90); % get the 90th percentile of the feature
        mwNorm_boxPlotSort{strainCtr,1} = strain;
        mwNorm_boxPlotSort{strainCtr,2} = nanmedian(mwNorm_featMetricVals); % get the median of the feature
        mwNorm_boxPlotSort{strainCtr,3} = prctile(mwNorm_featMetricVals,90); % get the 90th percentile of the feature
    end
    % update row counter
    rowCtr = rowCtr + numReps;
end

% sort strain order based on boxplot median or 90th percentile values
% by median
cluster_boxPlotSort = sortrows(cluster_boxPlotSort,2);
mw_boxPlotSort = sortrows(mw_boxPlotSort,2);
sw_boxPlotSort = sortrows(sw_boxPlotSort,2);
cluster_groupOrder_median = cluster_boxPlotSort(:,1);
mw_groupOrder_median = mw_boxPlotSort(:,1);
sw_groupOrder_median = sw_boxPlotSort(:,1);
if applySwNormalisation
    clusterNorm_boxPlotSort = sortrows(clusterNorm_boxPlotSort,2);
    mwNorm_boxPlotSort = sortrows(mwNorm_boxPlotSort,2);
    clusterNorm_groupOrder_median = clusterNorm_boxPlotSort(:,1);
    mwNorm_groupOrder_median = mwNorm_boxPlotSort(:,1);
end
% by 90th percentile
cluster_boxPlotSort = sortrows(cluster_boxPlotSort,3);
mw_boxPlotSort = sortrows(mw_boxPlotSort,3);
sw_boxPlotSort = sortrows(sw_boxPlotSort,3);
cluster_groupOrder_90prc = cluster_boxPlotSort(:,1);
mw_groupOrder_90prc = mw_boxPlotSort(:,1);
sw_groupOrder_90prc = sw_boxPlotSort(:,1);
if applySwNormalisation
    clusterNorm_boxPlotSort = sortrows(clusterNorm_boxPlotSort,3);
    mwNorm_boxPlotSort = sortrows(mwNorm_boxPlotSort,3);
    clusterNorm_groupOrder_90prc = clusterNorm_boxPlotSort(:,1);
    mwNorm_groupOrder_90prc = mwNorm_boxPlotSort(:,1);
end

%% perform ANOVA (to include p-values in figure titles)
clusterANOVA = anova1(cluster_boxplotValues,boxplotStrainNames,'off');
mwANOVA = anova1(mw_boxplotValues,boxplotStrainNames,'off');
swANOVA = anova1(sw_boxplotValues,boxplotStrainNames,'off');
if applySwNormalisation
    clusterNormANOVA = anova1(clusterNorm_boxplotValues,boxplotStrainNames,'off');
    mwNormANOVA = anova1(mwNorm_boxplotValues,boxplotStrainNames,'off');
end

%% plot and format
yLabel = ([feature  unit]);
isitHighlighted = strcmp(boxplotStrainNames, 'N2') + 2*strcmp(boxplotStrainNames, 'DA609') + 3*strcmp(boxplotStrainNames,'CB4856');

% boxplots sorted by median
% cluster boxplot
set(0,'CurrentFigure',cluster_boxplotFig_sortMedian)
boxplot(cluster_boxplotValues,boxplotStrainNames,'GroupOrder',cluster_groupOrder_median,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
set(0,'CurrentFigure',cluster_boxplotFig_sortMedian)
xlabel('stains')
ylabel(['cluster ' yLabel])
title ({['cluster ' feature ', sorted by median'] ['ANOVA p = ' num2str(clusterANOVA)]})
% mw boxplot
set(0,'CurrentFigure',mw_boxplotFig_sortMedian)
boxplot(mw_boxplotValues,boxplotStrainNames,'GroupOrder',mw_groupOrder_median,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
xlabel('stains')
set(0,'CurrentFigure',mw_boxplotFig_sortMedian)
ylabel(['multiworm ' yLabel])
title ({['multiworm ' feature ', sorted by median'] ['ANOVA p = ' num2str(mwANOVA)]})
% sw boxplot
set(0,'CurrentFigure',sw_boxplotFig_sortMedian)
boxplot(sw_boxplotValues,boxplotStrainNames,'GroupOrder',sw_groupOrder_median,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
xlabel('stains')
ylabel(['singleworm ' yLabel])
title ({['singleworm ' feature ', sorted by median'] ['ANOVA p = ' num2str(swANOVA)]})
if applySwNormalisation
    yLabelNorm = (['normalised ' feature]);
    % clusterNorm boxplot
    set(0,'CurrentFigure',clusterNorm_boxplotFig_sortMedian)
    boxplot(clusterNorm_boxplotValues,boxplotStrainNames,'GroupOrder',clusterNorm_groupOrder_median,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
    set(0,'CurrentFigure',clusterNorm_boxplotFig_sortMedian)
    xlabel('stains')
    ylabel(['cluster ' yLabelNorm])
    title ({['normalised cluster ' feature ', sorted by median'] ['ANOVA p = ' num2str(clusterNormANOVA)]})
    % mwNorm boxplot
    set(0,'CurrentFigure',mwNorm_boxplotFig_sortMedian)
    boxplot(mwNorm_boxplotValues,boxplotStrainNames,'GroupOrder',mwNorm_groupOrder_median,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
    xlabel('stains')
    ylabel(['multiworm ' yLabelNorm])
    title ({['normalised multiworm ' feature ', sorted by median'] ['ANOVA p = ' num2str(mwNormANOVA)]})
end

% boxplots sorted by 90th percentile
% cluster boxplot
set(0,'CurrentFigure',cluster_boxplotFig_sort90prc)
boxplot(cluster_boxplotValues,boxplotStrainNames,'GroupOrder',cluster_groupOrder_90prc,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
set(0,'CurrentFigure',cluster_boxplotFig_sort90prc)
xlabel('stains')
ylabel(['cluster ' yLabel])
title ({['cluster ' feature ', sorted by 90th percentile'] ['ANOVA p = ' num2str(clusterANOVA)]})
% mw boxplot
set(0,'CurrentFigure',mw_boxplotFig_sort90prc)
boxplot(mw_boxplotValues,boxplotStrainNames,'GroupOrder',mw_groupOrder_90prc,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
xlabel('stains')
ylabel(['multiworm ' yLabel])
title ({['multiworm ' feature ', sorted by 90th percentile'] ['ANOVA p = ' num2str(mwANOVA)]})
% sw boxplot
set(0,'CurrentFigure',sw_boxplotFig_sort90prc)
boxplot(sw_boxplotValues,boxplotStrainNames,'GroupOrder',sw_groupOrder_90prc,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
xlabel('stains')
ylabel(['singleworm ' yLabel])
title ({['singleworm ' feature ', sorted by 90th percentile'] ['ANOVA p = ' num2str(swANOVA)]})
if applySwNormalisation
    % clusterNorm boxplot
    set(0,'CurrentFigure',clusterNorm_boxplotFig_sort90prc)
    boxplot(clusterNorm_boxplotValues,boxplotStrainNames,'GroupOrder',clusterNorm_groupOrder_90prc,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
    set(0,'CurrentFigure',clusterNorm_boxplotFig_sort90prc)
    xlabel('stains')
    ylabel(['cluster ' yLabelNorm])
    title ({['normalised cluster ' feature ', sorted by 90th percentile'] ['ANOVA p = ' num2str(clusterNormANOVA)]})
    % mwNorm boxplot
    set(0,'CurrentFigure',mwNorm_boxplotFig_sort90prc)
    boxplot(mwNorm_boxplotValues,boxplotStrainNames,'GroupOrder',mwNorm_groupOrder_90prc,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
    xlabel('stains')
    ylabel(['multiworm ' yLabelNorm])
    title ({['normalised multiworm ' feature ', sorted by 90th percentile'] ['ANOVA p = ' num2str(mwNormANOVA)]})
end

%% save figures
if saveResults
    % figure names
    cluster_figurename_byMedian = ['figures/boxplots/' feature '_' metricText{metric} '_cluster_sortMedian'];
    mw_figurename_byMedian = ['figures/boxplots/' feature '_' metricText{metric} '_mw_sortMedian'];
    sw_figurename_byMedian = ['figures/boxplots/' feature '_' metricText{metric} '_sw_sortMedian'];
    cluster_figurename_by90prc = ['figures/boxplots/' feature '_' metricText{metric} '_cluster_sort90prc'];
    mw_figurename_by90prc = ['figures/boxplots/' feature '_' metricText{metric} '_mw_sort90prc'];
    sw_figurename_by90prc = ['figures/boxplots/' feature '_' metricText{metric} '_sw_sort90prc'];
    % append "_ns" to figure names if ANOVA result is not significant
    if clusterANOVA>=0.05
        cluster_figurename_byMedian = [cluster_figurename_byMedian '_ns'];
        cluster_figurename_by90prc = [cluster_figurename_by90prc '_ns'];
    end
    if mwANOVA>=0.05
        mw_figurename_byMedian = [mw_figurename_byMedian '_ns'];
        mw_figurename_by90prc = [mw_figurename_by90prc '_ns'];
    end
    if swANOVA>=0.05
        sw_figurename_byMedian = [sw_figurename_byMedian '_ns'];
        sw_figurename_by90prc = [sw_figurename_by90prc '_ns'];
    end
    % export figures
    exportfig(cluster_boxplotFig_sortMedian,[cluster_figurename_byMedian '.eps'],exportOptions)
    exportfig(mw_boxplotFig_sortMedian,[mw_figurename_byMedian '.eps'],exportOptions)
    exportfig(sw_boxplotFig_sortMedian,[sw_figurename_byMedian '.eps'],exportOptions)
    exportfig(cluster_boxplotFig_sort90prc,[cluster_figurename_by90prc '.eps'],exportOptions)
    exportfig(mw_boxplotFig_sort90prc,[mw_figurename_by90prc '.eps'],exportOptions)
    exportfig(sw_boxplotFig_sort90prc,[sw_figurename_by90prc '.eps'],exportOptions)
    if applySwNormalisation
        % figure names
        clusterNorm_figurename_byMedian = ['figures/boxplots/' feature 'Norm_' metricText{metric} '_cluster_sortMedian'];
        mwNorm_figurename_byMedian = ['figures/boxplots/' feature 'Norm_' metricText{metric} '_mw_sortMedian'];
        clusterNorm_figurename_by90prc = ['figures/boxplots/' feature 'Norm_' metricText{metric} '_cluster_sort90prc'];
        mwNorm_figurename_by90prc = ['figures/boxplots/' feature 'Norm_' metricText{metric} '_mw_sort90prc'];
        % append "_ns" to figure names if ANOVA result is not significant
        if clusterNormANOVA>=0.05
            clusterNorm_figurename_byMedian = [clusterNorm_figurename_byMedian '_ns'];
            clusterNorm_figurename_by90prc = [clusterNorm_figurename_by90prc '_ns'];
        end
        if mwNormANOVA>=0.05
            mwNorm_figurename_byMedian = [mwNorm_figurename_byMedian '_ns'];
            mwNorm_figurename_by90prc = [mwNorm_figurename_by90prc '_ns'];
        end
        % export figures
        exportfig(clusterNorm_boxplotFig_sortMedian,[clusterNorm_figurename_byMedian '.eps'],exportOptions)
        exportfig(mwNorm_boxplotFig_sortMedian,[mwNorm_figurename_byMedian '.eps'],exportOptions)
        exportfig(clusterNorm_boxplotFig_sort90prc,[clusterNorm_figurename_by90prc '.eps'],exportOptions)
        exportfig(mwNorm_boxplotFig_sort90prc,[mwNorm_figurename_by90prc '.eps'],exportOptions)
    end
end
    
%% export mapping variables
if export4Mapping
    % save median values
    mappingFileName = ['results/mapping/' feature '_' metricText{metric} '.mat'];
    if applySwNormalisation
        save(mappingFileName,'cluster_boxPlotSort','mw_boxPlotSort','sw_boxPlotSort','clusterNorm_boxPlotSort','mwNorm_boxPlotSort')
    else
        save(mappingFileName,'cluster_boxPlotSort','mw_boxPlotSort','sw_boxPlotSort')
    end
    % remove DA609 from mapping
    cluster_removeIdx = find(strcmp(cluster_boxPlotSort(:,1),'DA609'));
    clusterNoDA_boxPlotSort = vertcat(cluster_boxPlotSort(1:cluster_removeIdx-1,:),cluster_boxPlotSort(cluster_removeIdx+1:end,:));
    mw_removeIdx = find(strcmp(mw_boxPlotSort(:,1),'DA609'));
    mwNoDA_boxPlotSort = vertcat(mw_boxPlotSort(1:mw_removeIdx-1,:),mw_boxPlotSort(mw_removeIdx+1:end,:));
    sw_removeIdx = find(strcmp(sw_boxPlotSort(:,1),'DA609'));
    swNoDA_boxPlotSort = vertcat(sw_boxPlotSort(1:sw_removeIdx-1,:),sw_boxPlotSort(sw_removeIdx+1:end,:));
    if applySwNormalisation
        clusterNorm_removeIdx = find(strcmp(clusterNorm_boxPlotSort(:,1),'DA609'));
        clusterNormNoDA_boxPlotSort = vertcat(clusterNorm_boxPlotSort(1:clusterNorm_removeIdx-1,:),clusterNorm_boxPlotSort(clusterNorm_removeIdx+1:end,:));
        mwNorm_removeIdx = find(strcmp(mwNorm_boxPlotSort(:,1),'DA609'));
        mwNormNoDA_boxPlotSort = vertcat(mwNorm_boxPlotSort(1:mwNorm_removeIdx-1,:),mwNorm_boxPlotSort(mwNorm_removeIdx+1:end,:));
    end
    % export strain median and 90th percentile values of the user-defined metric (i.e. median) for mapping
    cluster_mappingFileName = ['results/mapping/' feature '_' metricText{metric} '_cluster.txt'];
    mw_mappingFileName = ['results/mapping/' feature '_' metricText{metric} '_mw.txt'];
    sw_mappingFileName = ['results/mapping/' feature '_' metricText{metric} '_sw.txt'];
    dlmcell(cluster_mappingFileName,clusterNoDA_boxPlotSort);
    dlmcell(mw_mappingFileName,mwNoDA_boxPlotSort);
    dlmcell(sw_mappingFileName,swNoDA_boxPlotSort);
    if applySwNormalisation
        clusterNorm_mappingFileName = ['results/mapping/' feature '_' metricText{metric} '_clusterNorm.txt'];
        mwNorm_mappingFileName = ['results/mapping/' feature '_' metricText{metric} '_mwNorm.txt'];
        dlmcell(clusterNorm_mappingFileName,clusterNormNoDA_boxPlotSort);
        dlmcell(mwNorm_mappingFileName,mwNormNoDA_boxPlotSort);
    end
end