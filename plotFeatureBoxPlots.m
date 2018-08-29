clear
close all

%% script takes saved features data (i.e. area, speed),
%% 1) calculates metrics of the features data for cluster, mw, and sw blob categories, and
%% 2) generates sorted box plots using user-defined metric (i.e. median, 90th percentile)

%% set parameter
feature = 'quirkiness';
metric = 4; %1 = min, 2 = 10th prc, 3 = 25th prc, 4 = median, 5 = 75th prc, 6 = 90th prc, 7 = max
saveResults = true;
export4Mapping = true;

%% set feature-specific parameter
if strcmp(feature,'area')
    unit = '\mum^2'; % micron squared
    applySwNormalisation = true;
elseif strcmp(feature,'perdurance')
    unit = 'frames elapsed at 25fps';
    applySwNormalisation = false; % this feature should not be normalised against single worm
elseif strcmp(feature,'speed')
    unit = '\mum/s';
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
    unit = '\mum';
    applySwNormalisation = true;
end

%% initialise
metricText = [{'min'};{'10th_percentile'};{'25th_percentile'};{'median'};{'75th_percentile'};{'90th_percentile'};{'max'}];
load(['results/' feature '.mat'])
cluster_boxplotFig = figure;
mw_boxplotFig = figure;
sw_boxplotFig = figure;
if applySwNormalisation
    clusterNorm_boxplotFig = figure;
    mwNorm_boxplotFig = figure;
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
try % calculate metric value only if it doesn't already exist, otherwise just load
    load(['results/featMetricVals/' featureNameSave '.mat'])
    % calculate total number of reps for all strains
    numStrains = numel(fieldnames(cluster_featureNorm));
    totalRepCtr = 0; % initialise total rep counter
    for strainCtr = 1:numStrains
        strainNames = fieldnames(cluster_featureNorm);
        numReps = numel(cluster_featureNorm.(strainNames{strainCtr}));
        totalRepCtr = totalRepCtr + numReps; % update total rep counter
    end
catch 
    % go through each strain
    numStrains = numel(fieldnames(cluster_featureNorm));
    totalRepCtr = 0; % initialise total rep counter
    for strainCtr = 1:numStrains
        strainNames = fieldnames(cluster_featureNorm);
        strain = strainNames{strainCtr}; % get strain name
        numReps = numel(cluster_featureNorm.(strain));
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
            cluster_repFeat = cluster_featureNorm.(strain){repCtr};
            mw_repFeat = mw_featureNorm.(strain){repCtr};
            sw_repFeat = sw_featureNorm.(strain){repCtr};
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
        save(['results/featMetricVals/' feature '.mat'],'cluster_featVals','mw_featVals','sw_featVals');
        if applySwNormalisation
            save(['results/featMetricVals/' feature 'Norm.mat'],'clusterNorm_featVals','mwNorm_featVals');
        end
    end
end

%% generate box plot for user-defined feature metric
% initialise
rowCtr = 1;
boxplotStrainNames = cell(totalRepCtr,1); % grouping variable
cluster_boxplotValues = NaN(totalRepCtr,1); % values to plot
cluster_boxPlotSort = cell(numStrains,2); % variable to sort grouping order by
mw_boxplotValues = NaN(totalRepCtr,1);
mw_boxPlotSort = cell(numStrains,2);
sw_boxplotValues = NaN(totalRepCtr,1);
sw_boxPlotSort = cell(numStrains,2);
if applySwNormalisation
    clusterNorm_boxplotValues = NaN(totalRepCtr,1);
    clusterNorm_boxPlotSort = cell(numStrains,2);
    mwNorm_boxplotValues = NaN(totalRepCtr,1);
    mwNorm_boxPlotSort = cell(numStrains,2);
end

% fill boxplot values and labels
for strainCtr = 1:numel(strainNames)
    strain = strainNames{strainCtr};
    numReps = numel(cluster_featureNorm.(strain));
    % get appropriate user-defined metric values
    cluster_featMetricVals = cluster_featVals.(strain)(:,metric);
    mw_featMetricVals = mw_featVals.(strain)(:,metric);
    sw_featMetricVals = sw_featVals.(strain)(:,metric);
    if applySwNormalisation
        clusterNorm_featMetricVals = clusterNorm_featVals.(strain)(:,metric);
        mwNorm_featMetricVals = mw_featVals.(strain)(:,metric);
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
    % get boxplot median for sorting
    cluster_boxPlotSort{strainCtr,1} = strain;
    cluster_boxPlotSort{strainCtr,2} = nanmedian(cluster_featMetricVals); % get the median of the feature
    mw_boxPlotSort{strainCtr,1} = strain;
    mw_boxPlotSort{strainCtr,2} = nanmedian(mw_featMetricVals); % get the median of the feature
    sw_boxPlotSort{strainCtr,1} = strain;
    sw_boxPlotSort{strainCtr,2} = nanmedian(sw_featMetricVals); % get the median of the feature
    if applySwNormalisation
        clusterNorm_boxPlotSort{strainCtr,1} = strain;
        clusterNorm_boxPlotSort{strainCtr,2} = nanmedian(clusterNorm_featMetricVals); % get the median of the feature
        mwNorm_boxPlotSort{strainCtr,1} = strain;
        mwNorm_boxPlotSort{strainCtr,2} = nanmedian(mwNorm_featMetricVals); % get the median of the feature
    end
    % update row counter
    rowCtr = rowCtr + numReps;
end

% sort strain order based on boxplot median values
cluster_boxPlotSort = sortrows(cluster_boxPlotSort,2);
mw_boxPlotSort = sortrows(mw_boxPlotSort,2);
sw_boxPlotSort = sortrows(sw_boxPlotSort,2);
cluster_groupOrder = cluster_boxPlotSort(:,1);
mw_groupOrder = mw_boxPlotSort(:,1);
sw_groupOrder = sw_boxPlotSort(:,1);
if applySwNormalisation
    clusterNorm_boxPlotSort = sortrows(clusterNorm_boxPlotSort,2);
    mwNorm_boxPlotSort = sortrows(mwNorm_boxPlotSort,2);
    clusterNorm_groupOrder = clusterNorm_boxPlotSort(:,1);
    mwNorm_groupOrder = mwNorm_boxPlotSort(:,1);
end

% plot and format
% yaxis label
metricNameYAxis = metricText{metric};
if metric == 2 | metric == 3 | metric == 5 | metric == 6
    metricNameYAxis = strrep(metricNameYAxis,'_','\_'); % modify ylabel display format so the p isn't lower case
end
if strcmp(unit,'')
    yLabel = ([metricNameYAxis ' ' feature]);
else
    yLabel = ([metricNameYAxis ' ' feature ' (' unit ')']);
end
% cluster boxplot
set(0,'CurrentFigure',cluster_boxplotFig)
boxplot(cluster_boxplotValues,boxplotStrainNames,'GroupOrder',cluster_groupOrder,'PlotStyle','compact','LabelVerbosity','all')
set(0,'CurrentFigure',cluster_boxplotFig)
xlabel('stains')
ylabel(['cluster ' yLabel])
% mw boxplot
set(0,'CurrentFigure',mw_boxplotFig)
boxplot(mw_boxplotValues,boxplotStrainNames,'GroupOrder',mw_groupOrder,'PlotStyle','compact','LabelVerbosity','all')
xlabel('stains')
ylabel(['multiworm ' yLabel])
% sw boxplot
set(0,'CurrentFigure',sw_boxplotFig)
boxplot(sw_boxplotValues,boxplotStrainNames,'GroupOrder',sw_groupOrder,'PlotStyle','compact','LabelVerbosity','all')
xlabel('stains')
ylabel(['singleworm ' yLabel])
if applySwNormalisation
    yLabel = (['normalised ' feature ' (' metricNameYAxis ')']);
    % clusterNorm boxplot
    set(0,'CurrentFigure',clusterNorm_boxplotFig)
    boxplot(clusterNorm_boxplotValues,boxplotStrainNames,'GroupOrder',clusterNorm_groupOrder,'PlotStyle','compact','LabelVerbosity','all')
    set(0,'CurrentFigure',clusterNorm_boxplotFig)
    xlabel('stains')
    ylabel(['cluster ' yLabel])
    % mwNorm boxplot
    set(0,'CurrentFigure',mwNorm_boxplotFig)
    boxplot(mwNorm_boxplotValues,boxplotStrainNames,'GroupOrder',mwNorm_groupOrder,'PlotStyle','compact','LabelVerbosity','all')
    xlabel('stains')
    ylabel(['multiworm ' yLabel])
end

% save figures
if saveResults
    cluster_figurename = ['figures/boxplots/' feature '_' metricText{metric} '_cluster'];
    mw_figurename = ['figures/boxplots/' feature '_' metricText{metric} '_mw'];
    sw_figurename = ['figures/boxplots/' feature '_' metricText{metric} '_sw'];
    exportfig(cluster_boxplotFig,[cluster_figurename '.eps'],exportOptions)
    exportfig(mw_boxplotFig,[mw_figurename '.eps'],exportOptions)
    exportfig(sw_boxplotFig,[sw_figurename '.eps'],exportOptions)
    if applySwNormalisation
        clusterNorm_figurename = ['figures/boxplots/' feature 'Norm_' metricText{metric} '_cluster'];
        mwNorm_figurename = ['figures/boxplots/' feature 'Norm_' metricText{metric} '_mw'];
        exportfig(clusterNorm_boxplotFig,[clusterNorm_figurename '.eps'],exportOptions)
        exportfig(mwNorm_boxplotFig,[mwNorm_figurename '.eps'],exportOptions)
    end
end

if export4Mapping
    %% export mapping variables
    % save median values
    mappingFileName = ['results/mapping/' feature '_' metricText{metric} '.mat'];
    save(mappingFileName,'cluster_boxPlotSort','mw_boxPlotSort','sw_boxPlotSort')
    if applySwNormalisation
        mappingFileName = ['results/mapping/' feature '_' metricText{metric} 'Norm.mat'];
        save(mappingFileName,'clusterNorm_boxPlotSort','mwNorm_boxPlotSort')
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
    % export strain median value (i.e. across replicates) of the user-defined metric (i.e. median or 90th percentile)
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