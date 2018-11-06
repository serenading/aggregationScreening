clear
close all

%% script takes saved features data (i.e. area, speed),
% 1) calculates metrics of the features data for cluster, mw, and sw blob categories, and
% 2) generates box plots using user-defined metric (i.e. median), sorted by median values (i.e. median of median, or median of 90th percentile), and
% 3) performs one-way ANOVA of the means and include p-values in the figure title, and
% 4) exports feature metric values (and removes DA609) for CeNDR GWAS mapping

%% set parameter
feature = 'area'; % specify feature as string. 'area','compactness','perimeter','quirkiness','solidity','speed','perdurance','perduranceSurvival','pcf','hc'
strainSet = 'all'; %'all'
metrics = [4 6]; % [4 6] to use median and 90th prc from each replicate. %1 = min, 2 = 10th prc, 3 = 25th prc, 4 = median, 5 = 75th prc, 6 = 90th prc, 7 = max
saveResults = false;
export4Mapping = false;
metricText = [{'min'};{'10percentile'};{'25percentile'};{'median'};{'75percentile'};{'90percentile'};{'max'}];

%% set feature-specific parameter
if strcmp(feature,'area')
    unit = ' (mm^2)'; % micron squared
    applySwNormalisation = true;
elseif strcmp(feature,'perdurance') | strcmp(feature,'perduranceSurvival')
    unit = ' (seconds elapsed)';
    applySwNormalisation = false; 
    if strcmp(feature,'perduranceSurvival')
        metricText = [{'100%'};{'90%'};{'75%'};{'50%'};{'25%'};{'10%'};{'0%'}];
    end
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
elseif strcmp(feature,'pcf')
    unit = '';
    applySwNormalisation = false;
    load('/Users/sding/Documents/AggScreening/results/pcf_all_sample10s_1pixel.mat')
    if size(pcf.DA609{1},1) == 15
        distBinWidth = 0.1; % in units of mm
        maxDist = 1.5; % in units of mm
        distBins = 0:distBinWidth:maxDist;
    else
        error('unknown binWidth for pcf data')
    end
end

%% prep work
addpath('auxiliary/')
if ~strcmp(feature,'pcf')
    load(['results/' feature '_all.mat'])
end
% rename loaded variables if necessary
if strcmp(feature,'perduranceSurvival')
    cluster_feature = cluster_perdDist;
    mw_feature = mw_perdDist;
    sw_feature = sw_perdDist;
elseif strcmp(feature,'pcf')
    cluster_feature = pcf; % for this feature "cluster/mw/sw" categoes do not matter %%%%%
    mw_feature = pcf; 
    sw_feature = pcf;
end

% create empty figures
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
        if strcmp(feature,'perduranceSurvival')
            numReps = 1; % survival data already pooled across different reps
        end
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
            if strcmp(feature,'perduranceSurvival')
                cluster_repFeat = cluster_feature.(strain);
                mw_repFeat = mw_feature.(strain);
                sw_repFeat = sw_feature.(strain);
            else
                cluster_repFeat = cluster_feature.(strain){repCtr};
                mw_repFeat = mw_feature.(strain){repCtr};
                sw_repFeat = sw_feature.(strain){repCtr};
            end
        end
        % apply swNormalisation as specified
        if applySwNormalisation
            swMedian = nanmedian(sw_repFeat);
            clusterNorm_repFeat = cluster_repFeat/swMedian;
            mwNorm_repFeat = mw_repFeat/swMedian;
        end
        % perform necessary feature calculations/conversions
        if strcmp(feature,'pcf')
            for frameCtr = 1:size(cluster_repFeat,2)
                % concatenate pcf values (value set equal to the median of the containing bin) for each frame
                pcfVal{frameCtr} = []; 
                for binCtr = 1:numel(distBins)-1
                    binVal = median([distBins(binCtr) distBins(binCtr+1)]); %
                    pcfVal{frameCtr} = [pcfVal{frameCtr} ones(1,round(cluster_repFeat(binCtr,frameCtr)/1e4))*binVal];
                end
            end
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
        if strcmp(feature,'perduranceSurvival') % overwrite values for cluster_featVals, mw_featVals, sw_featVals
            if ~isempty(cluster_feature.(strain))
                [ecdfy,ecdfx] = ecdf(cluster_feature.(strain));
                survivalFraction = 1-ecdfy;
                [~,idx] = min(abs(survivalFraction - 1)); % find the index for survival = 1 (100%)
                cluster_featVals.(strain)(repCtr,1) = ecdfx(idx); % find the time elapsed for 100% survival
                [~,idx] = min(abs(survivalFraction - 0.9));
                cluster_featVals.(strain)(repCtr,2) = ecdfx(idx);
                [~,idx] = min(abs(survivalFraction - 0.75));
                cluster_featVals.(strain)(repCtr,3) = ecdfx(idx);
                [~,idx] = min(abs(survivalFraction - 0.5)); % find the index for survival = 0.5 (50%)
                cluster_featVals.(strain)(repCtr,4) = ecdfx(idx); % find the time elapsed for 50% survival
                [~,idx] = min(abs(survivalFraction - 0.25));
                cluster_featVals.(strain)(repCtr,5) = ecdfx(idx);
                [~,idx] = min(abs(survivalFraction - 0.1)); % find the index for survival = 0.1 (10%) 
                cluster_featVals.(strain)(repCtr,6) = ecdfx(idx);% find the time elapsed for 10% survival
                [~,idx] = min(abs(survivalFraction - 0));
                cluster_featVals.(strain)(repCtr,7) = ecdfx(idx);
            end
            if ~isempty(mw_feature.(strain))
                [ecdfy,ecdfx] = ecdf(mw_feature.(strain));
                survivalFraction = 1-ecdfy;
                [~,idx] = min(abs(survivalFraction - 1)); % find the index for survival = 1 (100%)
                mw_featVals.(strain)(repCtr,1) = ecdfx(idx); % find the time elapsed for 100% survival
                [~,idx] = min(abs(survivalFraction - 0.9));
                mw_featVals.(strain)(repCtr,2) = ecdfx(idx);
                [~,idx] = min(abs(survivalFraction - 0.75));
                mw_featVals.(strain)(repCtr,3) = ecdfx(idx);
                [~,idx] = min(abs(survivalFraction - 0.5)); % find the index for survival = 0.5 (50%)
                mw_featVals.(strain)(repCtr,4) = ecdfx(idx); % find the time elapsed for 50% survival
                [~,idx] = min(abs(survivalFraction - 0.25));
                mw_featVals.(strain)(repCtr,5) = ecdfx(idx);
                [~,idx] = min(abs(survivalFraction - 0.1)); % find the index for survival = 0.1 (10%) 
                mw_featVals.(strain)(repCtr,6) = ecdfx(idx);% find the time elapsed for 10% survival
                [~,idx] = min(abs(survivalFraction - 0));
                mw_featVals.(strain)(repCtr,7) = ecdfx(idx);
            end
            if ~isempty(sw_feature.(strain))
                [ecdfy,ecdfx] = ecdf(sw_feature.(strain));
                survivalFraction = 1-ecdfy;
                [~,idx] = min(abs(survivalFraction - 1)); % find the index for survival = 1 (100%)
                sw_featVals.(strain)(repCtr,1) = ecdfx(idx); % find the time elapsed for 100% survival
                [~,idx] = min(abs(survivalFraction - 0.9));
                sw_featVals.(strain)(repCtr,2) = ecdfx(idx);
                [~,idx] = min(abs(survivalFraction - 0.75));
                sw_featVals.(strain)(repCtr,3) = ecdfx(idx);
                [~,idx] = min(abs(survivalFraction - 0.5)); % find the index for survival = 0.5 (50%)
                sw_featVals.(strain)(repCtr,4) = ecdfx(idx); % find the time elapsed for 50% survival
                [~,idx] = min(abs(survivalFraction - 0.25));
                sw_featVals.(strain)(repCtr,5) = ecdfx(idx);
                [~,idx] = min(abs(survivalFraction - 0.1)); % find the index for survival = 0.1 (10%) 
                sw_featVals.(strain)(repCtr,6) = ecdfx(idx);% find the time elapsed for 10% survival
                [~,idx] = min(abs(survivalFraction - 0));
                sw_featVals.(strain)(repCtr,7) = ecdfx(idx);
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


%% generate box plot inputs for user-defined feature metrics
% initialise
for metricCtr = 1:numel(metrics)
    metric = metrics(metricCtr);
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
    
    %% perform ANOVA (to include p-values in figure titles)
    clusterANOVA = anova1(cluster_boxplotValues,boxplotStrainNames,'off');
    mwANOVA = anova1(mw_boxplotValues,boxplotStrainNames,'off');
    swANOVA = anova1(sw_boxplotValues,boxplotStrainNames,'off');
    if applySwNormalisation
        clusterNormANOVA = anova1(clusterNorm_boxplotValues,boxplotStrainNames,'off');
        mwNormANOVA = anova1(mwNorm_boxplotValues,boxplotStrainNames,'off');
    end
    
    %% plot and format
    yLabel = ([feature ' ' metricText{metric} ' ' unit]);
    isitHighlighted = strcmp(boxplotStrainNames, 'N2') + 2*strcmp(boxplotStrainNames, 'DA609') + 3*strcmp(boxplotStrainNames,'CB4856');
    
    % cluster boxplot
    set(0,'CurrentFigure',cluster_boxplotFig)
    boxplot(cluster_boxplotValues,boxplotStrainNames,'GroupOrder',cluster_groupOrder,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
    set(0,'CurrentFigure',cluster_boxplotFig)
    xlabel('stains')
    ylabel(['cluster ' yLabel])
    if strcmp(feature,'perduranceSurvival')
        title (['cluster ' feature ' ' metricText{metric}])
    else
        title ({['cluster ' feature ' ' metricText{metric}] ['ANOVA p = ' num2str(clusterANOVA)]})
    end
    % mw boxplot
    set(0,'CurrentFigure',mw_boxplotFig)
    boxplot(mw_boxplotValues,boxplotStrainNames,'GroupOrder',mw_groupOrder,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
    set(0,'CurrentFigure',mw_boxplotFig)
    xlabel('stains')
    ylabel(['multiworm ' yLabel])
    if strcmp(feature,'perduranceSurvival')
        title (['multiworm ' feature ' ' metricText{metric}])
    else
        title ({['multiworm ' feature ' ' metricText{metric}] ['ANOVA p = ' num2str(mwANOVA)]})
    end
    % sw boxplot
    set(0,'CurrentFigure',sw_boxplotFig)
    boxplot(sw_boxplotValues,boxplotStrainNames,'GroupOrder',sw_groupOrder,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
    set(0,'CurrentFigure',sw_boxplotFig)
    xlabel('stains')
    ylabel(['singleworm ' yLabel])
    if strcmp(feature,'perduranceSurvival')
        title (['singleworm ' feature ' ' metricText{metric}])
    else
        title ({['singleworm ' feature ' ' metricText{metric}] ['ANOVA p = ' num2str(swANOVA)]})
    end
    if applySwNormalisation
        yLabelNorm = (['normalised ' feature]);
        % clusterNorm boxplot
        set(0,'CurrentFigure',clusterNorm_boxplotFig)
        boxplot(clusterNorm_boxplotValues,boxplotStrainNames,'GroupOrder',clusterNorm_groupOrder,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
        set(0,'CurrentFigure',clusterNorm_boxplotFig)
        xlabel('stains')
        ylabel(['cluster ' yLabelNorm])
        title ({['normalised cluster ' feature ' ' metricText{metric}] ['ANOVA p = ' num2str(clusterNormANOVA)]})
        % mwNorm boxplot
        set(0,'CurrentFigure',mwNorm_boxplotFig)
        boxplot(mwNorm_boxplotValues,boxplotStrainNames,'GroupOrder',mwNorm_groupOrder,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
        set(0,'CurrentFigure',mwNorm_boxplotFig)
        xlabel('stains')
        ylabel(['multiworm ' yLabelNorm])
        title ({['normalised multiworm ' feature ' ' metricText{metric}] ['ANOVA p = ' num2str(mwNormANOVA)]})
    end
    
    %% save figures
    if saveResults
        % figure names
        cluster_figurename = ['figures/boxplots/' feature '_' metricText{metric} '_cluster'];
        mw_figurename = ['figures/boxplots/' feature '_' metricText{metric} '_mw'];
        sw_figurename = ['figures/boxplots/' feature '_' metricText{metric} '_sw'];
        % append "_ns" to figure names if ANOVA result is not significant
        if ~strcmp(feature,'perduranceSurvival')
            if clusterANOVA>=0.05
                cluster_figurename = [cluster_figurename '_ns'];
            end
            if mwANOVA>=0.05
                mw_figurename = [mw_figurename '_ns'];
            end
            if swANOVA>=0.05
                sw_figurename = [sw_figurename '_ns'];
            end
        end
        % export figures
        exportfig(cluster_boxplotFig,[cluster_figurename '.eps'],exportOptions)
        exportfig(mw_boxplotFig,[mw_figurename '.eps'],exportOptions)
        exportfig(sw_boxplotFig,[sw_figurename '.eps'],exportOptions)
        if applySwNormalisation
            % figure names
            clusterNorm_figurename = ['figures/boxplots/' feature 'Norm_' metricText{metric} '_cluster'];
            mwNorm_figurename = ['figures/boxplots/' feature 'Norm_' metricText{metric} '_mw'];
            % append "_ns" to figure names if ANOVA result is not significant
            if clusterNormANOVA>=0.05
                clusterNorm_figurename = [clusterNorm_figurename '_ns'];
            end
            if mwNormANOVA>=0.05
                mwNorm_figurename = [mwNorm_figurename '_ns'];
            end
            % export figures
            exportfig(clusterNorm_boxplotFig,[clusterNorm_figurename '.eps'],exportOptions)
            exportfig(mwNorm_boxplotFig,[mwNorm_figurename '.eps'],exportOptions)
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
        % re-sort values by strain name and horizontally concatenate cluster/mw/sw/clusterNorm/mwNorm values
        cluster_boxPlotSort = sortrows(cluster_boxPlotSort,1);
        mw_boxPlotSort = sortrows(mw_boxPlotSort,1);
        sw_boxPlotSort = sortrows(sw_boxPlotSort,1);
        if applySwNormalisation
            clusterNorm_boxPlotSort = sortrows(clusterNorm_boxPlotSort,1);
            mwNorm_boxPlotSort = sortrows(mwNorm_boxPlotSort,1);
        end
        mappingFeatValExport = [cluster_boxPlotSort mw_boxPlotSort(:,2) sw_boxPlotSort(:,2)];
        if applySwNormalisation
            mappingFeatValExport = [mappingFeatValExport clusterNorm_boxPlotSort(:,2) mwNorm_boxPlotSort(:,2)];
        end
        % remove DA609 from mapping
        removeIdx = find(strcmp(mappingFeatValExport(:,1),'DA609')); % should be 17
        mappingFeatValExport = vertcat(mappingFeatValExport(1:removeIdx-1,:),mappingFeatValExport(removeIdx+1:end,:));
        % add heading
        headingText = {'strain', 'cluster','mw','sw'};
        if applySwNormalisation
            headingText = {'strain', 'cluster','mw','sw','clusterNorm','mwNorm'};
        end
        mappingFeatValExport = vertcat(headingText,mappingFeatValExport);
        % export for mapping
        mappingFeatValExportName = ['results/mapping/' feature '_' metricText{metric} '.txt'];
        dlmcell(mappingFeatValExportName,mappingFeatValExport);
    end
end