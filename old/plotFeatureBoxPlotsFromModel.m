clear
close all

%% set parameter
feature = 'ac'; % 'ac' or 'pcf'
strainSet = 'all'; %'all'
export4Mapping = false;
saveResults = false;

%% set feature-specific parameter
if strcmp(feature,'ac')
    modelFun = 'exp3'; % 'exp3', 'mexp3' or 'power2'
    resultFilePath = ['/Users/sding/Documents/AggScreening/results/modelFit/ac_all_fitmodel_' modelFun '_ac_sample0.32s_8pixel.mat'];
    if strcmp(modelFun, 'exp3')
        paramList = {'a','b','c','d','e'};
    elseif strcmp(modelFun, 'mexp3')
        paramList = {'a','b','c','d','e','f'};
    elseif strcmp(modelFun, 'power2')
        paramList = {'a','b','c'};
    end
elseif strcmp(feature,'pcf')
    modelFun = 'exp2'; % 'exp1' or 'exp2'
    resultFilePath = ['/Users/sding/Documents/AggScreening/results/modelFit/pcf_all_fitmodel_' modelFun '_pcf_sample10s_1pixel.mat'];
    if strcmp(modelFun, 'exp1')
        paramList = {'a','b'};
    elseif strcmp(modelFun, 'exp2')
        paramList = {'a','b','c','d'};
    end
end

%% prep work
addpath('auxiliary/')
load(resultFilePath)
numParams = numel(paramList);

% create empty figures
boxplotFig = figure;

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
    load(['results/featMetricVals/' feature '_' modelFun 'FitModel_all.mat'])
    % calculate total number of reps for all strains
    numStrains = numel(fieldnames(f));
    totalRepCtr = 0; % initialise total rep counter
    for strainCtr = 1:numStrains
        strainNames = fieldnames(f);
        numReps = numel(f.(strainNames{strainCtr}));
        totalRepCtr = totalRepCtr + numReps; % update total rep counter
    end
catch % calculate metric value only if it doesn't already exist
    numStrains = numel(fieldnames(f));
    totalRepCtr = 0; % initialise total rep counter
    for strainCtr = 1:numStrains
        strainNames = fieldnames(f);
        strain = strainNames{strainCtr}; % get strain name
        numReps = numel(f.(strain));
        totalRepCtr = totalRepCtr + numReps; % update total rep counter
        % preallocate variable to hold calculated feature metric values for that strain/rep
        modelParams_featVals.(strain) = NaN(numReps,numParams);
        % go through each replicate
        for repCtr = 1:numReps
            if ~isempty(f.(strain){repCtr})
                for paramCtr = 1:numParams
                    modelParams_featVals.(strain)(repCtr,paramCtr) = f.(strain){repCtr}.(paramList{paramCtr});
                end
            end
        end
    end
    % save calculated rep metrics
    if saveResults
        save(['results/featMetricVals/' feature '_' modelFun 'FitModel_all.mat'],'modelParams_featVals');
    end
end

%% initialise param feature export for mapping
if export4Mapping
    % pre-allocate export matrix for mapping
    mappingFeatValExport = cell(numStrains,numParams+1);
    mappingFeatValExport(:,1) = strainNames;
    headingText = cell(1,numParams+1);
    headingText{1,1} = 'strains';
end

%% generate box plots for model parameters
for paramCtr = 1:numParams
    % initialise
    rowCtr = 1;
    boxplotStrainNames = cell(totalRepCtr,1); % grouping variable
    boxplotValues = NaN(totalRepCtr,1); % values to plot
    boxPlotSort = cell(numStrains,2); % variable to sort grouping order by
    % fill boxplot values and labels
    for strainCtr = 1:numel(strainNames)
        strain = strainNames{strainCtr};
        numReps = numel(f.(strain));
        % get appropriate parameter values (i.e. a,b,or c of each replicate)
        featMetricVals = modelParams_featVals.(strain)(:,paramCtr);
        % populate boxplot variables
        boxplotStrainNames(rowCtr:rowCtr+numReps-1) = {strain};
        boxplotValues(rowCtr:rowCtr+numReps-1) = featMetricVals;
        % get boxplot median for sorting
        boxPlotSort{strainCtr,1} = strain;
        boxPlotSort{strainCtr,2} = nanmedian(featMetricVals); % get the median of the parameter value
        % update row counter
        rowCtr = rowCtr + numReps;
    end
    % sort strain order based on boxplot median values
    boxPlotSort = sortrows(boxPlotSort,2);
    groupOrder = boxPlotSort(:,1);
    %% perform ANOVA (to include p-values in figure titles)
    ANOVA = anova1(boxplotValues,boxplotStrainNames,'off');
    %% plot and format
    yLabel = ([feature ' param ' paramList{paramCtr}]);
    isitHighlighted = strcmp(boxplotStrainNames, 'N2') + 2*strcmp(boxplotStrainNames, 'DA609') + 3*strcmp(boxplotStrainNames,'CB4856');
    % boxplot
    set(0,'CurrentFigure',boxplotFig)
    boxplot(boxplotValues,boxplotStrainNames,'GroupOrder',groupOrder,'ColorGroup',isitHighlighted,'Colors','krbg','PlotStyle','compact','LabelVerbosity','all')
    set(0,'CurrentFigure',boxplotFig)
    xlabel('stains')
    ylabel(yLabel)
    title ({[feature ' param ' paramList{paramCtr}] ['ANOVA p = ' num2str(ANOVA)]})
    
    %% save figures
    if saveResults
        % figure names
        figurename = ['figures/boxplots/' feature '_' modelFun '_param_' paramList{paramCtr}];
        % append "_ns" to figure names if ANOVA result is not significant
            if ANOVA>=0.05
                figurename = [figurename '_ns'];
            end
        % export figures
        exportfig(boxplotFig,[figurename '.eps'],exportOptions)
    end
    %% export mapping variables
    if export4Mapping
        % re-sort values by strain name and horizontally concatenate parameter values
        boxPlotSort = sortrows(boxPlotSort,1);
        mappingFeatValExport(:,paramCtr+1) = boxPlotSort(:,2);
        headingText{paramCtr+1} = ['param_' paramList{paramCtr}];
        if saveResults
            % save median values (including DA609)
            mappingFileName = ['results/mapping/' feature '_' modelFun 'FitModel.mat'];
            save(mappingFileName,'mappingFeatValExport')
        end
    end
end

%% export mapping variables
if export4Mapping
    % remove DA609 from mapping
    removeIdx = find(strcmp(mappingFeatValExport(:,1),'DA609')); % should be 17
    assert(removeIdx == 17)
    mappingFeatValExport = vertcat(mappingFeatValExport(1:removeIdx-1,:),mappingFeatValExport(removeIdx+1:end,:));
    % add heading
    mappingFeatValExport = vertcat(headingText,mappingFeatValExport);
    % export for mapping
    if saveResults
        mappingFeatValExportName = ['results/mapping/' feature '_' modelFun 'FitModel.txt'];
        dlmcell(mappingFeatValExportName,mappingFeatValExport);
    end
end