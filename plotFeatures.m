clear
close all

%% script extracts blob features for each of the strains of interest, normalise them against corresponding single worm blob feature within each replicate, and plots a histogram for each strain

%% set parameters
% set analysis parameters
strainSet = 'divergent'; % 'controls','divergent','all'
feature = 'perdurance'; % specify feature as string. 'area','compactness','perimeter','quirkiness','solidity','speed','perdurance'
clusterArea = 4; % 4 by default
phaseRestrict = true; % phaseRestrict cuts out the first 15 min of each video
saveResults = true;
maxNumReplicates = 10;
applySwNormalisation = true; % by default, normalise mw and cluster data against sw from the same movie

% set plotting parameters
histogramNormalisation = 'pdf'; % 'pdf' by default. 'count' an option
yscale = 'linear'; % linear by default, 'log' an option
unit = ''; % nothing by default

% set feature-specific parameters
if strcmp(feature,'speed')
    histogramXLim = [0 5];
    unit = 'microns/s';
elseif strcmp(feature,'area') | strcmp(feature,'perimeter')
    yscale = 'log';
elseif strcmp(feature,'perdurance')
    yscale = 'log'; 
    applySwNormalisation = false;
    unit = 'frames elapsed (at 25fps)';
end

% set eps export options
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',30,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',25,...
    'LineWidth',3);

%% prep work
addpath('auxiliary/')
% load the strain names included in the specified strainSet
load(['strainsList/' strainSet '.mat'])
% get list of file names for each strain
[strainFileList,~,~] = getFileList(strains);
% generate colormap for plotting each strain
colorMap = distinguishable_colors(length(strains));
% create empty figures
if ~applySwNormalisation
    sw_featureFig = figure; hold on
    sw_featurePooledFig = figure; hold on
end
mw_featureFig = figure; hold on
mw_featurePooledFig = figure; hold on
cluster_featureFig = figure; hold on
cluster_featurePooledFig = figure; hold on
% create legend variable to hold strain name and experiment n numbers
legends = cell(size(strains));

%% go through each strain
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    filenames = strainFileList.([strain 'List_40']);
    % if there are many files, then subsample recordings without replacement
    if length(filenames)>maxNumReplicates
        fileInd = datasample(1:length(filenames),maxNumReplicates,'Replace',false);
    else
        fileInd = 1:length(filenames);
    end
    
    %% go through each recording
    for fileCtr = 1:length(fileInd)
        
        %% load data
        filename = filenames{fileInd(fileCtr)};
        trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        % features = h5read(strrep(filename,'skeletons','featuresN'),'/features_timeseries/');
        
        %% filter data
        % find single, multi, and cluster worms
        [singleWormLogInd,multiWormLogInd,clusterLogInd] = findClusters(trajData,blobFeats,clusterArea);
        % generate logical index for phase restriction (cuts out the first 15 min of each 45 min video)
        if phaseRestrict
            frameLogInd = trajData.frame_number>frameRate*60*15;
            singleWormLogInd = singleWormLogInd & frameLogInd;
            multiWormLogInd = multiWormLogInd & frameLogInd;
            clusterLogInd = clusterLogInd & frameLogInd;
        end
        % filter out manually labeled bad entries
        if isfield(trajData,'worm_label')
            singleWormLogInd = singleWormLogInd & trajData.worm_label~=3;
            multiWormLogInd = multiWormLogInd & trajData.worm_label~=3;
            clusterLogInd = clusterLogInd & trajData.worm_label~=3;
        end
        % filter out any blobs that do not persist for more than 1 second
        tempBlobLogInd = findTempBlobs(trajData,blobFeats,frameRate);
        singleWormLogInd = singleWormLogInd & ~tempBlobLogInd;
        multiWormLogInd = multiWormLogInd & ~tempBlobLogInd;
        clusterLogInd = clusterLogInd & ~tempBlobLogInd;
        
        %% analyse features
        % insert function to calculate features if absent
        if strcmp(feature,'speed')
            blobFeats.speed = calculateBlobSpeed(trajData, blobFeats,frameRate);
        elseif strcmp(feature,'perdurance')
            blobFeats.perdurance = trajData.worm_index_joined; % load variable, calculate feature later
        end
        % read features
        sw_feature.(strain){fileCtr} = blobFeats.(feature)(singleWormLogInd);
        mw_feature.(strain){fileCtr} = blobFeats.(feature)(multiWormLogInd);
        cluster_feature.(strain){fileCtr} = blobFeats.(feature)(clusterLogInd);
        % further feature calculation if necessary
        if strcmp(feature,'perdurance')
            [sw_feature.(strain){fileCtr},mw_feature.(strain){fileCtr},cluster_feature.(strain){fileCtr}] = ...
                calculatePerdurance(sw_feature,mw_feature,cluster_feature,strain,fileCtr);
        end
        % normalise features from this movie with sw features from this movie (if appropriate); store value for threshold box plot later
        if ~applySwNormalisation % just rename variable for consistency even if no normalisation is applied
            mw_featureNorm.(strain){fileCtr} = mw_feature.(strain){fileCtr};
            cluster_featureNorm.(strain){fileCtr} = cluster_feature.(strain){fileCtr};
            sw_featureNorm.(strain){fileCtr} = sw_feature.(strain){fileCtr};
        else % normalise using single worm data
            mw_featureNorm.(strain){fileCtr} = mw_feature.(strain){fileCtr}/nanmedian(sw_feature.(strain){fileCtr});
            cluster_featureNorm.(strain){fileCtr} = cluster_feature.(strain){fileCtr}/nanmedian(sw_feature.(strain){fileCtr});
        end
        % update strain n number for figure legend
        legends{strainCtr} = [strain ', n=' num2str(length(mw_featureNorm.(strain)))];
        
        %% plot individual experiments
        set(0,'CurrentFigure',mw_featureFig)
        histogram(mw_featureNorm.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
        set(0,'CurrentFigure',cluster_featureFig)
        histogram(cluster_featureNorm.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
        if ~applySwNormalisation
            set(0,'CurrentFigure',sw_featureFig)
            histogram(sw_featureNorm.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
        end
    end
    
    %% combine the replicates and plot
    mw_featureNormPooled.(strain) = mw_featureNorm.(strain){1};
    cluster_featureNormPooled.(strain) = cluster_featureNorm.(strain){1};
    sw_featureNormPooled.(strain) = sw_featureNorm.(strain){1};
    for fileCtr = 2:length(mw_featureNorm.(strain))
        mw_featureNormPooled.(strain) = vertcat(mw_featureNormPooled.(strain),mw_featureNorm.(strain){fileCtr});
        cluster_featureNormPooled.(strain) = vertcat(cluster_featureNormPooled.(strain),cluster_featureNorm.(strain){fileCtr});
        sw_featureNormPooled.(strain) = vertcat(sw_featureNormPooled.(strain),sw_featureNorm.(strain){fileCtr});
    end
    set(0,'CurrentFigure',mw_featurePooledFig)
    histogram(mw_featureNormPooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    set(0,'CurrentFigure',cluster_featurePooledFig)
    histogram(cluster_featureNormPooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    if ~applySwNormalisation
        set(0,'CurrentFigure',sw_featurePooledFig)
        histogram(sw_featureNormPooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    end
end

%% format and export histograms

% mw plot
set(0,'CurrentFigure',mw_featureFig)
set(gca, 'YScale', yscale)
legend(legends,'Location','eastoutside')
title(['multiworm ' feature])
if ~applySwNormalisation
    xlabel([feature ' ' unit])
else
    xlabel(['relative ' feature])
end
if strcmp(histogramNormalisation,'pdf')
    ylabel('probability')
elseif strcmp(histogramNormalisation,'count')
    ylabel('count')
end
if exist('histogramXLim')
    xlim(histogramXLim)
end
figurename = ['figures/' feature '_' strainSet '_mw_' yscale];
if saveResults
    exportfig(mw_featureFig,[figurename '.eps'],exportOptions)
end

% mw pooled plot
set(0,'CurrentFigure',mw_featurePooledFig)
set(gca, 'YScale', yscale)
legend(legends,'Location','eastoutside')
title(['multiworm ' feature])
if ~applySwNormalisation
    xlabel([feature ' ' unit])
else
    xlabel(['relative ' feature])
end
if strcmp(histogramNormalisation,'pdf')
    ylabel('probability')
elseif strcmp(histogramNormalisation,'count')
    ylabel('count')
end
if exist('histogramXLim')
    xlim(histogramXLim)
end
figurename = ['figures/' feature '_' strainSet '_mw_Pooled_' yscale];
if saveResults
    exportfig(mw_featurePooledFig,[figurename '.eps'],exportOptions)
end

% cluster plot
set(0,'CurrentFigure',cluster_featureFig)
set(gca, 'YScale', yscale)
legend(legends,'Location','eastoutside')
title(['cluster ' feature])
if ~applySwNormalisation
    xlabel([feature ' ' unit])
else
    xlabel(['relative ' feature])
end
if strcmp(histogramNormalisation,'pdf')
    ylabel('probability')
elseif strcmp(histogramNormalisation,'count')
    ylabel('count')
end
if exist('histogramXLim')
    xlim(histogramXLim)
end
figurename = ['figures/' feature '_' strainSet '_cluster_' yscale];
if saveResults
    exportfig(cluster_featureFig,[figurename '.eps'],exportOptions)
end

% cluster pooled plot
set(0,'CurrentFigure',cluster_featurePooledFig)
set(gca, 'YScale', yscale)
legend(legends,'Location','eastoutside')
title(['cluster ' feature])
if ~applySwNormalisation
    xlabel([feature ' ' unit])
else
    xlabel(['relative ' feature])
end
if strcmp(histogramNormalisation,'pdf')
    ylabel('probability')
elseif strcmp(histogramNormalisation,'count')
    ylabel('count')
end
if exist('histogramXLim')
    xlim(histogramXLim)
end
figurename = ['figures/' feature '_' strainSet '_cluster_Pooled_' yscale];
if saveResults
    exportfig(cluster_featurePooledFig,[figurename '.eps'],exportOptions)
end

if ~applySwNormalisation
    % single worm plot
    set(0,'CurrentFigure',sw_featureFig)
    set(gca, 'YScale', yscale)
    legend(legends,'Location','eastoutside')
    title(['single worm ' feature])
    xlabel([feature ' ' unit])
    if strcmp(histogramNormalisation,'pdf')
        ylabel('probability')
    elseif strcmp(histogramNormalisation,'count')
        ylabel('count')
    end
    if exist('histogramXLim')
        xlim(histogramXLim)
    end
    figurename = ['figures/' feature '_' strainSet '_sw_' yscale];
    if saveResults
        exportfig(sw_featureFig,[figurename '.eps'],exportOptions)
    end
    
    % single worm pooled plot
    set(0,'CurrentFigure',sw_featurePooledFig)
    set(gca, 'YScale', yscale)
    legend(legends,'Location','eastoutside')
    title(['single worm ' feature])
    xlabel(feature)
    if strcmp(histogramNormalisation,'pdf')
        ylabel('probability')
    elseif strcmp(histogramNormalisation,'count')
        ylabel('count')
    end
    if exist('histogramXLim')
        xlim(histogramXLim)
    end
    figurename = ['figures/' feature '_' strainSet '_sw_Pooled_' yscale];
    if saveResults
        exportfig(sw_featurePooledFig,[figurename '.eps'],exportOptions)
    end
end