clear
close all

%% script extracts blob features for each of the strains of interest, normalise them against corresponding single worm blob feature within each replicate, and plots a histogram for each strain

%% set parameters
% set analysis parameters
strainSet = 'controls'; % 'controls','divergent','all'
feature = 'speed'; % specify feature as string
clusterArea = 4;
phaseRestrict = true; % phaseRestrict cuts out the first 15 min of each video
saveResults = false;
maxNumReplicates = 100;

% set plotting parameters
histogramNormalisation = 'pdf'; % 'pdf' by default. 'count' an option
yscale = 'linear'; % 'linear' by default. 'log' an option

% set feature-specific parameters
if strcmp(feature,'speed')
    histogramXLim = [0 5];
elseif strcmp(feature,'area') | strcmp(feature,'perimeter')
    yscale = 'log';
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
% load the strain names included in the specified strainSet
load(['strainsList/' strainSet '.mat'])
% get list of file names for each strain
[strainFileList,~,~] = getFileList(strains);
% generate colormap for plotting each strain
colorMap = distinguishable_colors(length(strains));
% create empty figures
mw_featureFig = figure; hold on
mw_featurePooledFig = figure; hold on
cluster_featureFig = figure; hold on
cluster_featurePooledFig = figure; hold on
% create legend variable to hold strain name and experiment n numbers
legends = cell(size(strains));

%% go through each strain
for strainCtr = 1:length(strains)
    filenames = strainFileList.([strains{strainCtr} 'List_40']);
    % if there are many files, then subsample recordings without replacement
    if length(filenames)>maxNumReplicates
        fileInd = datasample(1:length(filenames),maxNumReplicates,'Replace',false);
    else
        fileInd = 1:length(filenames);
    end
    
    %% go through each recording
    for fileCtr = 1:length(fileInd)
        
        %% load data
        filename = filenames{fileCtr}
        if ~contains (filename, 'Agg_16.1')
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
            end
            % read features
            sw_feature.(strains{strainCtr}){fileCtr} = blobFeats.(feature)(singleWormLogInd);
            mw_feature.(strains{strainCtr}){fileCtr} = blobFeats.(feature)(multiWormLogInd);
            cluster_feature.(strains{strainCtr}){fileCtr} = blobFeats.(feature)(clusterLogInd);
            % normalise area and perimeter from this movie with sw features from this movie; store value for threshold box plot later
            mw_featureNorm.(strains{strainCtr}){fileCtr} = mw_feature.(strains{strainCtr}){fileCtr}/nanmedian(sw_feature.(strains{strainCtr}){fileCtr});
            cluster_featureNorm.(strains{strainCtr}){fileCtr} = cluster_feature.(strains{strainCtr}){fileCtr}/nanmedian(sw_feature.(strains{strainCtr}){fileCtr});
            % update strain n number for figure legend
            legends{strainCtr} = [strains{strainCtr} ', n=' num2str(length(mw_featureNorm.(strains{strainCtr})))];
            
            %% plot individual experiments
            set(0,'CurrentFigure',mw_featureFig)
            histogram(mw_featureNorm.(strains{strainCtr}){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
            set(0,'CurrentFigure',cluster_featureFig)
            histogram(cluster_featureNorm.(strains{strainCtr}){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
        end
    end
    
    %% combine the replicates and plot
    mw_featureNormPooled.(strains{strainCtr}) = mw_featureNorm.(strains{strainCtr}){1};
    cluster_featureNormPooled.(strains{strainCtr}) = cluster_featureNorm.(strains{strainCtr}){1};
    for fileCtr = 2:length(mw_featureNorm.(strains{strainCtr}))
        mw_featureNormPooled.(strains{strainCtr}) = vertcat(mw_featureNormPooled.(strains{strainCtr}),mw_featureNorm.(strains{strainCtr}){fileCtr});
        cluster_featureNormPooled.(strains{strainCtr}) = vertcat(cluster_featureNormPooled.(strains{strainCtr}),cluster_featureNorm.(strains{strainCtr}){fileCtr});
    end
    set(0,'CurrentFigure',mw_featurePooledFig)
    histogram(mw_featureNormPooled.(strains{strainCtr}),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    set(0,'CurrentFigure',cluster_featurePooledFig)
    histogram(cluster_featureNormPooled.(strains{strainCtr}),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
end

%% format and export histograms

% mw plot
set(0,'CurrentFigure',mw_featureFig)
set(gca, 'YScale', yscale)
legend(legends,'Location','eastoutside')
title(['multiworm ' feature])
xlabel(['relative ' feature])
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
xlabel(['relative ' feature])
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
xlabel(['relative ' feature])
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
xlabel(['relative ' feature])
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