clear
close all

%% script works with bright field wild isolate aggregation screening dataset to
%% extract or calculate single worm, multi worm, and cluster features for each of the strains of interest,
%% normalise mw and cluster features against corresponding sw feature within each replicate (if specified),
%% and plots a histogram (and any additional feature-specific plot) for each strain.

%% set parameters
% set analysis parameters
strainSet = 'controls'; % 'controls','divergent','all'
feature = 'area'; % specify feature as string. 'area','compactness','perimeter','quirkiness','solidity','speed','perdurance'
saveResults = false;
maxNumReplicates = 2; % controls have up to 60 reps, divergents up to 15 reps, all other strains up to 5 reps.
plotIndividualReps = false;

% set default parameters
clusterArea = 4; % 4 by default
phaseRestrict = true; % phaseRestrict cuts out the first 15 min of each video
histogramNormalisation = 'pdf'; % 'pdf' by default. 'count' an option

% set feature-specific parameters
if strcmp(feature,'area')
    unit = '\mum^2'; % micron squared
    applySwNormalisation = true;
    yscale = 'log';
elseif strcmp(feature,'perdurance')
    unit = 's';
    applySwNormalisation = false; % this feature should not be normalised against single worm
    yscale = 'log';
    histogramXLim = [0 600]; %[0 15000];
elseif strcmp(feature,'speed')
    unit = '\mum/s';
    applySwNormalisation = true;
    histogramXLimNorm = [0 5];
    histogramXLim = [0 100];
    yscale = 'log';
elseif strcmp(feature,'compactness')
    unit = '';
    applySwNormalisation = true;
    yscale = 'log';
elseif strcmp(feature,'quirkiness')
    unit = '';
    yscale = 'log';
    applySwNormalisation = true;
elseif strcmp(feature,'solidity')
    unit = '';
    applySwNormalisation = true;
    yscale = 'log';
elseif strcmp(feature,'perimeter')
    unit = '\mum';
    applySwNormalisation = true;
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
addpath('auxiliary/')
% load the strain names included in the specified strainSet
load(['strainsList/' strainSet '.mat'])
% get list of file names for each strain
[strainFileList,~,~] = getFileList(strains);
% generate colormap for plotting each strain
colorMap = distinguishable_colors(length(strains));
% create empty figures
mw_featurePooledFig = figure; hold on
cluster_featurePooledFig = figure; hold on
sw_featurePooledFig = figure; hold on
if plotIndividualReps
    mw_featureFig = figure; hold on
    cluster_featureFig = figure; hold on
    sw_featureFig = figure; hold on
end
if applySwNormalisation
    mwNorm_featurePooledFig = figure; hold on
    clusterNorm_featurePooledFig = figure; hold on
    if plotIndividualReps
        mwNorm_featureFig = figure; hold on
        clusterNorm_featureFig = figure; hold on
    end
end
if strcmp(feature,'perdurance')
    swPerduranceSurvivalCurveFig = figure; hold on
    mwPerduranceSurvivalCurveFig = figure; hold on
    clusterPerduranceSurvivalCurveFig = figure; hold on
end
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
    if strcmp(feature,'perdurance')
        sw_frameDist = [];
        mw_frameDist = [];
        cluster_frameDist = [];
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
        % feature calculation stage 1 (some features are calculated here, some below)
        if strcmp(feature,'speed')
            blobFeats.speed = calculateBlobSpeed(trajData, blobFeats,frameRate);
        elseif strcmp(feature,'perdurance')
            blobFeats.perdurance = trajData.worm_index_joined; % load variable, calculate feature later
        end
        % read features
        sw_feature.(strain){fileCtr} = blobFeats.(feature)(singleWormLogInd);
        mw_feature.(strain){fileCtr} = blobFeats.(feature)(multiWormLogInd);
        cluster_feature.(strain){fileCtr} = blobFeats.(feature)(clusterLogInd);
        % feature calculation stage 2 (some features are calculated here, some above)
        if strcmp(feature,'perdurance') % calculate perdurance in units of frames
            [sw_feature.(strain){fileCtr},mw_feature.(strain){fileCtr},cluster_feature.(strain){fileCtr},...
                sw_frameDist,mw_frameDist,cluster_frameDist] = calculatePerdurance...
                (blobFeats,trajData,singleWormLogInd,multiWormLogInd,clusterLogInd,sw_frameDist,mw_frameDist,cluster_frameDist);
            sw_feature.(strain){fileCtr} = sw_feature.(strain){fileCtr}/frameRate; % convert unit from frames to seconds
            mw_feature.(strain){fileCtr} = mw_feature.(strain){fileCtr}/frameRate;
            cluster_feature.(strain){fileCtr} = cluster_feature.(strain){fileCtr}/frameRate;
            sw_perdDist = sw_frameDist/frameRate;
            mw_perdDist = mw_frameDist/frameRate;
            cluster_perdDist = cluster_frameDist/frameRate;
        end
        % normalise features from this movie with sw features from this movie (if appropriate); store value for threshold box plot later
        if applySwNormalisation % normalise using single worm data
            medianSwFeat = nanmedian(sw_feature.(strain){fileCtr});
            mwNorm_feature.(strain){fileCtr} = mw_feature.(strain){fileCtr}/medianSwFeat;
            clusterNorm_feature.(strain){fileCtr} = cluster_feature.(strain){fileCtr}/medianSwFeat;
        end
        % update strain n number for figure legend
        legends{strainCtr} = [strain ', n=' num2str(length(mw_feature.(strain)))];
        
        %% plot individual experiments
        if plotIndividualReps
            set(0,'CurrentFigure',mw_featureFig)
            histogram(mw_feature.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
            set(0,'CurrentFigure',cluster_featureFig)
            histogram(cluster_feature.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
            set(0,'CurrentFigure',sw_featureFig)
            histogram(sw_feature.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
            if applySwNormalisation
                set(0,'CurrentFigure',mwNorm_featureFig)
                histogram(mwNorm_feature.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
                set(0,'CurrentFigure',clusterNorm_featureFig)
                histogram(clusterNorm_feature.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
            end
        end
    end
    
    %% combine the replicates and plot distribution
    mw_featurePooled.(strain) = mw_feature.(strain){1};
    cluster_featurePooled.(strain) = cluster_feature.(strain){1};
    sw_featurePooled.(strain) = sw_feature.(strain){1};
    if applySwNormalisation
        mwNorm_featurePooled.(strain) = mwNorm_feature.(strain){1};
        clusterNorm_featurePooled.(strain) = clusterNorm_feature.(strain){1};
    end
    for fileCtr = 2:length(mw_feature.(strain))
        mw_featurePooled.(strain) = vertcat(mw_featurePooled.(strain),mw_feature.(strain){fileCtr});
        cluster_featurePooled.(strain) = vertcat(cluster_featurePooled.(strain),cluster_feature.(strain){fileCtr});
        sw_featurePooled.(strain) = vertcat(sw_featurePooled.(strain),sw_feature.(strain){fileCtr});
        if applySwNormalisation
            mwNorm_featurePooled.(strain) = vertcat(mwNorm_featurePooled.(strain),mwNorm_feature.(strain){fileCtr});
            clusterNorm_featurePooled.(strain) = vertcat(clusterNorm_featurePooled.(strain),clusterNorm_feature.(strain){fileCtr});
        end
    end
    set(0,'CurrentFigure',mw_featurePooledFig)
    histogram(mw_featurePooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    set(0,'CurrentFigure',cluster_featurePooledFig)
    histogram(cluster_featurePooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    set(0,'CurrentFigure',sw_featurePooledFig)
    histogram(sw_featurePooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    if applySwNormalisation
        set(0,'CurrentFigure',mwNorm_featurePooledFig)
        histogram(mwNorm_featurePooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
        set(0,'CurrentFigure',clusterNorm_featurePooledFig)
        histogram(clusterNorm_featurePooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    end
    if strcmp(feature,'perdurance')
        swPerdDist.(strain) = sw_perdDist;
        mwPerdDist.(strain) = mw_perdDist;
        clusterPerdDist.(strain) = cluster_perdDist;
    end
end

%% save variables
if saveResults
    save(['results/' feature '_' strainSet '.mat'],'sw_feature','mw_feature','cluster_feature')
end

%% format and export histograms (and other feature-specific plots)

%% pooled plots: pooled across different replicates
% sw pooled plot
set(0,'CurrentFigure',sw_featurePooledFig)
set(gca, 'YScale', yscale)
legend(legends,'Location','eastoutside')
title(['single worm ' feature])
xlabel([feature ' (' unit ')'])
if strcmp(histogramNormalisation,'pdf')
    ylabel('probability')
elseif strcmp(histogramNormalisation,'count')
    ylabel('count')
end
if exist('histogramXLim')
    xlim(histogramXLim)
end
if strcmp(feature,'area')
    ax = gca;
    ax.XAxis.Exponent = 2;
end
figurename = ['figures/' feature '_' strainSet '_sw_Pooled_' yscale];
if saveResults
    exportfig(sw_featurePooledFig,[figurename '.eps'],exportOptions)
end

% mw pooled plot
set(0,'CurrentFigure',mw_featurePooledFig)
set(gca, 'YScale', yscale)
legend(legends,'Location','eastoutside')
title(['multiworm ' feature])
xlabel([feature ' (' unit ')'])
if strcmp(histogramNormalisation,'pdf')
    ylabel('probability')
elseif strcmp(histogramNormalisation,'count')
    ylabel('count')
end
if exist('histogramXLim')
    xlim(histogramXLim)
end
if strcmp(feature,'area')
    ax = gca;
    ax.XAxis.Exponent = 3;
end
figurename = ['figures/' feature '_' strainSet '_mw_Pooled_' yscale];
if saveResults
    exportfig(mw_featurePooledFig,[figurename '.eps'],exportOptions)
end

% cluster pooled plot
set(0,'CurrentFigure',cluster_featurePooledFig)
set(gca, 'YScale', yscale)
legend(legends,'Location','eastoutside')
title(['cluster ' feature])
xlabel([feature ' (' unit ')'])
if strcmp(histogramNormalisation,'pdf')
    ylabel('probability')
elseif strcmp(histogramNormalisation,'count')
    ylabel('count')
end
if exist('histogramXLim')
    xlim(histogramXLim)
end
if strcmp(feature,'area')
    ax = gca;
    ax.XAxis.Exponent = 3;
end
figurename = ['figures/' feature '_' strainSet '_cluster_Pooled_' yscale];
if saveResults
    exportfig(cluster_featurePooledFig,[figurename '.eps'],exportOptions)
end

if applySwNormalisation
    % mwNorm pooled plot
    set(0,'CurrentFigure',mwNorm_featurePooledFig)
    set(gca, 'YScale', yscale)
    legend(legends,'Location','eastoutside')
    title(['multiworm ' feature])
    xlabel(['relative ' feature])
    if strcmp(histogramNormalisation,'pdf')
        ylabel('probability')
    elseif strcmp(histogramNormalisation,'count')
        ylabel('count')
    end
    if exist('histogramXLimNorm')
        xlim(histogramXLimNorm)
    end
    figurename = ['figures/' feature '_' strainSet '_mwNorm_Pooled_' yscale];
    if saveResults
        exportfig(mwNorm_featurePooledFig,[figurename '.eps'],exportOptions)
    end
    
    % clusterNorm pooled plot
    set(0,'CurrentFigure',clusterNorm_featurePooledFig)
    set(gca, 'YScale', yscale)
    legend(legends,'Location','eastoutside')
    title(['cluster ' feature])
    xlabel(['relative ' feature])
    if strcmp(histogramNormalisation,'pdf')
        ylabel('probability')
    elseif strcmp(histogramNormalisation,'count')
        ylabel('count')
    end
    if exist('histogramXLimNorm')
        xlim(histogramXLimNorm)
    end
    figurename = ['figures/' feature '_' strainSet '_clusterNorm_Pooled_' yscale];
    if saveResults
        exportfig(clusterNorm_featurePooledFig,[figurename '.eps'],exportOptions)
    end
end

%% individual replicate plots
if plotIndividualReps
    % sw plot
    set(0,'CurrentFigure',sw_featureFig)
    set(gca, 'YScale', yscale)
    legend(legends,'Location','eastoutside')
    title(['single worm ' feature])
    xlabel([feature ' (' unit ')'])
    if strcmp(histogramNormalisation,'pdf')
        ylabel('probability')
    elseif strcmp(histogramNormalisation,'count')
        ylabel('count')
    end
    if exist('histogramXLim')
        xlim(histogramXLim)
    end
    if strcmp(feature,'area')
        ax = gca;
        ax.XAxis.Exponent = 2;
    end
    figurename = ['figures/' feature '_' strainSet '_sw_' yscale];
    if saveResults
        exportfig(sw_featureFig,[figurename '.eps'],exportOptions)
    end
    
    % mw plot
    set(0,'CurrentFigure',mw_featureFig)
    set(gca, 'YScale', yscale)
    legend(legends,'Location','eastoutside')
    title(['multiworm ' feature])
    xlabel([feature ' (' unit ')'])
    if strcmp(histogramNormalisation,'pdf')
        ylabel('probability')
    elseif strcmp(histogramNormalisation,'count')
        ylabel('count')
    end
    if exist('histogramXLim')
        xlim(histogramXLim)
    end
    if strcmp(feature,'area')
        ax = gca;
        ax.XAxis.Exponent = 3;
    end
    figurename = ['figures/' feature '_' strainSet '_mw_' yscale];
    if saveResults
        exportfig(mw_featureFig,[figurename '.eps'],exportOptions)
    end
    
    % cluster plot
    set(0,'CurrentFigure',cluster_featureFig)
    set(gca, 'YScale', yscale)
    legend(legends,'Location','eastoutside')
    title(['cluster ' feature])
    xlabel([feature ' (' unit ')'])
    if strcmp(histogramNormalisation,'pdf')
        ylabel('probability')
    elseif strcmp(histogramNormalisation,'count')
        ylabel('count')
    end
    if exist('histogramXLim')
        xlim(histogramXLim)
    end
    if strcmp(feature,'area')
        ax = gca;
        ax.XAxis.Exponent = 3;
    end
    figurename = ['figures/' feature '_' strainSet '_cluster_' yscale];
    if saveResults
        exportfig(cluster_featureFig,[figurename '.eps'],exportOptions)
    end
    
    if applySwNormalisation
        % mwNorm plot
        if applySwNormalisation
            set(0,'CurrentFigure',mwNorm_featureFig)
            set(gca, 'YScale', yscale)
            legend(legends,'Location','eastoutside')
            title(['multiwormNorm ' feature])
            xlabel(['relative ' feature])
            if strcmp(histogramNormalisation,'pdf')
                ylabel('probability')
            elseif strcmp(histogramNormalisation,'count')
                ylabel('count')
            end
            if exist('histogramXLimNorm')
                xlim(histogramXLimNorm)
            end
            figurename = ['figures/' feature '_' strainSet '_mwNorm_' yscale];
            if saveResults
                exportfig(mwNorm_featureFig,[figurename '.eps'],exportOptions)
            end
        end
        
        % clusterNorm plot
        set(0,'CurrentFigure',clusterNorm_featureFig)
        set(gca, 'YScale', yscale)
        legend(legends,'Location','eastoutside')
        title(['cluster ' feature])
        xlabel(['relative ' feature])
        if strcmp(histogramNormalisation,'pdf')
            ylabel('probability')
        elseif strcmp(histogramNormalisation,'count')
            ylabel('count')
        end
        if exist('histogramXLimNorm')
            xlim(histogramXLimNorm)
        end
        figurename = ['figures/' feature '_' strainSet '_clusterNorm_' yscale];
        if saveResults
            exportfig(clusterNorm_featureFig,[figurename '.eps'],exportOptions)
        end
    end
end

%% feature-specific plots
% perdurance cumulative survival figures
if strcmp(feature,'perdurance')
    % go through each strain to plot
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        
        set(0,'CurrentFigure',swPerduranceSurvivalCurveFig)
        [ecdfy,ecdfx] = ecdf(swPerdDist.(strain)); % frameRate is 25fps for this dataset
        plot(ecdfx,1-ecdfy,'Color',colorMap(strainCtr,:)) % gives a smoother curve than the survival function
        %ecdf(swFrameDist.(strain),'function','survivor','alpha',0.01,'bounds','on')
        hold on
        
        set(0,'CurrentFigure',mwPerduranceSurvivalCurveFig)
        [ecdfy,ecdfx] = ecdf(mwPerdDist.(strain));
        plot(ecdfx,1-ecdfy,'Color',colorMap(strainCtr,:)) % gives a smoother curve than the survival function
        %ecdf(mwFrameDist.(strain),'function','survivor','alpha',0.01,'bounds','on')
        hold on
        
        set(0,'CurrentFigure',clusterPerduranceSurvivalCurveFig)
        [ecdfy,ecdfx] = ecdf(clusterPerdDist.(strain));
        plot(ecdfx,1-ecdfy,'Color',colorMap(strainCtr,:)) % gives a smoother curve than the survival function
        %ecdf(clusterFrameDist.(strain),'function','survivor','alpha',0.01,'bounds','on')
        hold on
    end
    
    % format survival plots
    set(0,'CurrentFigure',swPerduranceSurvivalCurveFig)
    set(gca, 'YScale', yscale)
    title(['single worm ' feature ' survival'])
    xlabel('time elapsed (s)')
    ylabel('remaining proportion')
    xlim(histogramXLim)
    legend(legends,'Location','eastoutside')
    figurename = ['figures/' feature '_' strainSet '_swSurvival_'];
    if saveResults
        exportfig(swPerduranceSurvivalCurveFig,[figurename '.eps'],exportOptions)
    end
    
    set(0,'CurrentFigure',mwPerduranceSurvivalCurveFig)
    set(gca, 'YScale', yscale)
    title(['multi worm ' feature ' survival'])
    xlabel('time elapsed (s)')
    ylabel('remaining proportion')
    xlim(histogramXLim/2)
    legend(legends,'Location','eastoutside')
    figurename = ['figures/' feature '_' strainSet '_mwSurvival_'];
    if saveResults
        exportfig(mwPerduranceSurvivalCurveFig,[figurename '.eps'],exportOptions)
    end
    
    set(0,'CurrentFigure',clusterPerduranceSurvivalCurveFig)
    set(gca, 'YScale', yscale)
    title(['cluster ' feature ' survival'])
    xlabel('time elapsed (s)')
    ylabel('remaining proportion')
    xlim(histogramXLim/2)
    legend(legends,'Location','eastoutside')
    figurename = ['figures/' feature '_' strainSet '_clusterSurvival_'];
    if saveResults
        exportfig(clusterPerduranceSurvivalCurveFig,[figurename '.eps'],exportOptions)
    end
    
    % save variable
    if saveResults
        save(['results/perduranceSurvival_' strainSet '.mat'],'swFrameDist','mwFrameDist','clusterFrameDist')
    end
end