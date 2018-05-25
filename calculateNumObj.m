close all
clear

%% script extracts perimeter and area multi-worm blob features for each of the strains of interest, normalise them against corresponding single worm blob feature within each replicate, and plots a histogram for each strain

%% set parameters
strainSet = 'controls'; % 'controls','divergent','all'
clusterArea = 4;
phaseRestrict = false; % phaseRestrict cuts out the first 15 min of each video
poolReplicates = false;
saveResults = false;
maxNumReplicates = 3;

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',30,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',25,...
    'LineWidth',3);

addpath('auxiliary/')

% load the strain names included in the specified strainSet
load(['strainsList/' strainSet '.mat'])
% get list of file names for each strain
[strainFileList] = getFileList(strains);
% generate colormap for plotting each strain
colorMap = distinguishable_colors(length(strains));

numObjFig = figure; hold on

%% go through each strain
for strainCtr = 1:length(strains)
    filenames = strainFileList.([strains{strainCtr} 'List_40']);
    % if there are many files, then subsample recordings
    if length(filenames)>maxNumReplicates
        fileInd = datasample(1:length(filenames),maxNumReplicates,'Replace',false);
    else
        fileInd = 1:length(filenames);
        maxNumReplicates = length(filenames);
    end
    
    %% go through each recording
    for fileCtr = 1:length(fileInd)
        numObj.(strains{strainCtr})= cell(length(fileInd),1);
        
        %% load data
        filename = filenames{fileCtr}
        if  ~contains(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_16.1')...
                & ~contains(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/agg_11.3')...
                & ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_11.3_180218/11.3_7_eca259_none_a1_Set0_Pos0_Ch3_18022018_170514_skeletons.hdf5')...
                & ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_2.3_180119/2.3_6_nic261_f6_Set0_Pos0_Ch1_19012018_145323_skeletons.hdf5')...
            trajData = h5read(filename,'/trajectories_data');
            blobFeats = h5read(filename,'/blob_features');
            % skelData = h5read(filename,'/skeleton'); % in pixels
            frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
            % features = h5read(strrep(filename,'skeletons','featuresN'),'/features_timeseries/');
            
            %% generate logical index for cluster and single worms
%             [singleWormLogInd,multiWormLogInd,clusterLogInd] = findClusters(trajData,blobFeats,clusterArea);
%             
%             % apply phase restriction (cuts out the first 15 min of each 45 min video)
%             if phaseRestrict
%                 frameLogInd = trajData.frame_number>frameRate*60*15;
%                 singleWormLogInd = singleWormLogInd & frameLogInd;
%                 multiWormLogInd = multiWormLogInd & frameLogInd;
%                 clusterLogInd = clusterLogInd & frameLogInd;
%             end
%             
%             % filter out manually labeled bad entries
%             if isfield(trajData,'worm_label')
%                 singleWormLogInd = singleWormLogInd & trajData.worm_label~=3;
%                 multiWormLogInd = multiWormLogInd & trajData.worm_label~=3;
%                 clusterLogInd = clusterLogInd & trajData.worm_label~=3;
%             end
%             
%             % filter out any blobs that do not persist for more than 1 second
%             tempBlobLogInd = findTempBlobs(trajData,blobFeats,frameRate);
%             singleWormLogInd = singleWormLogInd & ~tempBlobLogInd;
%             multiWormLogInd = multiWormLogInd & ~tempBlobLogInd;
%             clusterLogInd = clusterLogInd & ~tempBlobLogInd;
            
            %% analyse features
            % calculate num of obj
            for frameCtr = 1:max(trajData.frame_number)
                frameLogInd = trajData.frame_number == frameCtr;
                numObj.(strains{strainCtr}){fileCtr}=trajData.worm_index_joined(frameLogInd);
            end
%             % separate blobSpeeds based on cluster status
%             blobSpeed.(strains{strainCtr}){fileCtr} = trajData.blobSpeed(clusterLogInd);
%             swblobSpeed.(strains{strainCtr}){fileCtr} = trajData.blobSpeed(singleWormLogInd);
%             % normalise area and perimeter from this movie with sw features from this movie; store value for threshold box plot later
%             blobSpeedNorm.(strains{strainCtr}){fileCtr} = blobSpeed.(strains{strainCtr}){fileCtr}/nanmedian(swblobSpeed.(strains{strainCtr}){fileCtr});
            
            % plot individual experiments
            if ~poolReplicates
                set(0,'CurrentFigure',numObjFig)
                histogram(numObj.(strains{strainCtr}){fileCtr},'BinWidth',frameRate,'Normalization','count','DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
            end
        end
    end
    
    % combine the replicates and plot
    if poolReplicates
        numObjPooled.(strains{strainCtr}) = numObj.(strains{strainCtr}){1}; % initialise new pooling variable with the first replicate
        for fileCtr = 2:length(numObj.(strains{strainCtr}))
            numObjPooled.(strains{strainCtr}) = vertcat(numObjPooled.(strains{strainCtr}),numObj.(strains{strainCtr}){fileCtr});
        end
        set(0,'CurrentFigure',numObjFig)
        histogram(numObjPooled.(strains{strainCtr}),'BinWidth',frameRate,'Normalization','count','DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:))
    end
end

% format and export histogram
set(0,'CurrentFigure',numObjFig)
set(gca, 'YScale', 'linear')
if poolReplicates
    legend(strains,'Location','eastoutside')
end
title('number of tracked objects')
xlabel('frame number')
ylabel('probability')
if poolReplicates
    figurename = ['figures/numObj_' strainSet '_' num2str(maxNumReplicates) 'reps_pooled'];
else
    figurename = ['figures/numObj_' strainSet '_' num2str(maxNumReplicates) 'reps'];
end
if saveResults
    exportfig(numObjFig,[figurename '.eps'],exportOptions)
    %eps2pdf([figurename '.eps'])
    %system(['epstopdf ' figurename '.eps']);
    %system(['rm ' figurename '.eps']);
end

% save data
%     if saveResults
%         save('results/areaPerimeter.mat','perimeterNorm','areaNorm')
%         savefig(perimeterFig,'figures/perimeter.fig')
%         savefig(areaFig,'figures/area.fig')
%     end