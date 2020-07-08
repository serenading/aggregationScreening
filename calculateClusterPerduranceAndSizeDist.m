close all
clear

%% Script extracts cluster perdurance and cluster size distribution statistics and optionally generates diagnostic plots for the divergent panel
% author: @serenading. July 2020

addpath('auxiliary/')

% TODO: what stats to extract for cluster size distribution?

%% Specify analysis parameters

% which HD is plugged in right now?
whichHD = 1; % Scalar. Specify 1 or 2. Specify 0 if both 1 and 2 are plugged in.
% set which feature extraction timestamp to use (doesn't matter as we are calculating our own feats)
extractStamp = '20200519_153722';
% which worm density?
wormNum = 5;
% some file attributes
frameRate = 25;
pixelToMicron = 10; % 1 pixel is 10 microns
% cluster size distriution parameters
timeWindowStarts = [0:10:40]; % starting minute of time windows
windowDuration = 5; % time window duration in minutes
% plot for divergent strains
plot4Divergent = true;

%% load feature table
if wormNum == 40 
    %% find which files are on which HD
    featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
    % load divergent strain names (these are copied onto HD1, the rest onto HD2)
    load('strainsList/divergent.mat','strains');
    % initialise onHD1 variable to say if the file is on HD1
    onHD1 = false(size(featureTable,1),1);
    % go through the divergent strains to set their onHD1 status to true
    for strainCtr = 1:numel(strains)
        strain = strains{strainCtr};
        strainLogInd = strcmp(featureTable.strain_name,strain);
        onHD1(strainLogInd) = true;
    end
    % append onHD1 variable to featureTable
    featureTable.onHD1 = onHD1;
    %% Go through each file that is present on this HD to extract features
    if whichHD ==1
        fileInd = find(featureTable.onHD1);
    elseif whichHD ==2
        fileInd = find(~featureTable.onHD1);
    elseif whichHD ==0 % both HDs plugged in
        fileInd = 1:numel(featureTable.onHD1);
    end
elseif wormNum ==5 
    %% all files are on HD1
    assert(whichHD == 1, ['HD' num2str(whichHD) ' is plugged in, but all five worm data are on HD1.']);
    featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/fiveWormFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
    fileInd = 1:size(featureTable,1);
end
%% Initialise and preallocate
% filenames
filenames = repmat({''},numel(fileInd),1);
% strainnames
strainnames = repmat({''},numel(fileInd),1);
% cell array to hold perdurance vectors
perduranceC = cell(numel(fileInd),1);
% cell array to hold cluster size vectors
clusterSizeC = cell(numel(fileInd),numel(timeWindowStarts));
% legend for cluster size distribuion analysis
legendtext_cs = repmat({''},numel(timeWindowStarts),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conduct cluster analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for fileCtr = 1:numel(fileInd)
    disp(['Extracting features for File ' num2str(fileCtr) ' out of ' num2str(numel(fileInd)) ])
    
    %% Load data
    % get featuresN filename from Results folder
    filename = [strrep(strrep(char(featureTable.dirname(fileInd(fileCtr))),'/Volumes/behavgenom_archive$/Serena/','/Volumes/Ashur DT2/'),'MaskedVideos','Results'),...
        '/', strrep(char(featureTable.filename(fileInd(fileCtr))),'.hdf5','_featuresN.hdf5')];
    % replace HD path only if both HD are plugged in and 1 is plugged in before 2.
    if whichHD == 0 && ~featureTable.onHD1(fileInd(fileCtr))
        filename = strrep(filename,'DT2/','DT2 1/');
    end
    % record filename
    filenames{fileCtr} = char(featureTable.filename(fileInd(fileCtr)));
    % record strain names
    strainnames{fileCtr} = char(featureTable.strain_name(fileInd(fileCtr)));
    % load tracking data from featuresN.hdf5
    trajData = h5read(filename,'/trajectories_data');
    blobFeats = h5read(filename,'/blob_features');
    
    %% initialise for cluster perdurance analysis
    wormInd = unique(trajData.worm_index_joined);
    n_worms = numel(wormInd);
    perdurance = NaN(n_worms,1);
    
    %% cluster perdurance analysis
    % go through each worm
    for wormCtr = 1:n_worms
        wormIdx = wormInd(wormCtr);
        % get the frames in which this worm is present
        rowLogInd = trajData.worm_index_joined == wormIdx;
        wormFrames = trajData.frame_number(rowLogInd);
        % check that the frame numbers for this worm are continuous
        assert(nnz(diff(wormFrames)~=1)==0,['Frame numbers are not consecutive for worm index ' num2str(wormIdx) '.'])
        % record perdurance vector of this worm index
        perdurance(wormCtr) = numel(wormFrames);
    end
    % convert unit from frames to seconds
    perdurance = perdurance/frameRate;
    % add to perdurance cell array
    perduranceC{fileCtr} = perdurance;
    
    %% cluster size distribution analysis
    % get valid frames from each time window
    for twCtr = 1:numel(timeWindowStarts)
        validFrameRange = frameRate*60*[timeWindowStarts(twCtr),timeWindowStarts(twCtr)+windowDuration];
        validFrames = validFrameRange(1):validFrameRange(end);
        % get cluster sizes from valid frames
        rowLogInd = ismember(trajData.frame_number,validFrames);
        clusterSizes = trajData.area(rowLogInd);
        % convert cluster size to mm^2
        clusterSizes = clusterSizes*pixelToMicron^2/1e6;
        % add to cluster size cell array
        clusterSizeC{fileCtr,twCtr} = clusterSizes;
        if fileCtr ==1 % only need to get the legend once
            legendtext_cs{twCtr} = [num2str(timeWindowStarts(twCtr)) '-' num2str(timeWindowStarts(twCtr)+windowDuration) ' min'];
        end
    end
end

%% save intermediate results
timestamp=datetime('today','Format','yyyyMMdd');
if wormNum == 40
    save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/newFeatsToAdd/clusterAnalysis_' char(timestamp) '.mat'],'perduranceC','clusterSizeC','filenames','strainnames','-v7.3');
elseif wormNum == 5
    save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/newFeatsToAdd/clusterAnalysis_' char(timestamp) '.mat'],'perduranceC','clusterSizeC','filenames','strainnames','-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate further stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cluster perdurance cumulative survival analysis
% initialise
cluster_perdurance_90prc = NaN(size(filenames));
cluster_perdurance_50prc = NaN(size(filenames));
cluster_perdurance_10prc = NaN(size(filenames));
% calculate time to 90%, 50%, 10% cluster survival (in seconds)
for fileCtr = 1:numel(filenames)
    [ecdfy,ecdfx] = ecdf(perduranceC{fileCtr});
    cluster_perdurance_90prc(fileCtr) = ecdfx(interp1(ecdfy,1:length(ecdfy),1-0.9,'nearest'));
    cluster_perdurance_50prc(fileCtr) = ecdfx(interp1(ecdfy,1:length(ecdfy),1-0.5,'nearest'));
    cluster_perdurance_10prc(fileCtr) = ecdfx(interp1(ecdfy,1:length(ecdfy),1-0.1,'nearest'));
end

%% save new feature table and append new features to the latest featureTable
newFeatureTable = table(filenames,cluster_perdurance_90prc,cluster_perdurance_50prc,cluster_perdurance_10prc);
newFeatureTable.Properties.VariableNames(1) = "filename"; % rename variable from filenames to filename to be consistent with other analyses
if wormNum == 40
    save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/newFeatsToAdd/fortyWormFeaturesTable_' extractStamp '_new_' char(timestamp) '.mat'],'newFeatureTable')
elseif wormNum == 5
    save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/newFeatsToAdd/fiveWormFeaturesTable_' extractStamp '_new_' char(timestamp) '.mat'],'newFeatureTable')
end
appendFeatsToFeatureTable(newFeatureTable,wormNum,extractStamp);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot analysis results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot4Divergent
    
    % initialise
    load('strainsList/divergent.mat','strains');
    clusterPerduranceFig = figure; hold on
    clusterSizeDistributionFig = figure; hold on
    plotcolors_cp = linspecer(numel(strains));
    plotcolors_cs = linspecer(numel(timeWindowStarts));
    
    % go through each strain to pool results
    for strainCtr = 1:numel(strains)
        strain = strains{strainCtr};
        strainInd = find(strcmp(strainnames,strain));
        % pre-allocate to hold pooled values at the strain level
        strainPerdurance = [];
        for twCtr = 1:numel(timeWindowStarts)
            strainClusterSize{twCtr} =[];
        end
        % go through each replicate of the strain
        for repCtr = 1:numel(strainInd)
            strainIdx = strainInd(repCtr);
            % pool cluster perdurance
            strainPerdurance = vertcat(strainPerdurance,perduranceC{strainIdx});
            % pool cluster size distribution for each time window
            for twCtr = 1:numel(timeWindowStarts)
                strainClusterSize{twCtr} = vertcat(strainClusterSize{twCtr},clusterSizeC{strainIdx,twCtr});
            end
        end
        % cluster perdurance cumulative survival analysis
        [ecdfy,ecdfx] = ecdf(strainPerdurance);
        
        % plot cluster perdurance cumulative survival
        set(0,'CurrentFigure',clusterPerduranceFig)
        plot(ecdfx,1-ecdfy,'Color',plotcolors_cp(strainCtr,:)) % gives a smoother curve than the survival function
        
        % plot cluster size distribution
        set(0,'CurrentFigure',clusterSizeDistributionFig)
        subplot(3,5,strainCtr); hold on
        for twCtr = 1:numel(timeWindowStarts)
            histogram(strainClusterSize{twCtr},'Normalization','pdf','DisplayStyle','stairs','EdgeColor',plotcolors_cs(twCtr,:))
        end
        title(strain)
        % format plot
        xlabel('cluster size (mm^2)')
        ylabel('pdf')
        if strainCtr ==1
            legend(legendtext_cs)
        end
    end
    
    % format figures
    % cluster perdurance plot
    set(0,'CurrentFigure',clusterPerduranceFig)
    xlabel('time elapsed (s)')
    ylabel('remaining proportion')
    set(gca, 'YScale', 'log')
    legend(strains)
    
    % cluster size distribution plot
    allaxes = findall(clusterSizeDistributionFig,'type','axes');
    linkaxes(allaxes)
    xlim(allaxes,[0.03,0.3]) % zoom in on areas corresponding to 1-6 worms
    ylim(allaxes,[0,5])
    
end