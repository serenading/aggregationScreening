clear
close all

%% script extracts features for each of the strains of interest based on downsampled pixel data.

%% set parameterslegends
% set analysis parameters
strainSet = 'all'; % 'controls','divergent','all'
feature = 'pcf'; % specify feature as string. 'pcf' (pair correlation function), 'hc'(hierarchical clustering), 'gf'(giant fluctuation).
maxNumReplicates = 60; % controls have up to 60 reps, divergents up to 15 reps, all other strains up to 5 reps.
sampleFrameEveryNSec = 10; % 10 works
sampleEveryNPixel = 16; % 16 works
saveResults = false;
plotIndividualReps = false;
showFrame = false; % true if monitoring script running: display current downsized masked binary image as script runs
makeDownSampledVid = false; % true if not monitoring script running: generate downsampled video to check afterwards that no obvious non-worm pixels are kept for analysis

% set feature-specific parameters
if strcmp(feature,'hc')
    featVarName = 'branchHeights';
    linkageMethod = 'single'; % 'single' (preferred), 'average','centroid','complete','median','weighted'
    yscale = 'linear';
    xLim = [0 4];
    yLabel = 'probability';
    xLabel = 'inter-wormpixel distance (mm)';
    figTitle = [linkageMethod ' linkage'];
elseif strcmp(feature,'pcf')
    featVarName = 'pcf';
    distBinWidth = 0.1; % in units of mm
    maxDist = 1.5; % in units of mm
    distBins = 0:distBinWidth:maxDist;
    yscale = 'linear';
    xLim = [0 1.2];
    yLabel = 'positional correlation g(r)';
    xLabel = 'distance r (mm)';
    figTitle = 'pair correlation function';
elseif strcmp(feature,'gf')
end

% set default parameters
useIntensityMask = true;
useOnFoodMask = true;
useMovementMask = true;
phaseRestrict = true; % phaseRestrict cuts out the first 15 min of each video
histogramNormalisation = 'pdf'; % 'pdf' by default. 'count' an option
dims = [2048 2048]; % can be read by the following but slow: fileInfo = h5info(maskedVideoFileName); dims = fileInfo.Datasets(2).Dataspace.Size; %[2048,2048,num]
pixelToMicron = 10; % 10 microns per pixel, read by pixelsize = double(h5readatt(filename,'/trajectories_data','microns_per_pixel?))

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
featurePooledFig = figure; hold on
if plotIndividualReps
    featureFig = figure; hold on
end
if showFrame
    sampleFrameFig = figure;
end
% empty the downsampledVideos directory
if makeDownSampledVid
    delete /Users/sding/Documents/AggScreening/downsampledVideos/*.avi 
end

% create legend variable to hold strain name and experiment n numbers
legends = cell(size(strains));
lineHandles = NaN(numel(strains),1);

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
    % update strain n number for figure legend
    legends{strainCtr} = [strain ', n=' num2str(length(fileInd))];
    %% go through each recording
    for fileCtr = 1:length(fileInd)
        %% load data
        filename = filenames{fileInd(fileCtr)};
        trajData = h5read(filename,'/trajectories_data');
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        maskedVideoFileName = strrep(strrep(filename,'Results','MaskedVideos'),'_skeletons.hdf5','.hdf5');
        foodContourCoords = h5read(filename,'/food_cnt_coord');
        %% initialise video if making one
        if makeDownSampledVid
            videoName = strsplit(maskedVideoFileName,'/');
            videoName = ['/Users/sding/Documents/AggScreening/downsampledVideos/' videoName{end}(1:end-5)];
            video = VideoWriter([videoName '.avi']);
            video.FrameRate = 60/sampleFrameEveryNSec*3; % 1 sec video is 3 min actual video (i.e. 180x speed)
            open(video)
        end
        %% generate onfood binary mask
        onFoodMask = poly2mask(foodContourCoords(1,:),foodContourCoords(2,:),dims(1),dims(2));
        % dilate mask to include pixels immediately outside the mask
        structuralElement = strel('disk',64); % dilate foodpatch by 64 pixels = 640 microns, about half a worm length
        onFoodMaskDilate = imdilate(onFoodMask,structuralElement);
        % get the overall area of food patch in micron squared
        overallArea = nnz(onFoodMaskDilate)*pixelToMicron^2;
        %% sample frames
        if phaseRestrict
            startFrameNum = frameRate*60*15; % cuts out the first 15 minutes
        else
            startFrameNum = 0;
        end
        endFrameNum = max(trajData.frame_number);
        numFrames = floor(numel(startFrameNum:endFrameNum)/frameRate/sampleFrameEveryNSec);
        sampleFrames = datasample(startFrameNum:endFrameNum,numFrames);
        %% feature-specific set up
        if strcmp(feature,'hc')
            branchHeights.(strain){fileCtr}= cell(numFrames,1); % pre-allocate cells to hold branch height values
        elseif strcmp(feature,'pcf')
            pcf.(strain){fileCtr} = NaN(length(distBins) - 1,numFrames);
        end
        %% go through each frame to generate downsampled frames
        maskedImageStack = NaN(numel(1:sampleEveryNPixel:dims(1)),numel(1:sampleEveryNPixel:dims(2)),numFrames);
        originalImageStack = NaN(numel(1:sampleEveryNPixel:dims(1)),numel(1:sampleEveryNPixel:dims(2)),numFrames);
        for frameCtr = 1:numFrames
            % load the frame
            imageFrame = h5read(maskedVideoFileName,'/mask',[1,1,double(sampleFrames(frameCtr))],[dims(1),dims(2),1]);
            imageFrame = rot90(fliplr(imageFrame)); % flip and rotate the image to match the orientation as displayed in Tierpsy Tracker so it works with food contours
            maskedImage = imageFrame;
            % apply various masks to get binary image of worm/nonworm pixels
            if useIntensityMask
                % generate binary segmentation based on black/white contrast
                maskedImage = maskedImage>0 & maskedImage<70;
            end
            if useOnFoodMask
                % generate binary segmentation based on on/off mask
                maskedImage = maskedImage & onFoodMaskDilate;
            end
            % downsample masked image
            downsampleMaskedImage = maskedImage(1:sampleEveryNPixel:dims(1),1:sampleEveryNPixel:dims(2));
            % add downsampled image to image stack
            maskedImageStack(:,:,frameCtr) = downsampleMaskedImage;
            % generate downsampled, unmasked frame for side by side masking comparison
            if showFrame | makeDownSampledVid
                intensityMaskedImageFrame = maskedImage>0 & maskedImage<70;
                originalImageStack(:,:,frameCtr) = intensityMaskedImageFrame(1:sampleEveryNPixel:dims(1),1:sampleEveryNPixel:dims(2));
            end
        end
        % generate no movement mask
        movementMask = std(maskedImageStack,0,3)>0;
        for frameCtr = 1:numFrames
            if useMovementMask
                maskedImageStack(:,:,frameCtr) = maskedImageStack(:,:,frameCtr) & movementMask;
            end
            % display frame
            if showFrame
            set(0,'CurrentFigure',sampleFrameFig)
            imshow([originalImageStack(:,:,frameCtr) maskedImageStack(:,:,frameCtr)])
            end
            % write frame to video
            if makeDownSampledVid
                writeVideo(video,[originalImageStack(:,:,frameCtr) maskedImageStack(:,:,frameCtr)])
            end
            % calculate feature
            N = nnz(maskedImageStack(:,:,frameCtr));% downsampleBinaryImage);
            if N>1 % need at least two worms in frame
                [x,y] = find(maskedImageStack(:,:,frameCtr));% downsampleBinaryImage);
                pairDists = pdist([x y]*sampleEveryNPixel*pixelToMicron/1000); % pairDist in mm;
                if strcmp(feature,'hc')
                    clustTree = linkage(pairDists,linkageMethod);
                    branchHeights.(strain){fileCtr}{frameCtr} = clustTree(:,3);
                    % dendrogram(clustTree,0,'Reorder',optimalleaforder(clustTree,pairDists));
                elseif strcmp(feature,'pcf')
                    pcf.(strain){fileCtr}(:,frameCtr) = histcounts(pairDists,distBins,'Normalization','count'); % radial distribution function
                    pcf.(strain){fileCtr}(:,frameCtr) = pcf.(strain){fileCtr}(:,frameCtr)'.*overallArea ...
                        ./(pi*(distBins(2:end).^2 - (distBins(2:end) - distBinWidth).^2)*N*(N-1)/2); % normalisation by N(N-1)/2 as pdist doesn't double-count pairs
                end
            end
        end
        % pool data from frames
        if strcmp(feature,'hc')
            branchHeights.(strain){fileCtr} = vertcat(branchHeights.(strain){fileCtr}{:});
        end
        % close video
        if makeDownSampledVid
            close(video)
        end
        if plotIndividualReps
            set(0,'CurrentFigure',featureFig)
            if strcmp(feature, 'hc')
                histogram(branchHeights.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:),'BinWidth',0.2)
            elseif strcmp(feature,'pcf')
                boundedline(distBins(2:end)-distBinWidth/2,nanmean(pcf.(strain){fileCtr},2),...
                    [nanstd(pcf.(strain){fileCtr},0,2) nanstd(pcf.(strain){fileCtr},0,2)]./sqrt(nnz(sum(~isnan(pcf.(strain){fileCtr}),2))),...
                    'alpha',featureFig.Children,'cmap',colorMap(strainCtr,:))
            end
        end
    end
    % pool data from multiple files
    if strcmp(feature,'hc')
        branchHeights_pooled.(strain) = vertcat(branchHeights.(strain){:});
    elseif strcmp(feature,'pcf')
        pcf_pooled.(strain) = horzcat(pcf.(strain){:});
    end
    % plot features
    set(0,'CurrentFigure',featurePooledFig)
    if strcmp(feature,'hc')
        histogram(branchHeights_pooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:),'BinWidth',0.2)
    elseif strcmp(feature,'pcf')
        [lineHandles(strainCtr), ~] = boundedline(distBins(2:end)-distBinWidth/2,nanmean(pcf_pooled.(strain),2),...
            [nanstd(pcf_pooled.(strain),0,2) nanstd(pcf_pooled.(strain),0,2)]./sqrt(nnz(sum(~isnan(pcf_pooled.(strain)),2))),...
            'alpha',featurePooledFig.Children,'cmap',colorMap(strainCtr,:));
    end
end

%% save variable
if saveResults
    save(['results/' feature '_' strainSet '_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel.mat'],featVarName)
end

% format and export
set(featurePooledFig,'PaperUnits','centimeters')
set(0,'CurrentFigure',featurePooledFig)
if strcmp(feature,'pcf')
    legend(featurePooledFig.Children,lineHandles,legends,'Location','eastoutside')
else
    legend(legends,'Location','eastoutside')
end
set(gca, 'YScale', yscale)
if exist('xLim') == 1
    xlim(xLim)
end
xlabel(xLabel)
ylabel(yLabel)
title(figTitle,'FontWeight','normal')
if saveResults
    figurename = ['figures/' feature '_' strainSet '_' featVarName '_Pooled_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel'];
    exportfig(featurePooledFig,[figurename '.eps'],exportOptions)
end

if plotIndividualReps
    set(featureFig,'PaperUnits','centimeters')
    set(0,'CurrentFigure',featureFig)
    legend(legends,'Location','eastoutside')
    set(gca, 'YScale', yscale)
    if exist('xLim') == 1
        xlim(xLim)
    end
    xlabel(xLabel)
    ylabel(yLabel)
    title(figTitle)
    if saveResults
        figurename = ['figures/' feature '_' strainSet '_' featVarName '_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel'];
        exportfig(featureFig,[figurename '.eps'],exportOptions)
    end
end