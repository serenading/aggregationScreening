clear
close all

%% script works with bright field wild isolate aggregation screening dataset (40 worms) to
% extract features for each of the strains of interest based on downsampled pixel data (options to show downsampled frame and make video);
% plot a distribution for each strain; 
% and (optionally) fit a function to the distribution (i.e. exponential decay function of pcf)

%% set parameterslegends
% set analysis parameters
strainSet = 'divergent'; % 'controls','divergent','all'
feature = 'ac'; % specify feature as string. 'pcf' (pair correlation function), 'hc'(hierarchical clustering), 'ac' (auto-correlation),'gf'(giant fluctuation).
maxNumReplicates = 60; % controls have up to 60 reps, divergents up to 15 reps, all other strains up to 5 reps.
saveResults = false;
plotIndividualReps = true;
showFrame = true; % true if monitoring script running: display current downsized masked binary image as script runs
makeDownSampledVid = false; % true if not monitoring script running: generate downsampled video to check afterwards that no obvious non-worm pixels are kept for analysis
fitModel = false;
    
% set feature-specific parameters 
% investigate using switch case instead of if strcmp
if strcmp(feature,'hc')
    sampleFrameEveryNSec = 10;
    sampleEveryNPixel = 8; % 4,8,16 for hc
    featVarName = 'branchHeights';
    linkageMethod = 'single'; % 'single' (preferred), 'average','centroid','complete','median','weighted'
    histogramNormalisation = 'pdf'; % 'pdf' (preferred) or 'count'
    xLim = [0 4];
    figTitle = [linkageMethod ' linkage'];
    xLabel = 'inter-wormpixel distance (mm)';
    if strcmp(histogramNormalisation,'pdf')
        yLabel = 'probability density function';
    else
        yLabel = histogramNormalisation;
    end
    if sampleEveryNPixel == 16
        yscale = 'linear';
        histBinWidth = 0.2;
    elseif sampleEveryNPixel == 8
        yscale = 'log';
        histBinWidth = 0.1;
    elseif sampleEveryNPixel == 4
        yscale = 'log';
        histBinWidth = 0.1;
    end
elseif strcmp(feature,'pcf')
    sampleFrameEveryNSec = 10;
    sampleEveryNPixel = 1;
    featVarName = 'pcf';
    distBinWidth = 0.1; % in units of mm
    maxDist = 1.5; % in units of mm
    distBins = 0:distBinWidth:maxDist;
    yscale = 'linear';
    xLim = [0 1.2];
    yLabel = 'positional correlation g(r)';
    xLabel = 'distance r (mm)';
    figTitle = 'pair correlation function';
    if fitModel
        modelFun = {'exp1'}; % 'exp1' or 'exp2'
    end
elseif strcmp(feature,'ac')
    sampleFrameEveryNSec = 0.32; % use multiples of 0.04 (= 1 frame/s at 25 fps). 0.32 = 3 frames/s
    sampleEveryNPixel = 8;
    featVarName = 'ac';
    maxLag = 300; % maximum lag time in seconds; currently do not set maxLag > 900s, as script extracts frames from twice the duration so the final frame has a full lag time.
    yscale = 'linear';
    yLabel = 'correlation coefficient';
    xLabel = 'lag (frames)';
    figTitle = 'video auto correlation';
    xLim = [0 maxLag/sampleFrameEveryNSec];
    if fitModel
        modelFun = {'mexp3'}; % 'mexp3','exp3' or 'power2'
        if contains(modelFun, 'mexp3') % exp3 with offset
            mexp3 = fittype('a*exp(b*x)+c*exp(d*x)+(1-a-c-f)*exp(e*x)+f');
            fitOptions = fitoptions(mexp3);
            fitOptions.Lower = [0 -0.2 0 -0.2 -0.2 0.05];
            fitOptions.Upper = [1 0 1 0 0 0.3];
            fitOptions.StartPoint = [0.5 -0.2 0.2 -0.02 -0.002 0.1];
        elseif contains(modelFun,'exp3') % triple exponential
            exp3 = fittype('a*exp(b*x)+c*exp(d*x)+(1-a-c)*exp(e*x)');
            fitOptions = fitoptions(exp3);
            fitOptions.Lower = [0 -0.2 0 -0.2 -0.2];
            fitOptions.Upper = [1 0 1 0 0];
            fitOptions.StartPoint = [0.33 -0.2 0.33 -0.02 -0.002];
        elseif contains(modelFun, 'snexp2')
            snexp2 = fittype('a*exp(-b*x^f)+c*exp(-d*x)+(1-a-c)'); % Stretched exponential and Normal exponential plus offset
            fitOptions = fitoptions(snexp2);
            fitOptions.Lower = [0 0 0 0 0];
            fitOptions.Upper = [1 0.2 1 0.2 0.7];
            fitOptions.StartPoint = [0.35 0.1 0.5 0.1 0.5];
        end
    end
elseif strcmp(feature,'gf')
end

% set default parameters
useIntensityMask = true;
useOnFoodMask = true;
useMovementMask = true;
phaseRestrict = true; % phaseRestrict cuts out the first 15 min of each video
pixelToMicron = 10; % 10 microns per pixel, read by pixelsize = double(h5readatt(filename,'/trajectories_data','microns_per_pixel?))
dims = [2048 2048]; % can be read by the following but slow: fileInfo = h5info(maskedVideoFileName); dims = fileInfo.Datasets(2).Dataspace.Size; %[2048,2048,num]

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
[strainFileList,~] = getFileList(strains);
% generate colormap for plotting each strain
colorMap = distinguishable_colors(length(strains));
% create empty figures
featurePooledFig = figure; hold on
if plotIndividualReps
    featureFig = figure; hold on
    if fitModel
        modelFig = figure; hold on
    end
end
if showFrame
    sampleFrameFig = figure;
end
if fitModel
    modelPooledFig = figure; hold on
end
% empty the downsampledVideos directory
if makeDownSampledVid
    delete /Users/sding/Documents/AggScreening/downsampledVideos/*.avi
end

% create legend variable to hold strain name and experiment n numbers
legends = cell(size(strains));
if strcmp(feature,'pcf')
    lineHandles = NaN(numel(strains),1);
end

%% calculate features only if they haven't already been calculated and stored
savedFileName = ['/Users/sding/Documents/AggScreening/results/' feature '_' strainSet '_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel.mat'];
if exist(savedFileName,'file') == 2
    load(savedFileName,featVarName); 
else % calculate features only if saved values don't exist
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
            [strainCtr, fileCtr]
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
            onFoodMask = poly2mask(foodContourCoords(2,:),foodContourCoords(1,:),dims(1),dims(2)); % transposes foodContourCoords to match image coordinate system
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
            if strcmp(feature,'ac')
                endFrameNum = maxLag*frameRate*2 + startFrameNum; % sample twice as many frames to generate maskedImageStack
                if endFrameNum > max(trajData.frame_number)
                    endFrameNum = max(trajData.frame_number);
                end
                numFrames = floor(numel(startFrameNum:endFrameNum)/frameRate/sampleFrameEveryNSec);
                sampleFrames = startFrameNum:sampleFrameEveryNSec*frameRate:endFrameNum;
            else
                endFrameNum = max(trajData.frame_number);
                numFrames = floor(numel(startFrameNum:endFrameNum)/frameRate/sampleFrameEveryNSec);
                sampleFrames = datasample(startFrameNum:endFrameNum,numFrames);
            end
            %% feature-specific set up
            if strcmp(feature,'hc')
                branchHeights.(strain){fileCtr}= cell(numFrames,1); % pre-allocate cells to hold branch height values
            elseif strcmp(feature,'pcf')
                pcf.(strain){fileCtr} = NaN(length(distBins) - 1,numFrames);
            end
            %% go through each frame to generate downsampled frames
            tic
            maskedImageStack = true(numel(1:sampleEveryNPixel:dims(1)),numel(1:sampleEveryNPixel:dims(2)),numFrames);
            if showFrame | makeDownSampledVid
                originalImageStack = NaN(numel(1:sampleEveryNPixel:dims(1)),numel(1:sampleEveryNPixel:dims(2)),numFrames);
            end
            for frameCtr = 1:numFrames
                % load the frame
                imageFrame = h5read(maskedVideoFileName,'/mask',[1,1,double(sampleFrames(frameCtr))],[dims(1),dims(2),1]);
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
                    intensityMaskedImageFrame = imageFrame>0 & imageFrame<70;
                    originalImageStack(:,:,frameCtr) = intensityMaskedImageFrame(1:sampleEveryNPixel:dims(1),1:sampleEveryNPixel:dims(2));
                end
            end
            % generate no movement mask
            movementMask = std(maskedImageStack,0,3)>0;
            % pre-allocate 2D maskedImage stack and frame standard deviation matrix
            if strcmp(feature,'ac')
                maskedImageStack2D = NaN(size(maskedImageStack,1)*size(maskedImageStack,2),size(maskedImageStack,3)); % npixels by time
                numStartingFrames = floor(numFrames/2);
                frameStd = NaN(1,numFrames);
            end
            % go through each frame
            for frameCtr = 1:numFrames
                % apply movement mask
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
                % check that the frame isn't all 1's by preallocation default
                assert (nnz(maskedImageStack(:,:,frameCtr)) < size(maskedImageStack(:,:,frameCtr),1)*size(maskedImageStack(:,:,frameCtr),2),...
                    ['Frame ' num2str(frameCtr) ' has all true pixels by pre-allocation default. Something is wrong'])
                % write frame to 2D maskedImageStack
                if strcmp(feature,'ac')
                    currentImage = maskedImageStack(:,:,frameCtr);
                    maskedImageStack2D(:,frameCtr) = currentImage(:) - mean(currentImage(:));
                    % calculate standard deviation
                    frameStd(frameCtr) = std(maskedImageStack2D(:,frameCtr));
                end
                % calculate feature
                if strcmp(feature,'hc')| strcmp(feature,'pcf')
                    N = nnz(maskedImageStack(:,:,frameCtr));
                    if N>1 % need at least two pixels in frame
                        [x,ydata] = find(maskedImageStack(:,:,frameCtr));
                        pairDists = pdist([x ydata]*sampleEveryNPixel*pixelToMicron/1000); % pairDist in mm;
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
            end
            % pool data from frames
            if strcmp(feature,'hc')
                branchHeights.(strain){fileCtr} = vertcat(branchHeights.(strain){fileCtr}{:});
            end
            % close video
            if makeDownSampledVid
                close(video)
            end
            toc
            tic
            % calculate feature
            if strcmp(feature,'ac')
                ac.(strain){fileCtr} = calculateImageAutocorrelation(maskedImageStack2D,frameStd);
            end
            toc
        end
    end
    %% save variable (only executed if new 'feature_all.mat' has been calculated)
    if saveResults
        if strcmp(strainSet,'all')
            save(['results/' feature '_' strainSet '_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel.mat'],featVarName)%,'-v7.3' for large files for ac feature)
        end
    end
end

%% plot features
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
        % update strain n number for figure legend
        legends{strainCtr} = [strain ', n=' num2str(length(fileInd))];
        % keep track of lag times
        if strcmp(feature,'ac')
            lagTimes(1,fileCtr) = size(ac.(strain){fileCtr},2);
        end
        if plotIndividualReps
            set(0,'CurrentFigure',featureFig)
            if strcmp(feature, 'hc')
                histogram(branchHeights.(strain){fileCtr},'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:),'BinWidth',histBinWidth)
            elseif strcmp(feature,'pcf')
                x = distBins(2:end)-distBinWidth/2; x = x';
                y = nanmean(pcf.(strain){fileCtr},2);
                % plot(x,y,'Color',colorMap(strainCtr,:));
                boundedline(distBins(2:end)-distBinWidth/2,nanmean(pcf.(strain){fileCtr},2),...
                    [nanstd(pcf.(strain){fileCtr},0,2) nanstd(pcf.(strain){fileCtr},0,2)]./sqrt(nnz(sum(~isnan(pcf.(strain){fileCtr}),2))),...
                    'alpha',featureFig.Children,'cmap',colorMap(strainCtr,:))
                % handle visibility off to remove shading to legend
                if fitModel
                    if ~isempty(y) & nnz(isnan(y))~=numel(y)
                        set(0,'CurrentFigure',modelFig)
                        for modelCtr = 1:length(modelFun)
                            [f.(strain){fileCtr},gof.(strain){fileCtr}] = fit(x,y,modelFun{modelCtr}); % fit model function
                            subplot(1,ceil(length(modelFun)/1),modelCtr); hold on
                            plot(x,y,'.','MarkerFaceColor',colorMap(strainCtr,:),'MarkerEdgeColor',colorMap(strainCtr,:));
                            plot(x,f.(strain){fileCtr}(x),'Color',colorMap(strainCtr,:));
                            title({modelFun{modelCtr},formula(f.(strain))})
                        end
                    end
                end
            elseif strcmp(feature, 'ac')
                y = nanmean(ac.(strain){fileCtr},1)';
                x = 0:numel(y)-1; x = x';
                plot(x,y,'Color',colorMap(strainCtr,:))
                if fitModel
                    if ~isempty(y) & nnz(isnan(y))~=numel(x)
                        set(0,'CurrentFigure',modelFig)
                        for modelCtr = 1:length(modelFun)
                            fitOptions.Weights = 1./sqrt(x+1);
                            if strcmp(modelFun{modelCtr},'exp3') % exp3 is a custom model, i.e. not a library function
                                [f.(strain){fileCtr},gof.(strain){fileCtr}] = fit(x,y,exp3,fitOptions); % fit model function
                            elseif strcmp(modelFun{modelCtr},'mexp3')
                                [f.(strain){fileCtr},gof.(strain){fileCtr}] = fit(x,y,mexp3,fitOptions);
                            elseif strcmp(modelFun{modelCtr},'snexp2') 
                                [f.(strain){fileCtr},gof.(strain){fileCtr}] = fit(x,y,snexp2,fitOptions); 
                            else % library functions
                                [f.(strain){fileCtr},gof.(strain){fileCtr}] = fit(x,y,modelFun{modelCtr});
                            end
                            subplot(1,ceil(length(modelFun)/1),modelCtr); hold on
                            plot(x,y,'.','MarkerFaceColor',colorMap(strainCtr,:),'MarkerEdgeColor',colorMap(strainCtr,:));
                            plot(x,f.(strain){fileCtr}(x),'Color',colorMap(strainCtr,:));
                            title(modelFun{modelCtr})
                        end
                    end
                end
            end
        end
    end
    % pool data from multiple files
    if strcmp(feature,'hc')
        branchHeights_pooled.(strain) = vertcat(branchHeights.(strain){:});
    elseif strcmp(feature,'pcf')
        pcf_pooled.(strain) = horzcat(pcf.(strain){:});
    elseif strcmp(feature,'ac')
        ac_pooled.(strain) = vertcat(ac.(strain){:});
    end
    % plot features
    set(0,'CurrentFigure',featurePooledFig)
    if strcmp(feature,'hc')
        histogram(branchHeights_pooled.(strain),'Normalization',histogramNormalisation,'DisplayStyle','stairs','EdgeColor',colorMap(strainCtr,:),'BinWidth',histBinWidth)
    elseif strcmp(feature,'pcf')
        x = distBins(2:end)-distBinWidth/2; x = x';
        y = nanmean(pcf_pooled.(strain),2);
        % plot(x,y,'Color',colorMap(strainCtr,:));
        [lineHandles(strainCtr), ~] = boundedline(distBins(2:end)-distBinWidth/2,nanmean(pcf_pooled.(strain),2),...
            [nanstd(pcf_pooled.(strain),0,2) nanstd(pcf_pooled.(strain),0,2)]./sqrt(nnz(sum(~isnan(pcf_pooled.(strain)),2))),...
            'alpha',featurePooledFig.Children,'cmap',colorMap(strainCtr,:));
        if fitModel
            set(0,'CurrentFigure',modelPooledFig)
            for modelCtr = 1:length(modelFun)
                [f_pooled.(strain),gof_pooled.(strain)] = fit(x,y,modelFun{modelCtr}); % fit model function
                subplot(1,ceil(length(modelFun)/1),modelCtr); hold on
                plot(x,y,'.','MarkerFaceColor',colorMap(strainCtr,:),'MarkerEdgeColor',colorMap(strainCtr,:));
                plot(x,f_pooled.(strain)(x),'Color',colorMap(strainCtr,:));
                title({modelFun{modelCtr},formula(f_pooled.(strain))})
            end
        end
    elseif strcmp(feature,'ac')
        y = nanmean(ac_pooled.(strain),1)';
        x = 0:numel(y)-1; x = x';
        plot(x,y,'Color',colorMap(strainCtr,:));
        if fitModel
            set(0,'CurrentFigure',modelPooledFig)
            fitOptions.Weights = 1./sqrt(x+1);
            for modelCtr = 1:length(modelFun)
                if strcmp(modelFun{modelCtr},'exp3') % exp3 is a custom model, i.e. not a library function
                    [f_pooled.(strain),gof_pooled.(strain)] = fit(x,y,exp3,fitOptions); % fit model function
                elseif strcmp(modelFun{modelCtr},'mexp3') 
                    [f_pooled.(strain),gof_pooled.(strain)] = fit(x,y,mexp3,fitOptions); 
                elseif strcmp(modelFun{modelCtr},'snexp2') 
                    [f_pooled.(strain),gof_pooled.(strain)] = fit(x,y,snexp2,fitOptions); 
                else % library functions
                    [f_pooled.(strain),gof_pooled.(strain)] = fit(x,y,modelFun{modelCtr});
                end
                subplot(1,ceil(length(modelFun)/1),modelCtr); hold on
                plot(x,y,'.','MarkerFaceColor',colorMap(strainCtr,:),'MarkerEdgeColor',colorMap(strainCtr,:));
                plot(x,f_pooled.(strain)(x),'Color',colorMap(strainCtr,:));
                title({modelFun{modelCtr},formula(f_pooled.(strain))})
            end
        end
    end
end

%% save model fitting variable
%if saveResults
    if fitModel
        if length(modelFun) == 1
            save(['results/modelFit/' feature '_' strainSet '_fitmodel_' modelFun{1} '_' featVarName '_Pooled_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel.mat'],'f_pooled','gof_pooled')
            if plotIndividualReps
                save(['results/modelFit/' feature '_' strainSet '_fitmodel_' modelFun{1} '_' featVarName '_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel.mat'],'f','gof')
            end
        end
    end
%end

%% format and export
set(featurePooledFig,'PaperUnits','centimeters')
set(0,'CurrentFigure',featurePooledFig)
if strcmp(feature,'pcf')
    % legend(featurePooledFig.Children,lineHandles,legends,'Location','eastoutside') % boundedline plot
    legend(legends,'Location','eastoutside') % model fitting plot
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

if fitModel
    set(modelPooledFig,'PaperUnits','centimeters')
    set(0,'CurrentFigure',modelPooledFig)
    xlabel(xLabel)
    ylabel(yLabel)
    if saveResults
        if length(modelFun) == 1
            figurename = ['figures/modelFit/' feature '_' strainSet '_fitmodel_' modelFun{1} '_' featVarName '_Pooled_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel'];
        else
            figurename = ['figures/modelFit/' feature '_' strainSet '_fitmodel' num2str(length(modelFun)) '_' featVarName '_Pooled_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel'];
        end
        exportOptions2 = struct('Format','eps2',...
            'Color','rgb',...
            'Width',30,...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',20,...
            'LineWidth',2);
        exportfig(modelPooledFig,[figurename '.eps'],exportOptions2)
    end
    if plotIndividualReps
        set(modelFig,'PaperUnits','centimeters')
        set(0,'CurrentFigure',modelFig)
        xlabel(xLabel)
        ylabel(yLabel)
        if saveResults
            if length(modelFun) == 1
                figurename = ['figures/modelFit/' feature '_' strainSet '_fitmodel_' modelFun{1} '_' featVarName '_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel'];
            else
                figurename = ['figures/modelFit/' feature '_' strainSet '_fitmodel' num2str(length(modelFun)) '_' featVarName '_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel'];
            end
            exportfig(modelFig,[figurename '.eps'],exportOptions2)
        end
    end
end