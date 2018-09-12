clear
close all

%% script extracts features for each of the strains of interest based on downsampled pixel data.

%% set parameters
% set analysis parameters
strainSet = 'controls'; % 'controls','divergent','all'
feature = 'hc'; % specify feature as string. 'pcf' (pair correlation function), 'hc'(hierarchical clustering), 'gf'(giant fluctuation).
maxNumReplicates =5; % controls have up to 60 reps, divergents up to 15 reps, all other strains up to 5 reps.
sampleFrameEveryNSec = 5;
sampleEveryNPixel = 8; % 8 or 16
saveResults = false;
makeDownSampledVid = false;
plotIndividualReps = false;

% set default parameters
phaseRestrict = true; % phaseRestrict cuts out the first 15 min of each video
histogramNormalisation = 'pdf'; % 'pdf' by default. 'count' an option
dims = [2048 2048];
pixelToMicron = 10; % 10 microns per pixel, read by pixelsize = double(h5readatt(filename,'/trajectories_data','microns_per_pixel?))

% set feature-specific parameters
if strcmp(feature,'hc')
    linkageMethod = 'single'; % 'single' (preferred), 'average','centroid','complete','median','weighted'
elseif strcmp(feature,'pcf')
elseif strcmp(feature,'gf')
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
featureFig = figure; hold on
sampleFrameFig = figure;
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
    % update strain n number for figure legend
    legends{strainCtr} = [strain ', n=' num2str(length(fileInd))];
    %% go through each recording
    for fileCtr = 1:length(fileInd)
        %% load data
        filename = filenames{fileInd(fileCtr)};
        trajData = h5read(filename,'/trajectories_data');
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        maskedVideoFileName = strrep(strrep(filename,'Results','MaskedVideos'),'_skeletons.hdf5','.hdf5');
        %fileInfo = h5info(maskedVideoFileName); % slow
        %dims = fileInfo.Datasets(2).Dataspace.Size; %[2048,2048,num]
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
        end
        %% go through each frame to generate a 3D stack of downsampled frames
        % preallocate 3D image stack
        frameStack = NaN(numel(1:sampleEveryNPixel:dims(1)), numel(1:sampleEveryNPixel:dims(2)), numFrames);
        for frameCtr = 1:numFrames
            % load the frame
            imageFrame = h5read(maskedVideoFileName,'/mask',[1,1,double(sampleFrames(frameCtr))],[dims(1),dims(2),1]);
            % generate binary segmentation based on black/white contrast
            binaryImage = imageFrame>0 & imageFrame<70;
            % downsample image
            downsampleBinaryImage = binaryImage(1:sampleEveryNPixel:dims(1),1:sampleEveryNPixel:dims(2));
            % add downsampled image to the stack
            frameStack(:,:,frameCtr) = downsampleBinaryImage;
        end
        % determine pixels that don't move from frame to frame and exclude those as artefacts
        stackStd = std(frameStack,0,3); % calculate standard deviation along the temporal/third dimension
        noMoveMask = stackStd>0; % get mask for moving pixels
        for frameCtr = 1:numFrames
            frameStack(:,:,frameCtr) = frameStack(:,:,frameCtr) & noMoveMask; % set nonmoving pixels to zero
            numPixelsRemoved = nnz(frameStack(:,:,frameCtr))-nnz(frameStack(:,:,frameCtr) & noMoveMask);
            if numPixelsRemoved >0
                disp(['frame ' num2str(frameCtr) ' has ' numPixelsRemoved ' pixels removed'])
            end
        end
        % calculate feature
        for frameCtr = 1:numFrames
            set(0,'CurrentFigure',sampleFrameFig)
            imshow(frameStack(:,:,frameCtr))
            if strcmp(feature,'hc')
                N = nnz(frameStack(:,:,frameCtr));
                if N>1 % need at least two worms in frame
                    [x,y] = find(frameStack(:,:,frameCtr));
                    pairDists = pdist([x y]*sampleEveryNPixel*pixelToMicron/1000); % pairDist in mm;
                    clustTree = linkage(pairDists,linkageMethod);
                    branchHeights.(strain){fileCtr}{frameCtr} = clustTree(:,3);
                    % dendrogram(clustTree,0,'Reorder',optimalleaforder(clustTree,pairDists));
                end
            end
        end
        % pool data from frames
        branchHeights.(strain){fileCtr} = vertcat(branchHeights.(strain){fileCtr}{:});
    end
    % pool data from multiple files
    branchHeights.(strain) = vertcat(branchHeights.(strain){:});
   % plot histogram of branch heights
    set(0,'CurrentFigure',featureFig)
    histogram(branchHeights.(strain),'Normalization',histogramNormalisation,'EdgeColor',colorMap(strainCtr,:),'DisplayStyle','stairs','BinWidth',0.2)
end

%% save variable
if saveResults
    save(['results/' feature '_' strainSet '_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel.mat'],'branchHeights')
end

% format and export
set(featureFig,'PaperUnits','centimeters')
set(0,'CurrentFigure',featureFig)
legend(legends,'Location','northeast')
set(gca, 'YScale', 'log')
xlim([0 10])
ylabel('probability')
xlabel('inter-wormpixel distance (mm)')
title([linkageMethod ' linkage'],'FontWeight','normal')
if saveResults
    figurename = ['figures/' feature '_' strainSet '_branchHeights_' linkageMethod '_sample' num2str(sampleFrameEveryNSec) 's_' num2str(sampleEveryNPixel) 'pixel'];
    exportfig(featureFig,[figurename '.eps'],exportOptions)
end