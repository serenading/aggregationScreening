clear
close all

%% script takes 5 worm imaging datasets for 3 control strains and plots them in shared PC space
% to assess the effect of diapause length, bleach prep, recording time of the day (approx. by run number), and camera number on data
% author: @serenading Oct 2019

%% set parameters
% set analysis parameters
strains = {'N2','DA609','CB4856'};
markerShapes = {'.','*','o'}; % marker shape to differentiate three strains in combined plots
wormNum = 5; % dataset exists for 5 and 40 worms.
useTierpsy256 = true; % true to use 256 features, false to use all (>4000) features
dropFeatThreshold = 0.2; % the maximum fraction of NaN values that a feature can have before being dropped
% imputeNaNThreshold = 0; 
drawPolygon = true; % connect related datapoints for easy visualisation

%% prep work
addpath('auxiliary/')
% load metadata
metadata = readtable('/Volumes/behavgenom_archive$/Serena/AggregationScreening/metadata_aggregationScreening.csv');
% load Tierpsy feature summary files and corresponding filenames
features_summary = readtable('/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/features_summary_tierpsy_plate_20191024_122847.csv');
filenames_summary = readtable('/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/filenames_summary_tierpsy_plate_20191024_122847.csv','Format', '%d%s%s');
% join the Tierpsy tables to match filenames with file_id. Required in case features were not extracted for any files.
combinedTierpsyTable = outerjoin(filenames_summary, features_summary,'MergeKeys', true);

%% use Tierpsy 256 features
if useTierpsy256
    % load the set of 256 features selected using based on classification accuracy on a set of mutant strains
    top256 = readtable('./auxiliary/tierpsy_256.csv','ReadVariableNames',0);
    featNames = top256.Var1;
    % curtail featNames to max 63 characters; anything beyond this is truncated inside combinedTierpsyTable and will not match up
    for featCtr = 1:numel(featNames)
        featNameLength = numel(featNames{featCtr});
        if featNameLength > 63
            featNames{featCtr} = featNames{featCtr}(1:63);
        end
    end
    % trim down feature matrix to contain just 256 features
    featMat = combinedTierpsyTable{:, featNames}; % this gives error when feat names do not match up perfectly
end

%% get info from metadata for files as specified by strain and worm number
% find logical indices for valid files
fileLogInd = [];
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    fileLogInd = [fileLogInd strcmp(metadata.strain_name,strain) & metadata.wormNum == wormNum & metadata.is_bad ==0];
end
fileLogIndAllStrains = logical(sum(fileLogInd,2)); % combine valid indices for all strains
numFiles = nnz(fileLogIndAllStrains);
% extract file name parts
dirname = metadata.dirname(fileLogIndAllStrains);
basename = metadata.basename(fileLogIndAllStrains);
strainname = metadata.strain_name(fileLogIndAllStrains);
% get days in diapause and run and bleach numbers for assessment
daysDiapause = metadata.daysDiapause(fileLogIndAllStrains);
runNum = metadata.run(fileLogIndAllStrains);
bleachNum = metadata.block(fileLogIndAllStrains);
camNum = metadata.channel(fileLogIndAllStrains);

%% go through each file name to grow combined features file indices
fileInd = [];
for fileCtr = 1:numFiles
    % get features file name
    filename = strrep(strrep(strcat(dirname{fileCtr},'/',basename{fileCtr}),'MaskedVideos','Results'),'.hdf5','_featuresN.hdf5');
    % get the file_id for this file (this is the "file_id" inside the "features_summary" and "filenames_summary" files and uses Python indexing)
    fileIdx = find(strcmp(combinedTierpsyTable.file_name,filename));
    % add to list of file indices
    fileInd = vertcat(fileInd,fileIdx);
end
assert(numel(fileInd) == numFiles)

%% get neccesary features for all files of interest
if useTierpsy256
    featMat = featMat(fileInd,:); % 256 features
else
    featMat = table2array(combinedTierpsyTable(fileInd,4:end)); % all features
end
lengths = combinedTierpsyTable.length_50th(fileInd); % length feature

%% analyze features
% drop features with too many NaN's
numNanFeat = sum(isnan(featMat),1);
nanFeatInd = find(numNanFeat>0); % get indices for features with NaN's
dropCtr = 0;
for featCtr = numel(nanFeatInd):-1:1
    featIdx = nanFeatInd(featCtr); % get feature index
    if numNanFeat(featIdx)>numFiles*dropFeatThreshold % if there are too many NaN values for this feature
        featMat = [featMat(:,1:nanFeatInd(featCtr)-1), featMat(:,nanFeatInd(featCtr)+1:end)]; % drop feature by concatenating matrix on either side of this feature
        dropCtr = dropCtr+1;
    end
end
disp([ num2str(dropCtr) ' out of  ' num2str(size(featMat,2)) ' features dropped due to too many NaN values'])

% impute nan values to global mean
featMeans = nanmean(featMat);
imputeCtr = 0;
for featCtr = 1:size(featMat, 2)
    nanInds = isnan(featMat(:, featCtr));
    if nnz(nanInds)>0
        featMat(nanInds, featCtr) = featMeans(featCtr);
        imputeCtr = imputeCtr+1;
    end
end
disp([ num2str(imputeCtr) ' out of  ' num2str(size(featMat,2)) ' features have NaN values imputed'])

% drop features with zero standard deviation
featStds = std(featMat);
zeroFeatStdInd = find(featStds == 0); % get indices for features with zero standard deviation
dropCtr = 0;
for featCtr = numel(zeroFeatStdInd):-1:1
    featMat = [featMat(:,1:zeroFeatStdInd(featCtr)-1), featMat(:,zeroFeatStdInd(featCtr)+1:end)]; % drop feature by concatenating matrix on either side of this feature
    dropCtr = dropCtr+1;
end
disp([ num2str(dropCtr) ' out of  ' num2str(size(featMat,2)) ' features dropped due to zero standard deviation'])

% z-normalise feature matrix
featMatNorm = normalize(featMat,1);

% do pca
[pc, score, ~, ~, explained] = pca(featMatNorm,'NumComponents',5);

%% plot first two PCs for a number of experimental variables to assess their effect on the data

% days in diapause plot: 2x2 plot, first panel all strains
diapauseFigure = figure;
subplot(2,2,1); hold on
days = unique(daysDiapause);
colorMap = hsv(numel(days));

for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    strainLogInd = strcmp(strainname, strain);
    for dayCtr = 1:numel(days)
        day = days(dayCtr);
        dayInd = daysDiapause == day;
        % combine strain and days
        plotLogInd = dayInd & strainLogInd;
        x = score(plotLogInd, 1);
        y = score(plotLogInd, 2);
        for plotCtr = [1, strainCtr+1]
            subplot(2,2,plotCtr); hold on
            plot(x, y, markerShapes{strainCtr}, 'MarkerSize', 12, 'Color', colorMap(dayCtr,:))
        end
        legends{dayCtr} = ['day ' num2str(day) ', n=' num2str(numel(x))];
    end
    legend(legends)
    title ([strain ', days in diapause'])
    xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
    ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
end
subplot(2,2,1)
title (['days in diapause'])
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])

% run number plot (to approximate time of the day/worm age): 2x2 plot, first panel all strains
runFig = figure; hold on
subplot(2,2,1); hold on
runs = unique(runNum);
colorMap = hsv(numel(runs));

for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    strainLogInd = strcmp(strainname, strain);
    skippedRunInd = ones(numel(runs),1); % keep track of run numbers that don't exist for a certain strain
    clear plotHandle
    for runCtr = 1:numel(runs)
        run = runs(runCtr);
        runInd = runNum == run;
        % combine strain and runs
        plotLogInd = runInd & strainLogInd;
        if nnz(plotLogInd)>0
            skippedRunInd(runCtr) = 0;
            x = score(plotLogInd, 1);
            y = score(plotLogInd, 2);
            for plotCtr = [1, strainCtr+1]
                subplot(2,2,plotCtr); hold on
                plotHandle(runCtr) = plot(x, y, markerShapes{strainCtr}, 'MarkerSize', 12, 'Color', colorMap(runCtr,:),...
                    'DisplayName',['run ' num2str(run) ', n=' num2str(numel(x))]);
                if drawPolygon & nnz(plotLogInd)>2 % need more than 2 points to draw convex hull shape
                    polygon = convhull(x, y);
                    patch(x(polygon), y(polygon), colorMap(runCtr,:),'FaceAlpha',0.3,'EdgeColor',colorMap(runCtr,:))
                end
            end
        end
    end
    legend(plotHandle(~skippedRunInd))
    title ([strain ', run number'])
    xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
    ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
end
subplot(2,2,1); hold on
title ('run number')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])

% bleach prep: 2x2 plot, first panel all strains
bleachFig = figure; hold on
subplot(2,2,1); hold on
bleaches = unique(bleachNum);
colorMap = hsv(numel(bleaches));

for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    strainLogInd = strcmp(strainname, strain);
    clear plotHandle
    for bleachCtr = 1:numel(bleaches)
        bleach = bleaches(bleachCtr);
        bleachInd = bleachNum == bleach;
        % combine strain and bleaches
        plotLogInd = bleachInd & strainLogInd;
        x = score(plotLogInd, 1);
        y = score(plotLogInd, 2);
        for plotCtr = [1, strainCtr+1]
            subplot(2,2,plotCtr); hold on
            plotHandle(bleachCtr) = plot(x, y, markerShapes{strainCtr}, 'MarkerSize', 12, 'Color', colorMap(bleachCtr,:),...
                'DisplayName',['bleach prep ' num2str(bleach) ', n=' num2str(numel(x))]);
            if drawPolygon & nnz(plotLogInd)>2 % need more than 2 points to draw convex hull shape
                polygon = convhull(x, y);
                patch(x(polygon), y(polygon), colorMap(bleachCtr,:),'FaceAlpha',0.3,'EdgeColor',colorMap(bleachCtr,:))
            end
        end
    end
    legend(plotHandle,'Location','eastoutside')
    title ([strain ', bleach prep'])
    xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
    ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
end
set(0,'CurrentFigure',bleachFig)
subplot(2,2,1);
title ([strain ', bleach prep'])
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])

% camera number plot: 2x2 plot, first panel all strains
camFig = figure; hold on
subplot(2,2,1); hold on
cams = unique(camNum);
colorMap = hsv(numel(cams));

for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    strainLogInd = strcmp(strainname, strain);
    skippedCamInd = ones(numel(cams),1); % keep track of camera numbers that don't exist for a certain strain
    clear plotHandle
    for camCtr = 1:numel(cams)
        cam = cams(camCtr);
        camInd = camNum == cam;
        % combine strain and cameras
        plotLogInd = camInd & strainLogInd;
        if nnz(plotLogInd)>0
            skippedCamInd(camCtr) = 0;
            x = score(plotLogInd, 1);
            y = score(plotLogInd, 2);
            for plotCtr = [1, strainCtr+1]
                subplot(2,2,plotCtr); hold on
                plotHandle(camCtr) = plot(x, y, markerShapes{strainCtr}, 'MarkerSize', 12, 'Color', colorMap(camCtr,:),...
                    'DisplayName',['cam ' num2str(cam) ', n=' num2str(numel(x))]);
                if drawPolygon & nnz(plotLogInd)>2 % need more than 2 points to draw convex hull shape
                    polygon = convhull(x, y);
                    patch(x(polygon), y(polygon), colorMap(camCtr,:),'FaceAlpha',0.3,'EdgeColor',colorMap(camCtr,:))
                end
            end
        end
    end
    legend(plotHandle(~skippedCamInd))
    title ([strain ', camera number'])
    xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
    ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
end
subplot(2,2,1); hold on
title ('camera number')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])

% % check length as a function of time of the day
% lengthsFig = figure; hold on
% for strainCtr = 1:numel(strains)
%     strain = strains{strainCtr};
%     strainLogInd = strcmp(strainname, strain);
%     group = NaN(numel(strainLogInd),1);
%     for runCtr = 1:numel(runs)
%         run = runs(runCtr);
%         runInd = runNum == run;
%         plotLogInd = runInd & strainLogInd;
%         group(plotLogInd) = run;
%     end
%     set(0,'CurrentFigure',lengthsFig)
%     subplot(1,3,strainCtr);
%     boxplot(lengths(strainLogInd),group(strainLogInd))
%     ylim([600 1050])
%     title (strain)
%     xlabel('run number')
%     ylabel('median length (microns)')
% end