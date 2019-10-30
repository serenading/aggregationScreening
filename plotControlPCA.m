clear
close all


%% script takes 5 worm imaging datasets for 3 control strains and plots them in shared PC space
% to assess the effect of strain, diapause length, bleach prep, recording time of the day (approx. by run number), and camera number on data
% author: @serenading Oct 2019

%% set parameters

% which strains and worm number
strains = {'N2','DA609','CB4856'};
wormNum = 5; % dataset exists for 5 and 40 worms.

% which features to use for PCA
useTierpsy256 = false; % true to use 256 features, false to use all (>4000) features
dropFeatThreshold = 0.2; % the maximum fraction of NaN values that a feature can have before being dropped

% how to plot
markerShapes = {'.','*','o'}; % marker shape to differentiate three strains in combined plots
drawPolygon = true; % connect related datapoints for easy visualisation

% how many PC's to use for manova test
useVariablePCs = false; % True: use hand-determined number of PC's for each test; False: use number of PC's necessary to explain a variance threshold set by minVar4PC
minVarianceExplained = 60; % minium variance threshold (in %) to determine the number of PC needed.

%% prep work

addpath('auxiliary/')
% load metadata
metadata = readtable('/Volumes/behavgenom_archive$/Serena/AggregationScreening/metadata_aggregationScreening.csv');
% load Tierpsy feature summary files and corresponding filenames
features_summary = readtable('/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/features_summary_tierpsy_plate_20191024_122847.csv');
filenames_summary = readtable('/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/filenames_summary_tierpsy_plate_20191024_122847.csv','Format', '%d%s%s');

% join the Tierpsy tables to match filenames with file_id. Required in case features were not extracted for any files.
combinedTierpsyTable = outerjoin(filenames_summary, features_summary,'MergeKeys', true);

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

% get days in diapause and run and bleach and camera numbers for assessment
strainName = metadata.strain_name(fileLogIndAllStrains);
daysDiapause = metadata.daysDiapause(fileLogIndAllStrains);
runNum = metadata.run(fileLogIndAllStrains);
bleachNum = metadata.block(fileLogIndAllStrains);
camNum = metadata.channel(fileLogIndAllStrains);

%% get features file indices as specified by strain and worm number

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

if useTierpsy256 % use Tierpsy256 feature set
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
    featMat = combinedTierpsyTable{fileInd, featNames}; 

else % use all features
    featMat = table2array(combinedTierpsyTable(fileInd,4:end)); 
end

% get length feature
lengths = combinedTierpsyTable.length_50th(fileInd); 

%% analyze features with PCA

% drop features with too many NaN's
featBefore = size(featMat,2);
numNanFeat = sum(isnan(featMat),1);
colsToDrop = numNanFeat > numFiles*dropFeatThreshold; % get logical index for features with too many NaN's
featMat = featMat(:,~colsToDrop);
disp([ num2str(nnz(colsToDrop)) ' out of  ' num2str(featBefore) ' features dropped due to too many NaN values'])

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
featBefore = size(featMat,2);
featStds = std(featMat);
colsToDrop = featStds == 0 ; % get logical index for features with too many NaN's
featMat = featMat(:,~colsToDrop);
disp([ num2str(nnz(colsToDrop)) ' out of  ' num2str(featBefore) ' features dropped due to zero standard deviation'])

% z-normalise feature matrix
featMatNorm = normalize(featMat,1);

% do pca
[pc, score, ~, ~, explained] = pca(featMatNorm);

% determine how many PC's are needed to explain a specied amount of variance in the data
explainedCumSum = cumsum(explained);
minNumPC = find(explainedCumSum > minVarianceExplained,1); % the number of PC needed
disp(['The first ' num2str(minNumPC) ' PC explain at least ' num2str(minVarianceExplained) '% of the variance in the data'])

%% plot first two PCs for a number of experimental variables to assess their effect on the data

% days in diapause plot: 2x2 plot, first panel all strains
diapauseFigure = figure;
subplot(2,2,1); hold on
days = unique(daysDiapause);
colorMap = hsv(numel(days));

for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    strainLogInd = strcmp(strainName, strain);
    for dayCtr = 1:numel(days)
        day = days(dayCtr);
        dayInd = daysDiapause == day;
        % combine strain and days
        plotLogInd = dayInd & strainLogInd;
        x = score(plotLogInd, 1);
        y = score(plotLogInd, 2);
        for plotCtr = [1, strainCtr+1]
            subplot(2,2,plotCtr); hold on
            plotHandle(dayCtr) = plot(x, y, markerShapes{strainCtr}, 'MarkerSize', 12, 'Color', colorMap(dayCtr,:),...
                'DisplayName',['day ' num2str(day) ', n=' num2str(numel(x))]);
            if drawPolygon & nnz(plotLogInd)>2 % need more than 2 points to draw convex hull shape
                polygon = convhull(x, y);
                patch(x(polygon), y(polygon), colorMap(dayCtr,:),'FaceAlpha',0.3,'EdgeColor',colorMap(dayCtr,:))
            end
        end
    end
    legend(plotHandle)
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
    strainLogInd = strcmp(strainName, strain);
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
    strainLogInd = strcmp(strainName, strain);
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
    strainLogInd = strcmp(strainName, strain);
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

%% perform manova using PC values as variables (note: the number of observations per group needs to exceed the number of variables)

% decide how many PCs to use as variables for manova
if useVariablePCs % individually decide for each test based on the number of available observations
    numPCs = [50, 40, 10, 7, 40]; 
else % use the number of PC's needed to explain a certain amount of variance
    numPCs = minNumPC*ones(1,5); % generate numPcs value set for all 5 tests
end
pcs2use = containers.Map({'strain','diapause','run','bleach','camera'},numPCs);

% strains (use all strains with minimum 50 observations)
[d,p] = manova1(score(:,1:pcs2use('strain')),strainName);
if d>0
    disp('Significant difference grouping by strain')
    d, p
end
% days in diapause (only use days 1-3 with minimum 42 observations)
manLogInd = daysDiapause<=3;
[d,p] = manova1(score(manLogInd,1:pcs2use('diapause')),daysDiapause(manLogInd));
if d>0
    disp('Significant difference grouping by days in diapause')
    d, p
end
% run number (only use runs 1-6 with 12 minimum observations)
manLogInd = runNum<=6;
[d,p] = manova1(score(manLogInd,1:pcs2use('run')),runNum(manLogInd));
if d>0
    disp('Significant difference grouping by run number')
    d, p
end
% bleach number (only use bleaches 1-15 with minimum 8 observations)
manLogInd = bleachNum<=15;
[d,p] = manova1(score(manLogInd,1:pcs2use('bleach')),bleachNum(manLogInd));
if d>0
    disp('Significant difference grouping by bleach number')
    d, p
end
% camera number (only use cameras 2,4,6 with minimum 42 observations)
manLogInd = camNum==2|camNum==4|camNum==6;
[d,p] = manova1(score(manLogInd,1:pcs2use('camera')),camNum(manLogInd));
if d>0
    disp('Significant difference grouping by camera number')
    d, p
end