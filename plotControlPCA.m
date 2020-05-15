%% script takes 5 worm imaging datasets for 3 control strains and plots them in shared PC space
% to assess the effect of strain, diapause length, bleach prep, recording time of the day (approx. by run number), and camera number on data
% author: @serenading May 2020

clear
close all

addpath('auxiliary/')

%% Set parameters

% which strains and worm number
strains = {'N2','DA609','CB4856'};
wormNums = 5; % 5, or, 40, or [5, 40]

% which features table to load
featExtractTimestamp = '20191024_122847'; %'20191024_122847' or '20181203_141111'
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features. =17

% which features to use for PCA
useTierpsy256 = false; % true to use 256 features, false to use all (>4000) features

% how to plot
markerShapes = {'.','*','o'}; % marker shape to differentiate three strains in combined plots
drawPolygon = true; % connect related datapoints for easy visualisation

% how many PC's to use for manova test
useVariablePCs = false; % True: use hand-determined number of PC's for each test; False: use number of PC's necessary to explain a variance threshold set by minVar4PC
minVarianceExplained = 60; % minium variance threshold (in %) to determine the number of PC needed.

%% Load features table
featureTable = readtable(['/Users/sding/Dropbox/aggScreening/results/fullFeaturesTable_' featExtractTimestamp '.csv'],'Delimiter',',','preserveVariableNames',true);

%% Get rows corresponding to specified strain and worm number
% compile logical index by strain
combinedStrainLogInd = false(size(featureTable,1),1);
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    strainLogInd = strcmp(featureTable.strain_name,strain);
    combinedStrainLogInd(strainLogInd) = true;
end
% compile logical index by worm number
combinedWormNumLogInd = false(size(featureTable,1),1);
for wormNumCtr = 1:numel(wormNums)
    wormNum = wormNums(wormNumCtr);
    wormNumLogInd = featureTable.wormNum == wormNum;
    combinedWormNumLogInd(wormNumLogInd) = true;
end
% combine to get row logical index by strain and number
rowLogInd = combinedStrainLogInd & combinedWormNumLogInd;

%% Get features matrix from features table for selected files (rows) and features (columns)
% optionally use the Tierpsy256 feature set
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
    featureMat = featureTable{rowLogInd, featNames};
else
    % use all features for feature matrix
    featureMat = table2array(featureTable(rowLogInd, n_nonFeatVar+1:end));
end

% get length feature separately
lengths = featureTable.length_50th(rowLogInd); 

%% analyze features with PCA

% pre-process feature matrix for PCA
featureMat = preprocess4PCA(featureMat);

% do pca
[pc, score, ~, ~, explained] = pca(featureMat);

% determine how many PC's are needed to explain a specied amount of variance in the data
explainedCumSum = cumsum(explained);
minNumPC = find(explainedCumSum > minVarianceExplained,1); % the number of PC needed
disp(['The first ' num2str(minNumPC) ' PC explain at least ' num2str(minVarianceExplained) '% of the variance in the data'])

% % see what's inside the first PC
% [feat,featInd] = sort(pc(:,1)); % PC1 
% if useTierpsy256
%     featNames(featInd)
% else
%     featureTable.Properties.VariableNames(featInd)'
% end

%% Get information from features table for the relevant files
% get days in diapause and run and bleach and camera numbers for assessment
strainName = featureTable.strain_name(rowLogInd);
daysDiapause = featureTable.daysDiapause(rowLogInd);
runNum = featureTable.run(rowLogInd);
bleachNum = featureTable.block(rowLogInd);
camNum = featureTable.channel(rowLogInd);

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