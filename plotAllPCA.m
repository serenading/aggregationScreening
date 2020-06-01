%% script performs PCA analysis with all features from all recordings (both 40 and 5 worms) in the full dataset, and plots in shared PC space
% author: @serenading May 2020

clear
close all

addpath('auxiliary/')

%% Set parameters
featExtractTimestamp = '20200519_153722'; %'20200519_153722','20191024_122847' or '20181203_141111'
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features. =17
feats2drop = {'path'};
strains2drop = {};

%% Load features table and extract features matrix
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' featExtractTimestamp '.csv'],'Delimiter',',','preserveVariableNames',true);
% Drop features and strains as specified
[featureTable, ~] = dropFeats(featureTable,feats2drop);
[featureTable,~] = dropStrains(featureTable,strains2drop);
featureMat = table2array(featureTable(:, n_nonFeatVar+1:end));
n_strains = numel(unique(featureTable.strain_name));
n_feats = size(featureMat,2);

%% Analyze features with PCA

% pre-process feature matrix for PCA
[featureMat,~] = preprocessFeatMat(featureMat);
% do pca
[pc, score, ~, ~, explained] = pca(featureMat);

%% Plot first two PCs and colour by different variables of interest

% general PCA plot
generalFig = figure; plot(score(:,1),score(:,2),'.')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
title(['PCA plot with 198 strains and ' num2str(size(featureMat,2)) ' feature set'])

% by worm number
wormNumFig = figure; hold on
fortyLogInd = featureTable.wormNum==40;
fiveLogInd = featureTable.wormNum==5;
plot(score(fortyLogInd,1),score(fortyLogInd,2),'r.')
plot(score(fiveLogInd,1),score(fiveLogInd,2),'b.')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
legend({'40 worms','5 worms'})
title(['PCA plot with 198 strains and ' num2str(size(featureMat,2)) ' feature set'])

% by strain
strainFig = figure; hold on
N2LogInd = strcmp(featureTable.strain_name,'N2');
CB4856LogInd = strcmp(featureTable.strain_name,'CB4856');
DA609LogInd = strcmp(featureTable.strain_name,'DA609');
otherStrainsLogInd = ~N2LogInd & ~CB4856LogInd & ~DA609LogInd;
plot(score(otherStrainsLogInd,1),score(otherStrainsLogInd,2),'y.')
plot(score(N2LogInd,1),score(N2LogInd,2),'c.')
plot(score(CB4856LogInd,1),score(CB4856LogInd,2),'r.')
plot(score(DA609LogInd,1),score(DA609LogInd,2),'b.')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
legend({'others','N2','CB4856','DA609',})
title(['PCA plot with 198 strains and ' num2str(size(featureMat,2)) ' feature set'])

%% 3D plot of the first three PCs
% by worm number
figure; 
scatter3(score(fortyLogInd,1),score(fortyLogInd,2),score(fortyLogInd,3),'r.')
hold on
scatter3(score(fiveLogInd,1),score(fiveLogInd,2),score(fiveLogInd,3),'b.')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
zlabel(['PC3 (' num2str(round(explained(3))) ')%'])
legend({'40 worms','5 worms'})
title(['PCA plot with ' num2str(n_strains) ' strains and ' num2str(n_feats) ' feature set'])

% by strain
figure;
scatter3(score(otherStrainsLogInd,1),score(otherStrainsLogInd,2),score(otherStrainsLogInd,3),'y.')
hold on
scatter3(score(N2LogInd,1),score(N2LogInd,2),score(N2LogInd,3),'c.')
scatter3(score(CB4856LogInd,1),score(CB4856LogInd,2),score(CB4856LogInd,3),'r.')
scatter3(score(DA609LogInd,1),score(DA609LogInd,2),score(DA609LogInd,3),'b.')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
zlabel(['PC3 (' num2str(round(explained(3))) ')%'])
legend({'others','N2','CB4856','DA609',})
title(['PCA plot with ' num2str(n_strains) ' strains and ' num2str(n_feats) ' feature set'])