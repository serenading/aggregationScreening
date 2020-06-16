%% script performs PCA analysis with all features from all recordings (both 40 and 5 worms) in the full dataset, and plots in shared PC space
% author: @serenading May 2020

clear
close all

addpath('auxiliary/')

%% Set parameters
featExtractTimestamp = '20200519_153722'; %'20200519_153722' (feat 3016),'20191024_122847' (feat 4548)

n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features. =17
classVar = {'strain_name','wormNum'}; 

strains2keep = {}; % Use all strains if cell left empty. {'all'} or {'divergent'} or {'controls'} or {'strain1', 'strain2'}. Cell array containing strains to keep for analysis. 
strains2drop = {}; % {'DA609','ECA252','LSJ1'} not in genoDM; Cell array containing strains to drop from analysis.
feats2keep = {}; % Use all features if left empty. {'Tierpsy_256'} or {'feat1','feat2'}. Cell array containing features to use for analysis. 
feats2drop = {}; % {'path'};

%% Load and process features table and extract features matrix
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' featExtractTimestamp '.csv'],'Delimiter',',','preserveVariableNames',true);
% filter featureTable based on specified strain and features
[featureTable, classLabels] = filterFeatureTable(featureTable,classVar,n_nonFeatVar,strains2keep,strains2drop,feats2keep,feats2drop);
% add retained labels to featureTable
for varCtr = 1:numel(classVar)
    var = classVar{varCtr};
    featureTable.(var) = classLabels.(var);
end
% extract featureMat
featureMat = table2array(featureTable(:, 1:end-numel(classVar)));
n_strains = numel(unique(featureTable.strain_name));

%% Analyze features with PCA
% pre-process feature matrix for PCA
[featureMat,~] = preprocessFeatMat(featureMat);
n_feats = size(featureMat,2);
% do pca
[pc, score, ~, ~, explained] = pca(featureMat);

%% Plot first two PCs and colour by different variables of interest

% by worm number
wormNumFig = figure; hold on
fortyLogInd = featureTable.wormNum==40;
fiveLogInd = featureTable.wormNum==5;
plot(score(fortyLogInd,1),score(fortyLogInd,2),'r.')
plot(score(fiveLogInd,1),score(fiveLogInd,2),'b.')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
legend({'40 worms','5 worms'})
title(['PCA plot with ' num2str(n_strains) ' and ' num2str(n_feats) ' features'])

% by strain
strainFig = figure; hold on
N2LogInd = strcmp(featureTable.strain_name,'N2');
CB4856LogInd = strcmp(featureTable.strain_name,'CB4856');
DA609LogInd = strcmp(featureTable.strain_name,'DA609');
NIC258LogInd = strcmp(featureTable.strain_name,'NIC258');
otherStrainsLogInd = ~N2LogInd & ~CB4856LogInd & ~DA609LogInd;
plot(score(otherStrainsLogInd,1),score(otherStrainsLogInd,2),'y.')
plot(score(N2LogInd,1),score(N2LogInd,2),'c.')
plot(score(CB4856LogInd,1),score(CB4856LogInd,2),'r.')
plot(score(DA609LogInd,1),score(DA609LogInd,2),'b.')
plot(score(NIC258LogInd,1),score(NIC258LogInd,2),'m*')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
legend({'others','N2','CB4856','DA609','NIC258'})
title(['PCA plot with ' num2str(n_strains) ' strains and ' num2str(n_feats) ' features'])

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
title(['PCA plot with ' num2str(n_strains) ' strains and ' num2str(n_feats) ' features'])

% by strain
figure;
scatter3(score(otherStrainsLogInd,1),score(otherStrainsLogInd,2),score(otherStrainsLogInd,3),'y.')
hold on
scatter3(score(N2LogInd,1),score(N2LogInd,2),score(N2LogInd,3),'c.')
scatter3(score(CB4856LogInd,1),score(CB4856LogInd,2),score(CB4856LogInd,3),'r.')
scatter3(score(DA609LogInd,1),score(DA609LogInd,2),score(DA609LogInd,3),'b.')
scatter3(score(NIC258LogInd,1),score(NIC258LogInd,2),score(NIC258LogInd,3),'m*')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
zlabel(['PC3 (' num2str(round(explained(3))) ')%'])
legend({'others','N2','CB4856','DA609','NIC258'})
title(['PCA plot with ' num2str(n_strains) ' strains and ' num2str(n_feats) ' features'])

%% Plot variance explained as a function of number of PCs
figure; hold on
% both densities PCA
plot(cumsum(explained),'-'); 
% 5 worm only PCA
fiveWormLogInd = featureTable.wormNum==5;
fiveWormFeatureMat = featureMat(fiveWormLogInd,:);
[pc, score, ~, ~, explained] = pca(fiveWormFeatureMat);
plot(cumsum(explained),'-');
% 40 worm only PCA
fortyWormLogInd = featureTable.wormNum==40;
fortyWormFeatureMat = featureMat(fortyWormLogInd,:);
[pc, score, ~, ~, explained] = pca(fortyWormFeatureMat);
plot(cumsum(explained),'-'); 
% format plot
xlim([0 500]); ylim([0 100]); xlabel('Number of PCs'); ylabel('Variance explained (%)'); legend({'All','5 worms','40 worms'})
if strcmp(feats2keep,'Tierpsy_256')
    xlim([0 80])
end