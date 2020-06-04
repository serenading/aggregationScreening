%% script performs PCA analysis with all features from all recordings (both 40 and 5 worms) in the full dataset, and plots in shared PC space
% author: @serenading May 2020

clear
close all

addpath('auxiliary/')

%% Set parameters
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features. =17
classVar = 'strain_name'; 

% set which feature extraction timestamp to use
featExtractTimestamp = '20200519_153722'; %'20200519_153722' (feat 3016),'20200511_162714' (feat 3016 three windows) or '20191024_122847' (feat 4548)
if strcmp(featExtractTimestamp,'20200511_162714')
    featExtractWindow = '1'; %'0','1','2'
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
else
    extractStamp = featExtractTimestamp;
end

% select which features and features to drop or keep
strains2keep = {'all'}; % Use all strains if cell left empty. {'all'} or {'divergent'} or {'controls'} or {'strain1', 'strain2'}. Cell array containing strains to keep for analysis. 
strains2drop = {'DA609','ECA252','LSJ1'}; % {'DA609','ECA252','LSJ1'} not in genoDM; Cell array containing strains to drop from analysis.
feats2keep = {'Tierpsy_256'}; % Use all features if left empty. {'Tierpsy_256'} or {'feat1','feat2'}. Cell array containing features to use for analysis. 
feats2drop = {}; % {'path'};

%% Load and process featureTable
% load
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);

% use only 5 worm data to compute phenotypic distances
fivewormLogInd = featureTable.wormNum == 5;
featureTable = featureTable(fivewormLogInd,:);

% filter featureTable based on specified strain and features
[featureTable, strainNames] = filterFeatureTable(featureTable,'strain_name',n_nonFeatVar,strains2keep,strains2drop,feats2keep,feats2drop);

%% Pre-process features matrix and turn back into table 
% split table into matrix and featNames
featureMat = table2array(featureTable);
featNames = featureTable.Properties.VariableNames;
% preprocess feature matrix: drop zero standard deviation, NaN's, z-normalise, etc. 
[featureMat,dropLogInd] = preprocessFeatMat(featureMat);
featNames = featNames(~dropLogInd);
% put the table back together
featureTable = array2table(featureMat,'VariableNames',featNames);
% Add classLabels to featureTable for classification
featureTable.strain_name = strainNames;

%% Average within each strain
uniqueStrains = unique(strainNames);
n_strains = numel(uniqueStrains);
n_feats = size(featureMat,2);

% initiate new feature matrix to hold mean values at the strain level
meanFeatureMat = NaN(n_strains,n_feats);
for strainCtr = 1:n_strains
    strain =  uniqueStrains(strainCtr);
    strainLogInd = strcmp(featureTable.strain_name,strain);
    meanFeatureMat(strainCtr,:) = mean(featureMat(strainLogInd,:),1);
end
% replace featureTable and feature matrix with the averaged one
clear featureTable
featureTable = array2table(meanFeatureMat,'VariableNames',featNames);
featureTable.strain_name = uniqueStrains;
featureMat = meanFeatureMat;

%% Load genotype distance matrix based on SNP variants
genoDM = readtable('/Users/sding/OneDrive - Imperial College London/aggScreening/source/Ce328_distance_matrix.csv','ReadVariableNames',true,'ReadRowNames',true);
assert(nnz(strcmp(genoDM.Properties.VariableNames, genoDM.Properties.RowNames')) == numel(genoDM.Properties.VariableNames),...
    'Genotype distance matrix column headings and row names do not match 100% before filtering for phenotyped strains.')

%% Match the number and order of geno-pheno strains
% the next line also sorts strain order in genoDM according to that of the pheno strain list, which is already in alphabetical order. 
genoDM = genoDM(uniqueStrains,uniqueStrains);
assert(nnz(strcmp(genoDM.Properties.RowNames,genoDM.Properties.VariableNames')) == numel(genoDM.Properties.VariableNames),...
    'Genotype distance matrix column headings and row names do not match 100% after filtering for phenotyped strains.')
assert(nnz(strcmp(genoDM.Properties.RowNames,featureTable.strain_name)) == numel(uniqueStrains),'Strain order is not the in GenoDM as in PhenoFeatureTable.')

%% Get distance matrices
% compute Phenotype distance matrix as a vector
phenoDM_v = pdist(featureMat); % pdist uses euclidean distance by default.
% convert genoDM into pdist output vector form
genoDM_v = squareform(table2array(genoDM),'tovector');
% get the scale factor between the two distance matrices
scaleFactor = median(phenoDM_v)/median(genoDM_v);
% get both matrices in squareform
phenoDM_s = squareform(phenoDM_v);
genoDM_s = squareform(genoDM_v);
% drop outlier genoDM values
genoDM_v2 = genoDM_v;
genoDM_v2(genoDM_v>0.1)=median(genoDM_v);
genoDM_s2 = squareform(genoDM_v2);

%% Compare geno vs. pheno distance matrices
% plot histograms of distance matrix values
figure;
subplot(2,3,1); histogram(phenoDM_v,'Normalization','probability'); xlim([0 60]); title('Phenotype DM vector'); xlabel('Distance'); ylabel('P');
subplot(2,3,2); histogram(genoDM_v*scaleFactor,'Normalization','probability'); title('Genotype DM vector'); xlabel('Scaled distance'); ylabel('P');
subplot(2,3,3); histogram(genoDM_v2*scaleFactor,'Normalization','probability'); xlim([0 60]); title('Genotype DM vector (no outliers)'); xlabel('Scaled distance'); ylabel('P');

% visualise matrices
subplot(2,3,4); imagesc(phenoDM_s); colorbar; caxis([0 50]); title('Phenotype Distance Matrix');
subplot(2,3,5); imagesc(genoDM_s*scaleFactor); colorbar; caxis([0 50]); title('Genotype Distance Matrix');
subplot(2,3,6); imagesc(genoDM_s2*scaleFactor); colorbar; caxis([0 50]); title('Genotype Distance Matrix (no outliers)');

% scatterplot of all distances
figure; 
subplot(1,2,1); scatter(genoDM_v*scaleFactor,phenoDM_v,'.'); xlabel('scaled genoDist'); ylabel('phenoDist');
subplot(1,2,2); scatter(genoDM_v2*scaleFactor,phenoDM_v,'.'); xlabel('scaled genoDist (no outliers)'); ylabel('phenoDist')

%%  Calculate correlation (Options: linear (Pearson) or monotonic (Spearman))
% Mantel test of matrice dissimilarity
[coeff, pval] = bramila_mantel(phenoDM_s,genoDM_s*scaleFactor,5000,'pearson');
disp(['Mantel test of matrix dissimilarity: p-value for is ' num2str(pval) ' and coefficient is ' num2str(coeff) '.'])
[coeff, pval] = bramila_mantel(phenoDM_s,genoDM_s2*scaleFactor,5000,'pearson');
disp(['Mantel test of matrix dissimilarity: p-value for is ' num2str(pval) ' (after removing outliers from genoDM)  and coefficient is ' num2str(coeff) '.'])

% Pearson and Spearman correlations of vectors
[coeff,pval] = corr(phenoDM_v', genoDM_v');
disp(['Pearson correlation of distance vectors: coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.']);
[coeff,pval] = corr(phenoDM_v', genoDM_v','Type','Spearman');
disp(['Spearman correlation of distance vectors: coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.']);
[coeff,pval] = corr(phenoDM_v', genoDM_v2');
disp(['Pearson correlation of distance vectors (after removing outliers from genoDist): coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.']);
[coeff,pval] = corr(phenoDM_v', genoDM_v2','Type','Spearman');
disp(['Spearman correlation of distance vectors (after removing outliers from genoDist): coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.']);

%% Agglomerative hierarchical clustering
phenoLink = linkage(phenoDM_v,'average');
genoLink = linkage(genoDM_v,'average');
genoLink2 = linkage(genoDM_v2,'average');

% dendrograms from hierarchical clustering (need to understand what the leaf node numbers mean)
figure; 
subplot(1,3,1);dendrogram(phenoLink); title('PhenoDM linkage');
subplot(1,3,2);dendrogram(genoLink); title('GenoDM linkage');
subplot(1,3,3);dendrogram(genoLink2); title('GenoDM linkage (no outliers');