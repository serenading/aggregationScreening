%% script performs PCA analysis with all features from all recordings (both 40 and 5 worms) in the full dataset, and plots in shared PC space
% author: @serenading May 2020

clear
close all

addpath('auxiliary/')

% TODO: rearrange phenotable to match the order of geno table.

%% Set parameters
featExtractTimestamp = '20200519_153722'; %'20200519_153722','20191024_122847' or '20181203_141111'
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features. =17
feats2drop = {'path','blob'};
strains2drop = {}; % {'DA609','ECA252','LSJ1'} are not in Ce328 geno tree.

%% Load features table and extract features matrix
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' featExtractTimestamp '.csv'],'Delimiter',',','preserveVariableNames',true);
% Drop features and strains as specified
[featureTable, ~] = dropFeats(featureTable,feats2drop);
[featureTable,~] = dropStrains(featureTable,strains2drop);
featureMat = table2array(featureTable(:, n_nonFeatVar+1:end));
featNames = featureTable.Properties.VariableNames(n_nonFeatVar+1:end);
n_feats = size(featureMat,2);
uniqueStrains = unique(featureTable.strain_name);
n_strains = numel(uniqueStrains);

%% Average within each strain
% initiate new feature matrix to hold mean values at the strain level
meanFeatureMat = NaN(n_strains,n_feats);
for strainCtr = 1:n_strains
    strain =  uniqueStrains(strainCtr);
    strainLogInd = strcmp(featureTable.strain_name,strain);
    meanFeatureMat(strainCtr,:) = mean(featureMat(strainLogInd,:),1);
end
% replace featureTable with the averaged one
clear featureTable
featureTable = array2table(meanFeatureMat,'VariableNames',featNames);
featureTable.strain_name = uniqueStrains;
featureMat = meanFeatureMat;

%% Match geno-pheno strains
% Load genotype distance matrix based on SNP variants
genoDM = readtable('/Users/sding/OneDrive - Imperial College London/aggScreening/Ce328_distance_matrix.csv','ReadVariableNames',true,'ReadRowNames',true);
assert(nnz(strcmp(genoDM.Properties.VariableNames, genoDM.Properties.RowNames')) == numel(genoDM.Properties.VariableNames),...
    'Genotype distance matrix column headings and row names do not match 100%.')
% Get the list of strains for which both genotype and phenotype data exist
genoStrains = unique(genoDM.Properties.VariableNames);
% Get strain logical index
genoStrainLogInd = false(1,numel(genoStrains));
phenoStrainLogInd = true(1,numel(uniqueStrains));
for phenoStrainCtr = 1:numel(uniqueStrains)
    strainIdx = find(ismember(genoStrains,uniqueStrains{phenoStrainCtr}));
    if isempty(strainIdx)
        phenoStrainLogInd(phenoStrainCtr) = false;
        disp([uniqueStrains{phenoStrainCtr} ' dropped from phenotype analysis.'])
    end
    genoStrainLogInd(strainIdx) = true;
end
% Trim geno
genoDM = genoDM(genoStrainLogInd,genoStrainLogInd);
% Trim pheno
featureTable = featureTable(phenoStrainLogInd,:);
featureMat = featureMat(phenoStrainLogInd,:);

%% Analyze features with PCA
% pre-process feature matrix for PCA
[featureMat,~] = preprocessFeatMat(featureMat);
% do pca
[pc, score, ~, ~, explained] = pca(featureMat);