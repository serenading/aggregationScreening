%% Script compares distance matrices  based on genotype (SNPs) and phenotype (Tierpsy features) across strains. 
% Useful option: makeMappingFile: outputs a .tsv file in the correct format for cegwas2-nf mapping

% author: @serenading June 2020

clear
close all

addpath('auxiliary/')

%% Set parameters
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features. =17
classVar = 'strain_name'; 
wormDensity = 5; % 5 or 40. wormDensities are not to be combined because feature values are averaged across different replicates for analysis. 
dropGenoOutlierStrains = true;
makeMappingFile = false;

% set which feature extraction timestamp to use
featExtractTimestamp = '20200519_153722'; %'20200519_153722' (feat 3016),'20200511_162714' (feat 3016 three windows) or '20191024_122847' (feat 4548)
if strcmp(featExtractTimestamp,'20200511_162714')
    featExtractWindow = '1'; %'0','1','2'
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
else
    extractStamp = featExtractTimestamp;
end

% select which features and features to drop or keep
strains2keep = {}; % Use all strains if cell left empty. {'all'} or {'divergent'} or {'controls'} or {'strain1', 'strain2'}. Cell array containing strains to keep for analysis. 
strains2drop = {'DA609','ECA252','LSJ1'}; % {'DA609','ECA252','LSJ1'} not in genoDM;  Cell array containing strains to drop from analysis.
if dropGenoOutlierStrains
    strains2drop = [strains2drop {'ECA36','ECA363','ECA396','XZ1514','XZ1516'}]; % {'ECA36','ECA363','ECA396','XZ1514','XZ1516'} are geno outlier strains.
end
feats2keep = {}; % Use all features if cell left empty. {'Tierpsy_256'} or {'feat1','feat2'}. Cell array containing features to use for analysis. 
feats2drop = {}; % {'path'};

%% Load and process featureTable
% load
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);

% use only data from specified worm density to compute phenotypic distances
wormNumLogInd = featureTable.wormNum == wormDensity;
featureTable = featureTable(wormNumLogInd,:);

% filter featureTable based on specified strain and features
[featureTable, classLabels] = filterFeatureTable(featureTable,classVar,n_nonFeatVar,strains2keep,strains2drop,feats2keep,feats2drop);

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
featureTable.(classVar) = classLabels;

%% Average within each strain
uniqueStrains = unique(classLabels);
n_strains = numel(uniqueStrains);
n_feats = size(featureMat,2);

% initiate new feature matrix to hold mean values at the strain level
meanFeatureMat = NaN(n_strains,n_feats);
for strainCtr = 1:n_strains
    strain =  uniqueStrains(strainCtr);
    strainLogInd = strcmp(featureTable.(classVar),strain);
    meanFeatureMat(strainCtr,:) = mean(featureMat(strainLogInd,:),1);
end
% replace featureTable and feature matrix with the averaged onewritetable(mappingTable,filename,'Delimiter','\t');
clear featureTable
featureTable = array2table(meanFeatureMat,'VariableNames',featNames);
featureTable.(classVar) = uniqueStrains;
featureMat = meanFeatureMat;

%% generate .tsv for GWAS mapping
if makeMappingFile
    % first column is strain names
    strain = uniqueStrains;
    mappingTable = table(strain);
    % the other columns are features
    mappingTable2 = array2table(featureMat);
    mappingTable2.Properties.VariableNames = featNames;
    % concatenate the tables
    mappingTable = horzcat(mappingTable,mappingTable2);
    % exclude DA609 from mapping
    nonDA609LogInd = ~strcmp(mappingTable.strain, 'DA609');
    mappingTable = mappingTable(nonDA609LogInd,:);
    % write to tsv
    filename = ['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/fiveWormMappingFile_' extractStamp char(feats2keep) '.txt'];
    writetable(mappingTable,filename,'Delimiter','\t');
    filename2 = strrep(filename,'.txt','.tsv');
    movefile(filename,filename2) % this is optionally really, as cegwas2-nf runs just fine with the tab-delimited .txt file.
end

%% Load genotype distance matrix based on SNP variants
genoDM = readtable('/Users/sding/OneDrive - Imperial College London/aggScreening/source/Ce328_distance_matrix.csv','ReadVariableNames',true,'ReadRowNames',true);
assert(nnz(strcmp(genoDM.Properties.VariableNames, genoDM.Properties.RowNames')) == numel(genoDM.Properties.VariableNames),...
    'Genotype distance matrix column headings and row names do not match 100% before filtering for phenotyped strains.')

%% Match the number and order of geno-pheno strains
% the next line also sorts strain order in genoDM according to that of the pheno strain list, which is already in alphabetical order. 
genoDM = genoDM(uniqueStrains,uniqueStrains);
assert(nnz(strcmp(genoDM.Properties.RowNames,genoDM.Properties.VariableNames')) == numel(genoDM.Properties.VariableNames),...
    'Genotype distance matrix column headings and row names do not match 100% after filtering for phenotyped strains.')
assert(nnz(strcmp(genoDM.Properties.RowNames,featureTable.(classVar))) == numel(uniqueStrains),'Strain order is not the in GenoDM as in PhenoFeatureTable.')

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
genoDM_v2(genoDM_v>0.1)=median(genoDM_v); % outlier strains: {'ECA36','ECA363','ECA396','XZ1514','XZ1516'}
genoDM_s2 = squareform(genoDM_v2);

%% Compare geno vs. pheno distance matrices
if dropGenoOutlierStrains
    nPlotCols = 2;
else 
    nPlotCols = 3;
end

% plot histograms of distance matrix values
figure;
subplot(2,nPlotCols,1);
histogram(phenoDM_v,'Normalization','probability','BinWidth',0.5); xlim([0 60]); ylim([0 0.06]); title('Phenotype DM vector'); xlabel('Distance'); ylabel('P');
subplot(2,nPlotCols,2);
histogram(genoDM_v*scaleFactor,'Normalization','probability','BinWidth',0.5); title('Genotype DM vector'); xlabel('Scaled distance'); ylabel('P');
if dropGenoOutlierStrains
    xlim([0 60]); ylim([0 0.06]);
else
    subplot(2,nPlotCols,3);
    histogram(genoDM_v2*scaleFactor,'Normalization','probability','BinWidth',0.5); xlim([0 60]);ylim([0 0.06]); title('Genotype DM vector (no outliers)'); xlabel('Scaled distance'); ylabel('P');
end

% visualise matrices
subplot(2,nPlotCols,nPlotCols+1);
imagesc(phenoDM_s); colorbar; caxis([0 50]); title('Phenotype Distance Matrix'); xlabel('Strains'); ylabel('Strains');
subplot(2,nPlotCols,nPlotCols+2);
imagesc(genoDM_s*scaleFactor); colorbar; caxis([0 50]); title('Genotype Distance Matrix'); xlabel('Strains'); ylabel('Strains');
if ~dropGenoOutlierStrains
    subplot(2,nPlotCols,nPlotCols+3);
    imagesc(genoDM_s2*scaleFactor); colorbar; caxis([0 50]); title('Genotype Distance Matrix (no outliers)'); xlabel('Strains'); ylabel('Strains');
end

% scatterplot of all distances
figure; 
if ~dropGenoOutlierStrains
    subplot(1,2,2);
    scatter(genoDM_v2*scaleFactor,phenoDM_v,'.'); xlabel('scaled genoDist (no outliers)'); ylabel('phenoDist'); xlim([0 60]); ylim([0 60]);
    subplot(1,2,1);
end
scatter(genoDM_v*scaleFactor,phenoDM_v,'.'); xlabel('scaled genoDist'); ylabel('phenoDist'); xlim([0 60]); ylim([0 60]);

%%  Calculate correlation (Options: linear (Pearson) or monotonic (Spearman))
% Mantel test of matrice dissimilarity
[coeff, pval] = bramila_mantel(phenoDM_s,genoDM_s*scaleFactor,5000,'pearson');
disp(['Mantel test of matrix dissimilarity: coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.'])
if ~dropGenoOutlierStrains
    [coeff, pval] = bramila_mantel(phenoDM_s,genoDM_s2*scaleFactor,5000,'pearson');
    disp(['Mantel test of matrix dissimilarity (after removing outliers from genoDM): coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.'])
end

% Pearson and Spearman correlations of vectors
[coeff,pval] = corr(phenoDM_v', genoDM_v');
disp(['Pearson correlation of distance vectors: coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.']);
[coeff,pval] = corr(phenoDM_v', genoDM_v','Type','Spearman');
disp(['Spearman correlation of distance vectors: coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.']);
if ~dropGenoOutlierStrains
    [coeff,pval] = corr(phenoDM_v', genoDM_v2');
    disp(['Pearson correlation of distance vectors (after removing outliers from genoDist): coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.']);
    [coeff,pval] = corr(phenoDM_v', genoDM_v2','Type','Spearman');
    disp(['Spearman correlation of distance vectors (after removing outliers from genoDist): coefficient is ' num2str(coeff) ' and p-value is ' num2str(pval) '.']);
end

%% Agglomerative hierarchical clustering
phenoLink = linkage(phenoDM_v,'average');
genoLink = linkage(genoDM_v,'average');
genoLink2 = linkage(genoDM_v2,'average');

% dendrograms from hierarchical clustering
figure; 
subplot(1,nPlotCols,1);dendrogram(phenoLink,'Orientation','right','Labels',uniqueStrains); title('PhenoDM linkage'); xlabel('Distance')
subplot(1,nPlotCols,2);dendrogram(genoLink,'Orientation','right','Labels',uniqueStrains); title('GenoDM linkage'); xlabel('Distance')
if ~dropGenoOutlierStrains
    subplot(1,nPlotCols,3);dendrogram(genoLink2,'Orientation','right','Labels',uniqueStrains); title('GenoDM linkage (no outliers)'); xlabel('Distance')
end