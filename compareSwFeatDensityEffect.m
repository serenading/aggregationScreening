%% This script analyses density effects in isolated (skeletonisable) worms, 
%% using automatically extracted Tierpsy features from 5 vs. 40 worm experiments.

% Significance is determined by two-sample t-test.
% (Tried using Wilcoxon rank sum test for features without normal distribution, but mixed results do not make much sense).

% author: serenading. May 2020

% issues:
% 1. Should also do PCA plot and colour by density to see separation by strain vs feats.

clear
close all

%% Set analysis parameters
bonCorr = true; % apply Bonferroni correction for multiple comparisons
featExtractTimestamp = '20191024_122847'; %'20191024_122847' or '20181203_141111'
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features. 17

%% Import features matrices and load strains list
featureTable = readtable(['/Users/sding/Dropbox/aggScreening/results/fullFeaturesTable_' featExtractTimestamp '.csv'],'Delimiter',',','preserveVariableNames',true);
load('strainsList/all.mat','strains')
n_feats = size(featureTable,2)-n_nonFeatVar;
n_strains = numel(strains);
featNames = featureTable.Properties.VariableNames(n_nonFeatVar+1:end);
featureMat = table2array(featureTable(:,n_nonFeatVar+1:end)); % featureMat is basically featureTable without the first 17 columns of non feature-value entries.

%% Conduct t-test and fill in p-values
if ~exist('/Users/sding/Dropbox/aggScreening/results/densityDependence_sw/strainByFeat_ttest_pvalues.mat')
    % Initiate new matrix to contain significant results
    swDensityEffectP = NaN(n_strains,n_feats); % double precision for p-values
    % Get the list of which strains are normally distributed
    load('/Users/sding/Dropbox/aggScreening/results/featuresDistribution/whichFeatNormal.mat','normalFeatNames','nonNormalFeatNames');
    % Go through strain by strain
    for strainCtr = 1:numel(strains)
        strain = char(strains(strainCtr));
        fiveLogInd = strcmp(featureTable.strain_name, strain) & featureTable.wormNum==5;
        fortyLogInd = strcmp(featureTable.strain_name, strain) & featureTable.wormNum==40;
        % Go through feature by feature
        for featCtr = 1: size(featureMat,2)
            fiveFeatVal = featureMat(fiveLogInd,featCtr);
            fortyFeatVal = featureMat(fortyLogInd,featCtr);
            [~,p] = ttest(fiveFeatVal,fortyFeatVal);
            % Save p-value
            swDensityEffectP(strainCtr,featCtr) = p;
        end
        % Display progress
        disp(['All features compared for ' num2str(strainCtr) ' out of ' num2str(n_strains) ' strains.'])
    end
    save('/Users/sding/Dropbox/aggScreening/results/densityDependence_sw/strainByFeat_ttest_pvalues.mat','swDensityEffectP','strains','featNames');
else
    load('/Users/sding/Dropbox/aggScreening/results/densityDependence_sw/strainByFeat_ttest_pvalues.mat','swDensityEffectP','strains','featNames')
end

%% Look into what's throwing significance...
%% At the strain level
% Apply threshold p-value levels
if bonCorr
    pThreshold = 0.05/n_feats;
else
    pThreshold = 0.05;
end
sigLogInd = swDensityEffectP<=pThreshold;
swDensitySigP = swDensityEffectP;
swDensitySigP(~sigLogInd)=NaN;
% Go through strain by strain
featsRatio = NaN(1,n_strains);
keyStrainCtr = 1;
clear keyStrains n_sigFeat_keyStrains
for strainCtr = 1:n_strains
    n_sigFeat = nnz(~isnan(swDensitySigP(strainCtr,:)));
    featsRatio(strainCtr) = n_sigFeat/n_feats;
    % Take note of strain if more than 10% of all features show density-dependent difference in feature values
    if n_sigFeat>0.1*n_feats
        keyStrains{keyStrainCtr} = strains{strainCtr};
        n_sigFeat_keyStrains(keyStrainCtr) = n_sigFeat;
        keyStrainCtr = keyStrainCtr+1;
    end
end
% Display and save key strain results
if keyStrainCtr>1
    disp([num2str(numel(keyStrains)) ' out of ' num2str(n_strains) ' strains show density-dependence in at least 10% of all features.'])
else
    disp('No strain show density-dependence in at least 10% of all features.')
end
%save(['/Users/sding/Dropbox/aggScreening/results/densityDependence_sw/keyStrains_p_' num2str(pThreshold) '.mat'],'keyStrains','n_sigFeat_keyStrains');
figure; histogram(featsRatio)
xlabel('proportion of features showing density dependence')
ylabel('number of strains')

%% At the feature level
% Apply threshold p-value levels
if bonCorr
    pThreshold = 0.05/n_strains;
else
    pThreshold = 0.05;
end
sigLogInd = swDensityEffectP<=pThreshold;
swDensitySigP = swDensityEffectP;
swDensitySigP(~sigLogInd)=NaN;
% Go through feature by feature
strainsRatio = NaN(1,n_feats);
keyFeatCtr = 1;
clear keyFeats n_sigStrains_keyFeats
for featCtr = 1:n_feats
    n_sigStrain = nnz(~isnan(swDensitySigP(:,featCtr)));
    strainsRatio(featCtr) = n_sigStrain/n_strains;
    % Take note of feature if more than 10% of all strains show density-dependent difference in feature values
    if n_sigStrain>0.1*n_strains
        keyFeats{keyFeatCtr} = featNames{featCtr};
        n_sigStrain_keyFeats(keyStrainCtr) = n_sigStrain;
        keyFeatCtr = keyFeatCtr+1;
    end
end
% Display and save key strain results
if keyFeatCtr>1
    disp([num2str(numel(keyFeats)) ' out of ' num2str(n_feats) ' features show density-dependence in at least 10% of all strains.'])
else
    disp('No feature shows density-dependence in at least 10% of all strains.')
end
figure; histogram(strainsRatio)
xlabel('proportion of strains showing density dependence')
ylabel('number of features')