%% This script takes automatically extracted Tierpsy features from 5 vs. 40 worm experiments, 
% compare each feature at both densities, and lists features that significantly differ between the two densities.
% Significance is determined by two-sample t-test and Wilcoxon rank sum
% tests, as some features have Gaussian distribution while others do not.
% author: serenading. May 2020

% issues: 
% 1. Should also do PCA plot and colour by density to see separation by strain vs feats. 
% 2. Also when using a reduced set of features the 10% feature will be a sensitive measurement.
% 3. Multiple comparison and correction?
% 4. Correct for size/speed? Probably already taken into account by feature normalisation

clear 
close all

pThreshold = 0.01;

%% Import features matrices and load strains list

featExtractTimestamp = '20191024_122847'; %'20191024_122847' or '20181203_141111'
featureTable = readtable(['/Users/sding/Dropbox/aggScreening/results/fullFeaturesTable_' featExtractTimestamp '.csv'],'Delimiter',',','preserveVariableNames',true);
load('strainsList/all.mat','strains')
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features
n_Feats = size(featureTable,2)-n_nonFeatVar;
n_strains = numel(strains);
featNames = featureTable.Properties.VariableNames(n_nonFeatVar+1:end);
featureMat = table2array(featureTable(:,n_nonFeatVar+1:end)); % featureMat is basically featureTable without the first 17 columns of non feature-value entries.

%% Conduct t-test and fill in p-values
if ~exist('/Users/sding/Dropbox/aggScreening/results/densityDependence_sw/strainByFeat_ttest_pvalues.mat')
    % Initiate new matrix to contain significant results
    swDensityEffectP = NaN(n_strains,n_Feats); % double precision for p-values
    % Go through strain by strain
    for strainCtr = 1:numel(strains)
        strain = char(strains(strainCtr));
        fiveLogInd = strcmp(featureTable.strain_name, strain) & featureTable.wormNum==5;
        fortyLogInd = strcmp(featureTable.strain_name, strain) & featureTable.wormNum==40;
        % Go through feature by feature
        for featCtr = 1: size(featureMat,2)
            % Do two-sample t-test
            fiveFeatVal = featureMat(fiveLogInd,featCtr);
            fortyFeatVal = featureMat(fortyLogInd,featCtr);
            [~,p] = ttest2(fiveFeatVal,fortyFeatVal);
            % Save p-value
            swDensityEffectP(strainCtr,featCtr) = p;
        end
        % Display progress
        disp(['all features compared for ' num2str(strainCtr) ' out of ' num2str(n_strains) ' strains.'])
    end
    save('/Users/sding/Dropbox/aggScreening/results/densityDependence_sw/strainByFeat_ttest_pvalues.mat','swDensityEffectP','strains','featNames');
else
    load('/Users/sding/Dropbox/aggScreening/results/densityDependence_sw/strainByFeat_ttest_pvalues.mat','swDensityEffectP','strains','featNames')
end

%% Apply threshold p-value levels
display(['Applying p-value threshold of ' num2str(pThreshold)])
sigLogInd = swDensityEffectP<=pThreshold;
swDensitySigP = swDensityEffectP;
swDensitySigP(~sigLogInd)=NaN;

%% Look into what's throwing significance 
% Go through strain by strain
keyStrainCtr = 1;
clear keyStrains n_sigFeat_keyStrains
for strainCtr = 1:n_strains
    n_sigFeat = nnz(~isnan(swDensitySigP(strainCtr,:)));
    % Take note of strain if more than 10% of all features show density-dependent difference in feature values
    if n_sigFeat>0.1*n_Feats
        keyStrains{keyStrainCtr} = strains{strainCtr};
        n_sigFeat_keyStrains(keyStrainCtr) = n_sigFeat;
        keyStrainCtr = keyStrainCtr+1;
    end
end
% Display and save key strain results
disp([num2str(numel(keyStrains)) ' out of ' num2str(n_strains) ' strains show density-dependence in at least 10% of all features.'])
save(['/Users/sding/Dropbox/aggScreening/results/densityDependence_sw/keyStrains_p_' num2str(pThreshold) '.mat'],'keyStrains','n_sigFeat_keyStrains');

% Go through feature by feature
keyFeatCtr = 1;
clear keyFeats n_sigStrains_keyFeats
for featCtr = 1:n_Feats
    n_sigStrain = nnz(~isnan(swDensitySigP(:,featCtr)));
    % Take note of feature if more than 10% of all strains show density-dependent difference in feature values
    if n_sigStrain>0.1*n_strains
        keyFeats{keyFeatCtr} = featNames{featCtr};
        n_sigStrain_keyFeats(keyStrainCtr) = n_sigStrain;
        keyFeatCtr = keyFeatCtr+1;
    end
end
% Display and save key strain results
disp([num2str(numel(keyFeats)) ' out of ' num2str(n_Feats) ' features show density-dependence in at least 10% of all strains.'])
save(['/Users/sding/Dropbox/aggScreening/results/densityDependence_sw/keyFeatures_p_' num2str(pThreshold) '.mat'],'keyFeats','n_sigStrain_keyFeats');