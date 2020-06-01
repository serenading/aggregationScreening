%% This script analyses density effects in isolated single (skeletonisable) worms, 
%% using automatically extracted Tierpsy features from 5 vs. 40 worm experiments.
% Optional: When "mixed" test is selected, significance is determined by two-sample t-test if the feature is Gaussian, and by ranksum test if the feature is not.
% keyStrains are strains with higher than featThresh fraction of features altered between the two densities;
% keyFeats are feats that are altered in higher than strainThresh fraction of strains between the two densities.

% author: serenading. May 2020

clear
close all

addpath('auxiliary/')

%% Set analysis parameters
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features. =17
featExtractTimestamp = '20200519_153722'; %'20200519_153722' (feat 3016),'20200511_162714' (feat 3016 three windows) or '20191024_122847' (feat 4548) or '20181203_141111' (feat 4548)
if strcmp(featExtractTimestamp,'20200511_162714')
    featExtractWindow = '1'; %'0','1','2'
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
else
    extractStamp = featExtractTimestamp;
end

dropFeatThreshold = 0.2; % the maximum fraction of NaN values that a feature can have before being dropped
feats2drop = {'path'}; % {} if none, or cell array of strings specifying features to drop {'path','blob'}.
featThresh = 0.1; % fraction of altered features to quality a strain as keyStrain
strainThresh = 0.1; % fraction of strains to quality an altered feature as keyFeat
whichTest = 't'; % 'mixed','t'
bonCorr = true; % apply Bonferroni correction for multiple comparisons

%% Import features matrices and load strains list
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
load('strainsList/all.mat','strains')
n_feats = size(featureTable,2)-n_nonFeatVar;
n_strains = numel(strains);
featNames = featureTable.Properties.VariableNames(n_nonFeatVar+1:end);
featureMat = table2array(featureTable(:,n_nonFeatVar+1:end)); % featureMat is basically featureTable without the first 17 columns of non feature-value entries.

%% Conduct t-test and fill in p-values
testValFilename = ['/Users/sding/OneDrive - Imperial College London/aggScreening/results/densityDependence_sw/strainByFeat_' whichTest 'test_pvalues_' extractStamp '.mat'];
if ~exist(testValFilename)
    % Initiate new matrix to contain significant results
    swDensityEffectP = NaN(n_strains,n_feats); % double precision for p-values
    % Get the list of which features have too many NaN and are normally distributed
    if strcmp(featExtractTimestamp,'20200511_162714')
        load('/Users/sding/OneDrive - Imperial College London/aggScreening/results/featuresDistribution/whichFeatNormalSWTest_20200519_153722.mat','normalFeatNames','nonNormalFeatNames','NaNFeatNames');
    else
        load(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/featuresDistribution/whichFeatNormalSWTest_' extractStamp '.mat'],'normalFeatNames','nonNormalFeatNames','NaNFeatNames');
    end
    % Go through strain by strain
    for strainCtr = 1:numel(strains)
        strain = char(strains(strainCtr));
        fiveLogInd = strcmp(featureTable.strain_name, strain) & featureTable.wormNum==5;
        fortyLogInd = strcmp(featureTable.strain_name, strain) & featureTable.wormNum==40;
        % Go through feature by feature
        for featCtr = 1: size(featureMat,2)
            fiveFeatVal = featureMat(fiveLogInd,featCtr);
            fortyFeatVal = featureMat(fortyLogInd,featCtr);
            if strcmp(whichTest,'mixed')
                if ismember(featNames(featCtr),normalFeatNames)
                    [~,p] = ttest2(fiveFeatVal,fortyFeatVal);
                elseif ismember(featNames(featCtr),nonNormalFeatNames)
                    p = ranksum(fiveFeatVal,fortyFeatVal);
                end
            elseif strcmp(whichTest,'t')
                if ~ismember(featNames(featCtr),NaNFeatNames)
                    [~,p] = ttest2(fiveFeatVal,fortyFeatVal);
                end
            elseif strcmp(whichTest,'ranksum')
                if ~ismember(featNames(featCtr),NaNFeatNames)
                    p = ranksum(fiveFeatVal,fortyFeatVal);
                end
            else
                error('Please specify a valid whichTest')
            end
            % Save p-value
            swDensityEffectP(strainCtr,featCtr) = p;
        end
        % Display progress
        disp(['All features compared for ' num2str(strainCtr) ' out of ' num2str(n_strains) ' strains.'])
    end
    save(testValFilename,'swDensityEffectP','strains','featNames');
else
    load(testValFilename,'swDensityEffectP','strains','featNames')
end

%% Drop specified features from analysis
if ~isempty(feats2drop)
    [featureTable, dropLogInd] = dropFeats(featureTable,feats2drop);
    dropLogInd = dropLogInd(n_nonFeatVar+1:end);
    assert(numel(dropLogInd) == size(swDensityEffectP,2));
    swDensityEffectP = swDensityEffectP(:,~dropLogInd);
    featNames = featNames(~dropLogInd);
    n_feats = numel(featNames);
    assert(size(featureTable,2)-n_nonFeatVar == size(swDensityEffectP,2));
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
    if n_sigFeat > strainThresh * n_feats
        keyStrains{keyStrainCtr} = strains{strainCtr};
        n_sigFeat_keyStrains(keyStrainCtr) = n_sigFeat;
        keyStrainCtr = keyStrainCtr+1;
    end
end
% Display and save key strain results
if keyStrainCtr>1
    disp([num2str(numel(keyStrains)) ' out of ' num2str(n_strains) ' strains show density-dependence in at least ' num2str(featThresh*100) '% of all features.'])
else
    disp(['No strain show density-dependence in at least ' num2str(featThresh*100) '% of all features.'])
end
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
    if n_sigStrain > strainThresh * n_strains
        keyFeats{keyFeatCtr} = featNames{featCtr};
        n_sigStrain_keyFeats(keyStrainCtr) = n_sigStrain;
        keyFeatCtr = keyFeatCtr+1;
    end
end
% Display and save key strain results
if keyFeatCtr>1
    disp([num2str(numel(keyFeats)) ' out of ' num2str(n_feats) ' features show density-dependence in at least ' num2str(strainThresh*100) '% of all strains.'])
else
    disp(['No feature shows density-dependence in at least ' num2str(strainThresh*100) '% of all strains.'])
end
figure; histogram(strainsRatio)
xlabel('proportion of strains showing density dependence')
ylabel('number of features')

%% Results: keyFeats that show density-dependence in a large number of strains include: 
%     {'motion_mode_forward_duration_50th'                             }
%     {'motion_mode_forward_fraction'                                  }
%     {'motion_mode_paused_fraction'                                   }
%     {'food_region_inside_duration_50th'                              }
%     {'food_region_edge_duration_50th'                                }
%     {'turn_intra_duration_50th'                                      }
%     {'turn_intra_frequency'                                          }
%     {'relative_to_head_base_radial_velocity_head_tip_w_backward_50th'}
%     {'d_speed_w_backward_50th'                                       }
%     {'orientation_food_edge_w_forward_IQR'                           }

%% Plot key feature distributions
keyFeats'
strains2plot = {'N2','DA609','CB4856'};

figure;
% go through each keyFeat
n_keyFeats = numel(keyFeats);
for keyFeatCtr = 1:n_keyFeats
    keyFeat = keyFeats{keyFeatCtr};
    keyFeatColIdx = find(strcmp(featureTable.Properties.VariableNames,keyFeat));
    % initiate
    groupCtr = 1;
    % add violin plot info
    % 40, all strains
    fortyRowLogInd = featureTable.wormNum==40;
    violinplotVals{:,groupCtr} = table2array(featureTable(fortyRowLogInd,keyFeatColIdx));
    groupCtr = groupCtr+1;
    % 5, all strains
    fiveRowLogInd = featureTable.wormNum==5;
    violinplotVals{:,groupCtr} = table2array(featureTable(fiveRowLogInd,keyFeatColIdx));
    groupCtr = groupCtr+1;
    % strain by strain
    for strainCtr = 1:numel(strains2plot)
        strain = strains2plot{strainCtr};
        % 40
        fortyStrainRowLogInd = strcmp(featureTable.strain_name,strain) & featureTable.wormNum==40;
        violinplotVals{:,groupCtr} = table2array(featureTable(fortyStrainRowLogInd,keyFeatColIdx));
        groupCtr = groupCtr+1;
        % 5
        fiveStrainRowLogInd = strcmp(featureTable.strain_name,strain) & featureTable.wormNum==5;
        violinplotVals{:,groupCtr} = table2array(featureTable(fiveStrainRowLogInd,keyFeatColIdx));
        groupCtr = groupCtr+1;
    end
    % plot boxplot/violin plot
    if n_keyFeats <10
        subplot(2,ceil(numel(keyFeats)/2),keyFeatCtr);
    else
        subplot(3,ceil(numel(keyFeats)/3),keyFeatCtr);
    end
    violin(violinplotVals,'facecolor',[1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1]);
    xticks([2,4,6,8])
    xticklabels(['all',strains2plot])
    legend off
    title(keyFeats{keyFeatCtr},'Interpreter','None')
end
legend({'40','Mean','Median','5'}) % add legend to final violin plot