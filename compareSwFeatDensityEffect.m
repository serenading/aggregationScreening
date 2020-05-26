%% This script analyses density effects in isolated single (skeletonisable) worms, 
%% using automatically extracted Tierpsy features from 5 vs. 40 worm experiments.
% Optional: Significance is determined by two-sample t-test if the feature is Gaussian, and by ranksum test if the feature is not.
% keyStrains are strains with higher than featThresh fraction of features altered between the two densities;
% keyFeats are feats that are altered in higher than strainThresh fraction of strains between the two densities.

% author: serenading. May 2020

%% Things to consider:
% 1. Classify strains based on 40 worms or 5 worms.
% 2. plot a few skeletons
% 3. PCA on features, colour by worm density

clear
close all

%% Set analysis parameters
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features. =17
featExtractTimestamp = '20200519_153722'; %'20200519_153722' (feat 3006),'20200511_162714' (feat 3006 three windows) or '20191024_122847' (feat 4548) or '20181203_141111' (feat 4548)
if strcmp(featExtractTimestamp,'20200511_162714')
    featExtractWindow = '2'; %'0','1','2'
    extractStamp = [featExtractTimestamp '_window_' featExtractWindow];
else
    extractStamp = featExtractTimestamp;
end

dropPathBlobFeatures = true;
featThresh = 0.1; % fraction of altered features to quality a strain as keyStrain
strainThresh = 0.1; % fraction of strains to quality an altered feature as keyFeat
whichTest = 't'; % 'mixed','t','ranksum';
bonCorr = true; % apply Bonferroni correction for multiple comparisons

%% Import features matrices and load strains list
featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
if dropPathBlobFeatures
    featureTable = dropPathBlobFeats(featureTable);
end
load('strainsList/all.mat','strains')
n_feats = size(featureTable,2)-n_nonFeatVar;
n_strains = numel(strains);
featNames = featureTable.Properties.VariableNames(n_nonFeatVar+1:end);
featureMat = table2array(featureTable(:,n_nonFeatVar+1:end)); % featureMat is basically featureTable without the first 17 columns of non feature-value entries.

%% Conduct t-test and fill in p-values
testValFilename = ['/Users/sding/OneDrive - Imperial College London/aggScreening/results/densityDependence_sw/strainByFeat_' whichTest 'test_pvalues' extractStamp '.mat'];
if dropPathBlobFeatures
    testValFilename = strrep(testValFilename,'.mat','_noPathBlob.mat');
end
if ~exist(testValFilename)
    % Initiate new matrix to contain significant results
    swDensityEffectP = NaN(n_strains,n_feats); % double precision for p-values
    if strcmp(whichTest,'mixed')
        % Get the list of which strains are normally distributed
        % this list does not work for feat 3006.
        load('/Users/sding/OneDrive - Imperial College London/aggScreening/results/featuresDistribution/whichFeatNormalSWTest.mat','normalFeatNames','nonNormalFeatNames');
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
                else
                    error(['The feature ' featNames{featCtr} ' is not pre-classified as having a normal distribution or not.'])
                end
            elseif strcmp(whichTest,'t')
                [~,p] = ttest2(fiveFeatVal,fortyFeatVal);
            elseif strcmp(whichTest,'ranksum')
                 p = ranksum(fiveFeatVal,fortyFeatVal);
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
%% other path-related and blob features are likely just sensitive to tracking at different densities. keyFeats total = 77. no keyStrain

%% Extract likely keyFeats by ignoring the ones that are likely feature calculation or tracking artifacts
keyFeatsLogInd = ~cellfun(@(x) contains(x,'path'), keyFeats) & ~cellfun(@(x) contains(x,'blob'), keyFeats);
keyFeats2check = keyFeats(keyFeatsLogInd);

%% Display mean values for the key features at each density
figure;
strains2check = {'N2','DA609','CB4856'};
% go through each keyFeat
for keyFeatCtr = 1:numel(keyFeats2check)
    keyFeat2check = keyFeats2check{keyFeatCtr};
    keyFeatColIdx = find(strcmp(featureTable.Properties.VariableNames,keyFeat2check));
    % initiate
    boxplotVals = [];
    densityLabels = [];
    strainLabels = [];
    groupCtr = 1;
    % add boxplot info
    fortyRowLogInd = featureTable.wormNum==40;
    fiveRowLogInd = featureTable.wormNum==5;
    boxplotVals = vertcat(boxplotVals,table2array(featureTable(fortyRowLogInd,keyFeatColIdx)));
    boxplotVals = vertcat(boxplotVals,table2array(featureTable(fiveRowLogInd,keyFeatColIdx)));
    densityLabels = vertcat(densityLabels, 40*ones(nnz(fortyRowLogInd),1));
    densityLabels = vertcat(densityLabels, 5*ones(nnz(fiveRowLogInd),1));
    strainLabels = vertcat(strainLabels, groupCtr*ones(nnz(fortyRowLogInd),1));
    groupCtr = groupCtr+1;
    strainLabels = vertcat(strainLabels, groupCtr*ones(nnz(fiveRowLogInd),1));
    groupCtr = groupCtr+1;
%     fortyFeatMean = mean(table2array(featureTable(fortyRowLogInd,keyFeatColIdx)));
%     fiveFeatMean = mean(table2array(featureTable(fiveRowLogInd,keyFeatColIdx)));
%     disp([sprintf('\n') 'In all strains, ' keyFeat2check ' is ' num2str(fortyFeatMean) ' at 40 worm density and ' num2str(fiveFeatMean) ' at 5 worm density.'])
    for strainCtr = 1:numel(strains2check)
        strain = strains2check{strainCtr};
        fortyStrainRowLogInd = strcmp(featureTable.strain_name,strain) & featureTable.wormNum==40;
        fiveStrainRowLogInd = strcmp(featureTable.strain_name,strain) & featureTable.wormNum==5;
%         fortyStrainFeatMean = mean(table2array(featureTable(fortyStrainRowLogInd,keyFeatColIdx)));
%         fiveStrainFeatMean = mean(table2array(featureTable(fiveStrainRowLogInd,keyFeatColIdx)));
%         disp(['In ' strain ', ' keyFeat2check ' is ' num2str(fortyStrainFeatMean) ' at 40 worm density and ' num2str(fiveStrainFeatMean) ' at 5 worm density.'])
        boxplotVals = vertcat(boxplotVals,table2array(featureTable(fortyStrainRowLogInd,keyFeatColIdx)));
        boxplotVals = vertcat(boxplotVals,table2array(featureTable(fiveStrainRowLogInd,keyFeatColIdx)));
        densityLabels = vertcat(densityLabels, 40*ones(nnz(fortyStrainRowLogInd),1));
        densityLabels = vertcat(densityLabels, 5*ones(nnz(fiveStrainRowLogInd),1));
        strainLabels = vertcat(strainLabels, groupCtr*ones(nnz(fortyStrainRowLogInd),1));
        groupCtr = groupCtr+1;
        strainLabels = vertcat(strainLabels, groupCtr*ones(nnz(fiveStrainRowLogInd),1));
        groupCtr = groupCtr+1;
    end
    % plot boxplot
    subplot(2,4,keyFeatCtr);
    boxplot(boxplotVals,strainLabels,'Labels',{'all_40','all_5','N2_40','N2_5','DA609_40','DA609_5','CB4856_40','CB4856_5'},...
        'PlotStyle','compact','BoxStyle','filled','ColorGroup',densityLabels);
    title(keyFeats2check{keyFeatCtr},'Interpreter','None')
end