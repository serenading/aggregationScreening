clear
close all

%% Script reads in featuresTable and plots stacked bar graphs for blob category and food region statistics to show comprising fractions of tracking indices.
% Author: @serenading. July 2020

%% Specify which dataset
wormNum = 5;

%% load
if wormNum == 40
    featureTablefilename = '/Users/sding/OneDrive - Imperial College London/aggScreening/results/fortyWorm/fortyWormFeaturesTable_20200519_153722_20200708_ft2891.csv';
elseif wormNum == 5
    featureTablefilename = '/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/fiveWormFeaturesTable_20200519_153722_20200708_ft3024.csv';
end
featureTable = readtable(featureTablefilename,'Delimiter',',','preserveVariableNames',true);
uniqueStrains = unique(featureTable.strain_name);

%% pre-allocate
wormFractionsMean = NaN(numel(uniqueStrains),4);
wormFractionsStd = NaN(numel(uniqueStrains),4);
foodRegionFractionsMean = NaN(numel(uniqueStrains),3);
foodRegionFractionsStd = NaN(numel(uniqueStrains),3);

%% go through each strain to get worm fractions and food region fractions

% each strain
for strainCtr = 1:numel(uniqueStrains)
    strain = uniqueStrains(strainCtr);
    strainLogInd = strcmp(featureTable.strain_name,strain);
    
    % worm category
    wormFractions = [featureTable.swFraction(strainLogInd),featureTable.pausedMwFraction(strainLogInd),featureTable.clusterFraction(strainLogInd),featureTable.tempBlobFraction(strainLogInd)];
    wormFractionsMean(strainCtr,:) = mean(wormFractions,1);
    wormFractionsSE(strainCtr,:) = std(wormFractions)/sqrt(nnz(strainLogInd));
    
    % food region
    foodRegionFractions = [featureTable.food_region_inside_fraction(strainLogInd),featureTable.food_region_edge_fraction(strainLogInd),featureTable.food_region_outside_fraction(strainLogInd)];
    foodRegionFractionsMean(strainCtr,:) = mean(foodRegionFractions,1);
    foodRegionFractionsSE(strainCtr,:) = std(foodRegionFractions)/sqrt(nnz(strainLogInd));
end

%% stacked bar plots to show proportions

% worm category
figure;bar(wormFractionsMean,'stacked')
% hold on;errorbar(cumsum(wormFractionsMean')',wormFractionsSE,'.k');
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'single worm','paused multiworm','cluster','temp'},'FontSize',15,'Location','northeastoutside')
title('tracking indices by worm category','FontSize',15)

% food region
figure;bar(foodRegionFractionsMean,'stacked')
% hold on;errorbar(cumsum(foodRegionFractionsMean')',foodRegionFractionsSE,'.k');
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylim([0 1])
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'onFood','foodEdge','offFood'},'FontSize',15,'Location','northeastoutside')
title('tracking indices by food region','FontSize',15)