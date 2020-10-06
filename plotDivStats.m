clear
close all

%% Script reads in featuresTable and plots stacked bar graphs for blob category and food region statistics to show comprising fractions of tracking indices.
% Author: @serenading. July 2020

%% Specify which dataset
wormNum = 5;
extractStamp = '20200519_153722';

%% load
[featureTable,pathname] = loadLatestFeatureTable(extractStamp,wormNum);
uniqueStrains = unique(featureTable.strain_name);

%% pre-allocate
% all worm category
wormFractionsMean = NaN(numel(uniqueStrains),4);
wormFractionsStd = NaN(numel(uniqueStrains),4);
% phase 1 worm category
phase1WormFractionsMean = NaN(numel(uniqueStrains),4);
% not phase 1 worm category
notPhase1WormFractionsMean = NaN(numel(uniqueStrains),4);
% all worms food region
foodRegionFractionsMean = NaN(numel(uniqueStrains),3);
foodRegionFractionsStd = NaN(numel(uniqueStrains),3);
% lone worms food region
swFoodRegionFractionsMean = NaN(numel(uniqueStrains),3);
% paused multiworm food region
pausedMwFoodRegionFractionsMean = NaN(numel(uniqueStrains),3);
% cluster food region
clusterFoodRegionFractionsMean = NaN(numel(uniqueStrains),3);
% phase 1 food region
phase1FoodRegionFractionsMean = NaN(numel(uniqueStrains),3);
% not phase 1 food region
notPhase1FoodRegionFractionsMean = NaN(numel(uniqueStrains),3);

%% go through each strain to get worm fractions and food region fractions

% each strain
for strainCtr = 1:numel(uniqueStrains)
    strain = uniqueStrains(strainCtr);
    strainLogInd = strcmp(featureTable.strain_name,strain);
    
    % all worm category
    wormFractions = [featureTable.sw_fraction(strainLogInd),featureTable.pausedMw_fraction(strainLogInd),featureTable.cluster_fraction(strainLogInd),featureTable.tempBlob_fraction(strainLogInd)];
    wormFractionsMean(strainCtr,:) = mean(wormFractions,1);
    wormFractionsSE(strainCtr,:) = std(wormFractions)/sqrt(nnz(strainLogInd));
        
    % phase 1 worm category
    phase1WormFractions = [featureTable.sw_phase1_fraction(strainLogInd),featureTable.pausedMw_phase1_fraction(strainLogInd),featureTable.cluster_phase1_fraction(strainLogInd),featureTable.tempBlob_phase1_fraction(strainLogInd)]...
        ./featureTable.phase1_fraction(strainLogInd);
    phase1WormFractionsMean(strainCtr,:) = mean(phase1WormFractions,1);
    
    % not phase 1 worm category
    notPhase1WormFractions = [featureTable.sw_notPhase1_fraction(strainLogInd),featureTable.pausedMw_notPhase1_fraction(strainLogInd),featureTable.cluster_notPhase1_fraction(strainLogInd),featureTable.tempBlob_notPhase1_fraction(strainLogInd)]...
        ./featureTable.notPhase1_fraction(strainLogInd);
    notPhase1WormFractionsMean(strainCtr,:) = mean(notPhase1WormFractions,1);
    
    % all worms food region
    foodRegionFractions = [featureTable.food_region_inside_fraction(strainLogInd),featureTable.food_region_edge_fraction(strainLogInd),featureTable.food_region_outside_fraction(strainLogInd)];
    foodRegionFractionsMean(strainCtr,:) = mean(foodRegionFractions,1);
    foodRegionFractionsSE(strainCtr,:) = std(foodRegionFractions)/sqrt(nnz(strainLogInd));
    
    % single worm food region 
    swFoodRegionFractions = [featureTable.sw_onFood_fraction(strainLogInd),featureTable.sw_foodEdge_fraction(strainLogInd),featureTable.sw_offFood_fraction(strainLogInd)]...
        ./featureTable.sw_fraction(strainLogInd);
    swFoodRegionFractionsMean(strainCtr,:) = mean(swFoodRegionFractions,1);
    
    % paused multiworm food region
    pausedMwFoodRegionFractions = [featureTable.pausedMw_onFood_fraction(strainLogInd),featureTable.pausedMw_foodEdge_fraction(strainLogInd),featureTable.pausedMw_offFood_fraction(strainLogInd)]...
        ./featureTable.pausedMw_fraction(strainLogInd);
    pausedMwFoodRegionFractionsMean(strainCtr,:) = mean(pausedMwFoodRegionFractions,1);
    
    % cluster food region
    clusterFoodRegionFractions = [featureTable.cluster_onFood_fraction(strainLogInd),featureTable.cluster_foodEdge_fraction(strainLogInd),featureTable.cluster_offFood_fraction(strainLogInd)]...
        ./featureTable.cluster_fraction(strainLogInd);
    clusterFoodRegionFractionsMean(strainCtr,:) = mean(clusterFoodRegionFractions,1);
    
    % phase 1 food region
    phase1FoodRegionFractions = [featureTable.onFood_phase1_fraction(strainLogInd),featureTable.foodEdge_phase1_fraction(strainLogInd),featureTable.offFood_phase1_fraction(strainLogInd)]...
        ./featureTable.phase1_fraction(strainLogInd);
    phase1FoodRegionFractionsMean(strainCtr,:) = mean(phase1FoodRegionFractions,1);
    
    % not phase 1 food rgion
    notPhase1FoodRegionFractions = [featureTable.onFood_notPhase1_fraction(strainLogInd),featureTable.foodEdge_notPhase1_fraction(strainLogInd),featureTable.offFood_phase1_fraction(strainLogInd)]...
        ./featureTable.notPhase1_fraction(strainLogInd);
    notPhase1FoodRegionFractionsMean(strainCtr,:) = mean(notPhase1FoodRegionFractions,1);
end

%% stacked bar plots to show proportions

% all worm category
figure;bar(wormFractionsMean,'stacked')
% hold on;errorbar(cumsum(wormFractionsMean')',wormFractionsSE,'.k');
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'single worm','paused multiworm','cluster','temp'},'FontSize',15,'Location','northeastoutside')
title('tracking indices by worm category','FontSize',15)

% phase 1 worm category
figure;bar(phase1WormFractionsMean,'stacked')
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'single worm','paused multiworm','cluster','temp'},'FontSize',15,'Location','northeastoutside')
title('phase 1 tracking indices by worm category','FontSize',15)

% not phase 1 worm category
figure;bar(notPhase1WormFractionsMean,'stacked')
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'single worm','paused multiworm','cluster','temp'},'FontSize',15,'Location','northeastoutside')
title('not phase 1 tracking indices by worm category','FontSize',15)

% all worms food region
figure;bar(foodRegionFractionsMean,'stacked')
% hold on;errorbar(cumsum(foodRegionFractionsMean')',foodRegionFractionsSE,'.k');
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylim([0 1])
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'onFood','foodEdge','offFood'},'FontSize',15,'Location','northeastoutside')
title('tracking indices by food region','FontSize',15)

% single worms food region
figure;bar(swFoodRegionFractionsMean,'stacked')
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylim([0 1])
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'onFood','foodEdge','offFood'},'FontSize',15,'Location','northeastoutside')
title('single worm tracking indices by food region','FontSize',15)

% paused multiworms food region
figure;bar(pausedMwFoodRegionFractionsMean,'stacked')
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylim([0 1])
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'onFood','foodEdge','offFood'},'FontSize',15,'Location','northeastoutside')
title('paused multiworm tracking indices by food region','FontSize',15)

% cluster food region
figure;bar(clusterFoodRegionFractionsMean,'stacked')
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylim([0 1])
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'onFood','foodEdge','offFood'},'FontSize',15,'Location','northeastoutside')
title('cluster tracking indices by food region','FontSize',15)

% phase 1 food region
figure;bar(phase1FoodRegionFractionsMean,'stacked')
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylim([0 1])
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'onFood','foodEdge','offFood'},'FontSize',15,'Location','northeastoutside')
title('phase 1 tracking indices by food region','FontSize',15)

% not phase 1 food region
figure;bar(notPhase1FoodRegionFractionsMean,'stacked')
xticks(1:numel(uniqueStrains))
xticklabels(uniqueStrains)
xtickangle(90)
ylim([0 1])
ylabel(['fraction of tracking indices for ' num2str(wormNum) ' worms'],'FontSize',15)
legend({'onFood','foodEdge','offFood'},'FontSize',15,'Location','northeastoutside')
title('not phase 1 tracking indices by food region','FontSize',15)