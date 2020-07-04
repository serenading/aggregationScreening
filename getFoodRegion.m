function [onFoodLogInd,foodEdgeLogInd,offFoodLogInd] = getFoodRegion(tsFeatures)

onFoodLogInd = tsFeatures.food_region == 1;
foodEdgeLogInd = tsFeatures.food_region == 0;
offFoodLogInd = tsFeatures.food_region == -1;