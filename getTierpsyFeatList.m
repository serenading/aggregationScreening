%% function generates a cell array with a list of Tierpsy features which contain a specified string;
% for (plate average features read by feature_stats = h5read(filename,'/features_stats/')

%% INPUT
% featString: string to be contained in the feature names
%% OUTPUT
% featList: cell array with a list of Tierpsy features
% featPos: 1xn matrix with a list of Tierpsy feature indices for reading from /features_stats.values

function [featList, featPos] = getTierpsyFeatList(featString)

% open sample features.N file
filename = '/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/agg_7.2_180205/7.2_3_ab1_78_Set0_Pos0_Ch6_05022018_110819_featuresN.hdf5';
% read features
features = h5read(strrep(filename,'skeletons','featuresN'),'/timeseries_data/');
feature_stats = h5read(strrep(filename,'skeletons','featuresN'),'/features_stats/');
feature_stats_name = feature_stats.name';
% find relevant features with plate average stats (there should be 89 "food" features, for example);
featPos = find(~cellfun('isempty',strfind(cellstr(feature_stats_name),featString)));
% generate cell array to contain feature names
featList = cell(1,numel(featPos));
% add feature name to array
for featCtr = 1:numel(featPos)
    pos = featPos(featCtr);
    featList{featCtr} = feature_stats_name(pos,:); % get the name of the food-related feature
end
