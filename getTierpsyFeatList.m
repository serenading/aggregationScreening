%% function generates a cell array with a list of Tierpsy features which contain a specified string;
% for (plate average features read by feature_stats = h5read(filename,'/features_stats/')

%% INPUT
% featString: string to be contained in the feature names
% featCSV: load from 'strainsList/featCSV.mat', which is identical to reading from the original features csv file as a Tierpsy output
%% OUTPUT
% featList: cell array with a list of Tierpsy features
% featPos: 1xn matrix with a list of Tierpsy feature indices for reading from /features_stats.values

function [featList, featPos] = getTierpsyFeatList(featString,featCSV)

% get feature stats name
feature_stats_name = {featCSV{1,1:end}};
if ~strcmp(featString,'Tierpsy_4548')
    if ~strcmp(featString,'Tierpsy_256')
        % find relevant features with plate average stats (there should be 89 "food" features, for example);
        featPos = find(~cellfun('isempty',strfind(cellstr(feature_stats_name),featString)));
    else
        featPos = NaN(256,1);
        % load file
        fid = fopen('strainsList/Tierpsy_256.txt');
        featnameString = 'initialise';
        fileCtr = 0;
        % read text file one line at a time
        while ischar(featnameString)
            % move onto the next line
            featnameString = fgetl(fid);
            fileCtr = fileCtr+1;
            featInd = find(~cellfun('isempty',strfind(cellstr(feature_stats_name),featnameString)));
            for featCtr = 1:numel(featInd)
                if cellfun('length',feature_stats_name(featInd(featCtr))) == numel(featnameString)
                    featPos(fileCtr,1) = featInd(featCtr);
                    break
                end
            end
        end
        fclose(fid);
    end
    % generate cell array to contain feature names
    featList = cell(numel(featPos),1);
    % add feature name to array
    for featCtr = 1:numel(featPos)
        feat = featPos(featCtr);
        featList{featCtr} = feature_stats_name(feat); % get the name of the food-related feature
    end
else
    featList = {feature_stats_name{2:end}}';
    featPos = [1:numel(featList)];
end