clear
close all

%% script works with bright field wild isolate aggregation screening dataset (40 and 5 worms) to compare single worm food-related Tierpsy features
% note: phase restriction is not possible with plate average stats as they are pre-calculated for the entire duration of the movie

%% set parameters
% set analysis parameters
strainSet = 'all'; % 'controls','divergent','all'
maxNumReplicates = 60; % controls have up to 60 reps, divergents up to 15 reps, all other strains up to 5 reps.
wormNums = {'40','5'};
featString = 'Tierpsy_4548'; %'food','blob','eigen','width','length','area','axis','speed','velocity','curvature','Tierpsy_256','Tierpsy_4548','forward','backward','paused'
saveResults = true;

%% prep work

addpath('auxiliary/')
% load the strain names included in the specified strainSet
load(['strainsList/' strainSet '.mat'])
% get list of file names for each strain
[strainFileList,~] = getFileList(strains);
% load features csv files
load('strainsList/featCSV.mat'); % much faster to load saved cell array than to read directly from the csv file
[~,~,filenameCSV] = xlsread('/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/filenames_summary_tierpsy_plate_20181203_141111.xlsx');
% generate a list of required Tierpsy features
[featList, featPos] = getTierpsyFeatList(featString,featCSV);

%% go through each wormNum
for numCtr = 1:length(wormNums)
    wormNum = wormNums{numCtr};
    %% go through each strain
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        filenames = strainFileList.([strain 'List_' wormNum]);
        % if there are many files, then subsample recordings without replacement
        if length(filenames)>maxNumReplicates
            fileInd = datasample(1:length(filenames),maxNumReplicates,'Replace',false);
        else
            fileInd = 1:length(filenames);
        end
        %% initialise
        allFeats.(strain) = cell(numel(featList),length(fileInd)+1);
        %% go through each recording
        for fileCtr = 1:length(fileInd)
            % get file name
            filename = strrep(filenames{fileInd(fileCtr)},'skeletons','featuresN');
            %% find feature value for that filename
            % (Note that a series of manuvoures are necessary because filenameCSV and featCSV
            % have different sizes in the first dimension, so need to use file_id to connect them.)
            % find corresponding index in the filenameCSV
            filenameCSVnames = {filenameCSV{:,2}};
            fileIdxNameCSV = find(~cellfun('isempty',strfind(cellstr(filenameCSVnames),filename)));
            % find file_id which connects filenameCSV to featCSV
            file_id = cell2mat(filenameCSV(fileIdxNameCSV,1));
            % find corresponding index in the featCSV
            file_ids = [featCSV{2:end,1}]; % remove top row to keep matrix double
            fileIdxFeatCSV = find(file_ids==file_id)+1; % +1 to add back the row so row indices match
            % add feature name and value to allFeats variable to be saved
            for featCtr = 1:numel(featList)
                allFeats.(strain){featCtr,1} = featList{featCtr}; % add feature name to the first column
                allFeats.(strain){featCtr,fileCtr+1} = featCSV{fileIdxFeatCSV,featPos(featCtr)}; % add the corresponding feature stat
            end
        end
    end
    %% save results
    if saveResults
        if strcmp(wormNum,'40')
            save(['results/TierpsyFeat_' featString '_' strainSet '.mat'],'allFeats')
        elseif strcmp(wormNum,'5')
            save(['results/5worm_TierpsyFeat_' featString '_' strainSet '.mat'],'allFeats')
        end
    end
end