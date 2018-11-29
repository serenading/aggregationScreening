clear
close all

%% script works with bright field wild isolate aggregation screening dataset (40 and 5 worms) to compare single worm food-related Tierpsy features
% note: phase restriction is not possible with plate average stats as they are pre-calculated for the entire duration of the movie

%% set parameters
% set analysis parameters
strainSet = 'all'; % 'controls','divergent','all'
maxNumReplicates = 60; % controls have up to 60 reps, divergents up to 15 reps, all other strains up to 5 reps.
wormNums = {'40','5'};
featString = 'width';
saveResults = true;

%% prep work
% addpath('auxiliary/')
% load the strain names included in the specified strainSet
load(['strainsList/' strainSet '.mat'])
% get list of file names for each strain
[strainFileList,~,~] = getFileList(strains);

%% generate a list of food-related Tierpsy features
[featList, featPos] = getTierpsyFeatList(featString);

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
        allFeats.(strain) = cell(numel(featList),length(filenames)+1);
        %% go through each recording
        for fileCtr = 1:length(fileInd)
            %% load data
            filename = strrep(filenames{fileInd(fileCtr)},'skeletons','featuresN');
            if ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_11.1_180216/11.1_6_eca259_a1_ju2522_94_Set0_Pos0_Ch5_16022018_144346_featuresN.hdf5')...
                & ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/agg_7.4_180209/7.4_1_ju778_7b_Set0_Pos0_Ch1_09022018_155743_featuresN.hdf5')...
                & ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_11.1_180216/11.1_6_qg557_8d_Set0_Pos0_Ch3_16022018_145004_featuresN.hdf5')...
                & ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/agg_1.2_180115/1.2_9_cb4856_oo_Set0_Pos0_Ch4_15012018_192913_featuresN.hdf5')...
                & ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_1.3_180116/1.3_6_eca246_8f_Set0_Pos0_Ch2_16012018_145913_featuresN.hdf5')...
                & ~strcmp(filename,'/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/agg_1.2_180115/1.2_5_ed3048_c4_Set0_Pos0_Ch4_15012018_141902_featuresN.hdf5')
                feature_stats = h5read(filename,'/features_stats');
                if numel(feature_stats.value) == 4688 % full features list should have 4688 features. only use feature files with the full list
                    for featCtr = 1:numel(featList)
                        allFeats.(strain){featCtr,1} = featList{featCtr}; % add feature name to the first column
                        allFeats.(strain){featCtr,fileCtr+1} = feature_stats.value(featPos(featCtr)); % add the corresponding feature stat to array
                    end
                else
                    disp([filename ' only has ' num2str(numel(feature_stats.value)) ' features'])
                end
            end
        end
    end
    if saveResults
        if strcmp(wormNum,'40')
            save(['results/TierpsyFeat_' featString '_' strainSet '.mat'])
        elseif strcmp(wormNum,'5')
            save(['results/5worm_TierpsyFeat_' featString '_' strainSet '.mat'],'allFeats')
        end
    end
end