clear
close all

%% script counts the average number of worms in each video and generates a file that lists this actual number against expected 40/5 annotation

load('strainsList/all.mat')
% get list of file names for each strain
[strainFileList] = getFileList(strains);
strainWormNum = cell(0,3); % create empty cell to hold wormNum information

%% go through each strain
for strainCtr = 160:length(strains)
    filenames = vertcat(strainFileList.([strains{strainCtr} 'List_40']), strainFileList.([strains{strainCtr} 'List_5']));
    %% go through each recording
    for fileCtr = 1:length(filenames)
        
        %% load data
        filename = filenames{fileCtr};
        if  ~contains(filename, '/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_16.')...
          & ~contains(filename, '/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_4.3')...
          & ~strcmp(filename, '/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_11.3_180218/11.3_7_eca259_none_a1_Set0_Pos0_Ch3_18022018_170514_skeletons.hdf5')...
          & ~strcmp(filename, '/Volumes/behavgenom_archive$/Serena/AggregationScreening/Results/Agg_2.3_180119/2.3_6_nic261_f6_Set0_Pos0_Ch1_19012018_145323_skeletons.hdf5')
        % fill in filename information
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        end
    end
end