clear
close all

%% script counts the average number of worms in each video and generates a file that lists this actual number against expected 40/5 annotation

load('strainsList/all.mat')
% get list of file names for each strain
[strainFileList] = getFileList(strains);
strainWormNum = cell(0,3); % create empty cell to hold wormNum information

%% go through each strain
for strainCtr = 1:length(strains)
    filenames = vertcat(strainFileList.([strains{strainCtr} 'List_40']), strainFileList.([strains{strainCtr} 'List_5']));
    % generate empty array to hold recording and wormNum info
    wormNum = cell(length(filenames),3);
    % fill in expected wormNum info
    wormNum(:,2) = num2cell(vertcat(40*ones(length(strainFileList.([strains{strainCtr} 'List_40'])),1), 5*ones(length(strainFileList.([strains{strainCtr} 'List_5'])),1)),2);
    %% go through each recording
    for fileCtr = 1:length(filenames)
        
        %% load data
        filename = filenames{fileCtr};
        % fill in filename information
        wormNum{fileCtr,1} = filename;
        trajData = h5read(filename,'/trajectories_data');
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        numWormPerFrame = length(trajData.frame_number)/max(trajData.frame_number);
        % fill in actual wormNum info
        wormNum{fileCtr,3} = numWormPerFrame;
    end
    strainWormNum = vertcat(strainWormNum,wormNum);
end