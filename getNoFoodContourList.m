clear
close all

%% script find the list of files with no food contour

strainSet = 'all'; % 'controls','divergent','all'
wormNum = 5;
load(['strainsList/' strainSet '.mat'])
[strainFileList,fortyNum,fiveNum] = getFileList(strains);
noFoodContourFiles_skel = cell(1000,1);
noFoodContourCtr_skel = 1;
noFoodContourFiles_feat = cell(1000,1);
noFoodContourCtr_feat = 1;
noSkeletonFiles = cell(1000,1);
noSkeletonCtr = 1;

addpath('auxiliary/')

%% go through each strain
for strainCtr = 1:length(strains)
    filenames = strainFileList.([strains{strainCtr} 'List_' num2str(wormNum)]);
    %% go through each recording
    for fileCtr = 1:length(filenames)
        %% load data
        filename = filenames{fileCtr};
%         try
%             foodContourCoords = h5read(filename,'/food_cnt_coord');
%         catch
%             noFoodContourFiles_skel{noFoodContourCtr_skel} = filename;
%             noFoodContourCtr_skel = noFoodContourCtr_skel+1;
%         end
        filename = strrep(filenames{fileCtr},'_skeletons','_featuresN');
        try
            foodContourCoords = h5read(filename,'/food_cnt_coord');
        catch
            noFoodContourFiles_feat{noFoodContourCtr_feat} = filename;
            noFoodContourCtr_feat = noFoodContourCtr_feat+1;
        end
    end
end
% remove empty cells
% noFoodContourFiles_skel = noFoodContourFiles_skel(~cellfun('isempty',noFoodContourFiles_skel));
noFoodContourFiles_feat = noFoodContourFiles_feat(~cellfun('isempty',noFoodContourFiles_feat));

% save file names
%  dlmcell(['noFoodContourFiles_skel_' num2str(wormNum) 'new.txt'],noFoodContourFiles_skel);
 dlmcell(['noFoodContourFiles_feat_' num2str(wormNum) 'new.txt'],noFoodContourFiles_feat);

% delete bad feature files
for fileCtr = 1:length(noFoodContourFiles_feat)
    fileCtr
    filename = noFoodContourFiles_feat{fileCtr};
    delete(filename)
end
    