%clear
%close all

if ~exist('strainsList/copiedFileList.mat')
    %         addpath('auxiliary/')
    %         copiedFileList = rdir('/Volumes/diskAshurDT/behavgenom_archiv/AggregationScreening/MaskedVideos/*/*.hdf5');
    %         save('strainsList/copiedFileList.mat','copiedFileList')
    % Load the name of copied files and analysis featureTable
    load('strainsList/copiedFileList.mat','copiedFileList')
    featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fullFeaturesTable_20200519_153722.csv'],'Delimiter',',','preserveVariableNames',true);
    copiedFileLogInd = false(size(featureTable,1),1);
    notUsedFileCtr = 1;
    
    % go through each copied file to get copiedFile logical index
    for fileCtr=1:length(copiedFileList)
        filename = copiedFileList(fileCtr).name;
        filenamesplit = strsplit(filename,'/');
        filename = filenamesplit{end};
        fileIdx = find(strcmp(filename,featureTable.filename));
        if isempty(fileIdx)
            notUsedFileList{notUsedFileCtr} = filename;
            notUsedFileCtr = notUsedFileCtr+1;
        end
        copiedFileLogInd(fileIdx) = true;
    end
    
    % get info about copied files
    copiedFileMetadata = featureTable(copiedFileLogInd,1:17);
    save('strainsList/copiedFileList.mat','copiedFileList','copiedFileMetadata')
else
    load('strainsList/copiedFileList.mat','copiedFileList','copiedFileMetadata')
end