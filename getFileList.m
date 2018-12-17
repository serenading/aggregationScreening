function [strainFileList,strainFileListMat] = getFileList(strains)

%% function takes a list of desired strain names and generates a struct containing lists of analysis files

%% INPUT: 
% strains: cell array containing separate strings for each of the desired strain i.e. strains = {'N2','DA609','CB4856'};

%% OUTPUT:
% strainFileList: a struct containing two fields for each strain, with each field containing a cell array 
% listing all the relevant file paths for 40 worm or 5 worm recordings.
% strainFileListMeta: a cell array containing strain name, worm number, and file path
% (identical to strainFileList, just different format i.e. matrix as opposed to struct)

% load metadata
[~,~,metadata] = xlsread('/Volumes/behavgenom$/Serena/aggregationScreeningDocs/aggregation_metadata.xlsx');
metadata = metadata(2:end,:); % remove the first row with labels

% preallocate cell arrays to hold path names
strainFileListMat = cell(2300,3); 
matCtr = 1;
for strainCtr = 1:length(strains)
    strainFileList.(genvarname([strains{strainCtr} 'List_40'])) = cell(80,1);
    strainFileList.(genvarname([strains{strainCtr} 'List_5'])) = cell(80,1);
    strainFileCtr.(genvarname([strains{strainCtr} 'List_40'])) = 1;
    strainFileCtr.(genvarname([strains{strainCtr} 'List_5'])) = 1;
end

% create a list of hdf5 file paths for each of the strains of interest
for fileCtr = 1:size(metadata,1)
    if metadata{fileCtr,3} == 0 % check the video doesn't have a bad flag (is_bad logical index is double class for some reason...)
        metaStrainName = metadata{fileCtr,7};
        for strainCtr = 1:length(strains)
            if strcmp(metaStrainName,strains{strainCtr}) % strain names match
                if  metadata{fileCtr,12} == 40 | metadata{fileCtr,12} == 5
                    % add info to strainFileListMeta array
                    strainFileListMat{matCtr,1} = metaStrainName;
                    strainFileListMat{matCtr,2} = num2str(metadata{fileCtr,12});
                    strainFileListMat{matCtr,3} = [strrep(metadata{fileCtr,2},'MaskedVideos','Results') '/' strrep(metadata{fileCtr,1},'.hdf5','_skeletons.hdf5')];
                    matCtr = matCtr+1; % update counter
                    % add info to strainFileList struct
                    if metadata{fileCtr,12} == 40
                        strainFileList.([strains{strainCtr} 'List_40']){strainFileCtr.([strains{strainCtr} 'List_40'])} =...
                            [strrep(metadata{fileCtr,2},'MaskedVideos','Results') '/' strrep(metadata{fileCtr,1},'.hdf5','_skeletons.hdf5')]; % add file path to list
                        strainFileCtr.([strains{strainCtr} 'List_40']) = strainFileCtr.([strains{strainCtr} 'List_40'])+1; % update counter
                    elseif metadata{fileCtr,12} == 5
                        strainFileList.([strains{strainCtr} 'List_5']){strainFileCtr.([strains{strainCtr} 'List_5'])} =...
                            [strrep(metadata{fileCtr,2},'MaskedVideos','Results') '/' strrep(metadata{fileCtr,1},'.hdf5','_skeletons.hdf5')]; % add file path to list
                        strainFileCtr.([strains{strainCtr} 'List_5']) = strainFileCtr.([strains{strainCtr} 'List_5'])+1; % update counter
                    end
                end
            end
        end
    end
end

% remove empty cells
strainFileListMat = strainFileListMat(~cellfun('isempty',strainFileListMat));
strainFileListMat = reshape(strainFileListMat,[numel(strainFileListMat)/3, 3]); % reshape
strainFileListMat = sort(strainFileListMat); % sort
for strainCtr = 1:length(strains)
    strainFileList.([strains{strainCtr} 'List_40']) = strainFileList.([strains{strainCtr} 'List_40'])(~cellfun('isempty',strainFileList.([strains{strainCtr} 'List_40'])));
    strainFileList.([strains{strainCtr} 'List_5']) = strainFileList.([strains{strainCtr} 'List_5'])(~cellfun('isempty',strainFileList.([strains{strainCtr} 'List_5'])));
end