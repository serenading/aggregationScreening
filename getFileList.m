function [strainFileList] = getFileList(strains)

%% function takes a list of desired strain names and generates a struct containing lists of analysis files

% INPUT: strains: cell array containing separate strings for each of the
% desired strain i.e. strains = {'N2','DA609','CB4856'};
% OUTPUT: strainFileList: a struct containing two fields for each strain,
% with each field containing a cell array listing all the relevant file
% paths for 40 worm or 5 worm recordings. 

% load metadata
[~,~,metadata] = xlsread('/Volumes/behavgenom_archive$/Serena/AggregationScreening/aggregation_data.xlsx','A2:J2213');
metaStrainList = {metadata{1:end,7}};

% preallocate cell arrays to hold path names
for strainCtr = 1:length(strains)
    strainFileList.(genvarname([strains{strainCtr} 'List_40'])) = cell(80,1);
    strainFileList.(genvarname([strains{strainCtr} 'List_5'])) = cell(80,1);
    strainFileCtr.(genvarname([strains{strainCtr} 'List_40'])) = 1;
    strainFileCtr.(genvarname([strains{strainCtr} 'List_5'])) = 1;
end

% create a list of hdf5 file paths for each of the strains of interest
for metaCtr = 1:length(metaStrainList)
    metaStrainName = metaStrainList{metaCtr};
    for strainCtr = 1:length(strains)
        if strcmp(metaStrainName,strains{strainCtr}) % strain names match
            if ~strcmp(strains{strainCtr},'DA609') & length(strsplit(metadata{metaCtr,1},'_')) == 9
                if mod(metadata{metaCtr,9},2)==1 % if camera pos is odd
                    strainFileList.([strains{strainCtr} 'List_40']){strainFileCtr.([strains{strainCtr} 'List_40'])} =...
                        [strrep(metadata{metaCtr,2},'MaskedVideos','Results') '/' strrep(metadata{metaCtr,1},'.hdf5','_skeletons.hdf5')]; % add file path to list
                    strainFileCtr.([strains{strainCtr} 'List_40']) = strainFileCtr.([strains{strainCtr} 'List_40'])+1; % update counter
                else
                    strainFileList.([strains{strainCtr} 'List_5']){strainFileCtr.([strains{strainCtr} 'List_5'])} =...
                        [strrep(metadata{metaCtr,2},'MaskedVideos','Results') '/' strrep(metadata{metaCtr,1},'.hdf5','_skeletons.hdf5')]; % add file path to list
                    strainFileCtr.([strains{strainCtr} 'List_5']) = strainFileCtr.([strains{strainCtr} 'List_5'])+1; % update counter
                end
            elseif strcmp(strains{strainCtr},'DA609') & length(strsplit(metadata{metaCtr,1},'_')) == 8
                if mod(metadata{metaCtr,9},2)==1 % if camera pos is odd
                    strainFileList.([strains{strainCtr} 'List_40']){strainFileCtr.([strains{strainCtr} 'List_40'])} =...
                        [strrep(metadata{metaCtr,2},'MaskedVideos','Results') '/' strrep(metadata{metaCtr,1},'.hdf5','_skeletons.hdf5')]; % add file path to list
                    strainFileCtr.([strains{strainCtr} 'List_40']) = strainFileCtr.([strains{strainCtr} 'List_40'])+1; % update counter
                else
                    strainFileList.([strains{strainCtr} 'List_5']){strainFileCtr.([strains{strainCtr} 'List_5'])} =...
                        [strrep(metadata{metaCtr,2},'MaskedVideos','Results') '/' strrep(metadata{metaCtr,1},'.hdf5','_skeletons.hdf5')];% add file path to list
                    strainFileCtr.([strains{strainCtr} 'List_5']) = strainFileCtr.([strains{strainCtr} 'List_5'])+1; % update counter
                end
                % ignoring the cases where file names contain two strains for now
            end
        end
    end
end

% remove empty cells
for strainCtr = 1:length(strains)
    strainFileList.([strains{strainCtr} 'List_40']) = strainFileList.([strains{strainCtr} 'List_40'])(~cellfun('isempty',strainFileList.([strains{strainCtr} 'List_40'])));
    strainFileList.([strains{strainCtr} 'List_5']) = strainFileList.([strains{strainCtr} 'List_5'])(~cellfun('isempty',strainFileList.([strains{strainCtr} 'List_5'])));
end