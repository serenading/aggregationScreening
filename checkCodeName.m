clear all

% script takes valid file names and that the two-character strain codes
% match strain names

%% load recording metadata
[~,~,metadata] = xlsread('/Volumes/behavgenom_archive$/Serena/AggregationScreening/aggregation_data.xlsx','A2:J2213');
metaNameList = {metadata{:,7}};
metaCodeList = {metadata{:,8}};
for codeCtr = 1:length(metaCodeList) % convert all numeric strain codes to strings
    if isa(metaCodeList{codeCtr},'double')
        metaCodeList{codeCtr} = num2str(metaCodeList{codeCtr});
    end
end

addpath('auxiliary/')

%% load strain code information
% load code reference
[~,~,strainCodeRefs] = xlsread('/Volumes/behavgenom_archive$/Serena/AggregationScreening/codes_CeNDR.xlsx','B1:B202');
% generate cell arrays containing 
strainCodes = cell(length(strainCodeRefs),2);
for strainCtr = 1:length(strainCodeRefs)
    strainCodeRef = strsplit(strainCodeRefs{strainCtr},' ');
    strainCodes{strainCtr,1} = strainCodeRef{1};
    strainCodes{strainCtr,2} = strainCodeRef{2};
end

%% check that the codes in the metadata are appropriate for the strain
% create a cell array to hold problematic file names
problemCodeNames = cell(100,5);
problemCodeNames{1,1} = 'refName';
problemCodeNames{1,2} = 'refCode';
problemCodeNames{1,3} = 'metaName';
problemCodeNames{1,4} = 'metaCode';
problemCodeNames{1,5} = 'metaBaseName';
problemCtr = 2;
% go through each recording file
for fileCtr = 1:length(metaNameList)
    metaName = metaNameList{fileCtr};
    if strcmp(metaName, 'NONE')
        refIdx = 1;
    else
        refIdx = find(ismember(strainCodes,metaName)); % find the index for the strain name in question
    end
    if strcmp(metaCodeList{fileCtr},strainCodes{refIdx,2}) % compare the metadata code to the reference code
        disp(['file ' num2str(fileCtr) ' ok'])
    else
        % record the incorrect files
        problemCodeNames{problemCtr,1} = strainCodes{refIdx,1};
        problemCodeNames{problemCtr,2} = strainCodes{refIdx,2};
        problemCodeNames{problemCtr,3} = metadata{fileCtr,7};
        problemCodeNames{problemCtr,4} = metadata{fileCtr,8};
        problemCodeNames{problemCtr,5} = metadata{fileCtr,1}; 
        problemCtr = problemCtr+1;
    end
end

%% save the list of problematic file names
dlmcell('strainList/problematicStrainCodes.txt',problemCodeNames);