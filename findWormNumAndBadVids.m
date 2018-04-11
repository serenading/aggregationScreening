clear

%% script 1) generates unique file identifier from the metadata file, 2) finds the number of worms in each recording from the metametadata file,
%% and 3) finds bad videos (either extra videos or 'NONE' videos). It then updates the metadata file with all these new information. 

%% load files
[~,~,metadata] = xlsread('/Volumes/behavgenom$/Serena/aggregationScreeningDocs/aggregation_metadata.xlsx');
[~,~,metametadata] = xlsread('/Volumes/behavgenom$/Serena/aggregationScreeningDocs/aggregation_metametadata.xlsx','A1:W2180');

% reformat data files
metadata = metadata(2:end,:); % remove the row with column labels

%% go through each row of the metadata and metametadata file to unique generate session number

metaFileID = cell(size(metadata,1),1);
for fileCtr = 1:size(metadata,1)
    % generate session number by combining block number and day number
    sessionN = [num2str(metadata{fileCtr,4}) '.' num2str(metadata{fileCtr,5})];
    % generate unique fileID in the format of sessionN_runN_camN
    metaFileID{fileCtr} = [sessionN '_' num2str(metadata{fileCtr,6}) '_' num2str(metadata{fileCtr,9})];
end

metametaFileID = cell(size(metametadata,1),1);
for fileCtr = 1:size(metametadata,1)
    % generate unique fileID in the format of sessionN_runN_camN
    metametaFileID{fileCtr} = [num2str(metametadata{fileCtr,16}) '_' num2str(metametadata{fileCtr,17}) '_' num2str(metametadata{fileCtr,19})];
end

%% search for worm number based on fileID

% create new variables to hold wormNum and to contain extra/bad video flag
wormNum = NaN(size(metadata,1),1);

% create a file to hold "extraFiles" i.e. bad recordings that we keep in case we want to look at them later (though unlikely)
extraFiles = cell(190,3);
extraFiles{1,1} = 'identifier';
extraFiles{1,2} = 'metaIndex';
extraFiles{1,3} = 'strainName';
extraFileCtr = 2;

% go through metadata recordings
for fileCtr = 1:size(metadata,1)
    thismetaFileID = metaFileID{fileCtr}; % get the unique identifier from the metadata file
    metametaIdx = find(strcmp(metametaFileID,thismetaFileID)); % find the index of the metameta file matching this identifier
    % check that there is only one file identifier match
    if length(metametaIdx) ==1 % exact match by identifier
        % check that the strain names match
        if strcmp(metadata{fileCtr,7},metametadata{metametaIdx,1})
            wormNum(fileCtr) = metametadata{metametaIdx,21};
        elseif ~strcmp(metadata{fileCtr,7},'NONE')
            warning(['file ' thismetaFileID ' with meta index ' num2str(fileCtr) ' and metameta index ' num2str(metametaIdx) ' has mismatched strain names'])
        end
    elseif length(metametaIdx) >1 % more than one matching identifier
        warning(['file ' thismetaFileID ' with meta index ' num2str(fileCtr) ' has ' num2str(length(metametaIdx)) ' identifier matches in metametadata: metametaIdx ' num2str(metametaIdx(1)) ' and '  num2str(metametaIdx(2))])
    elseif isempty(metametaIdx) % no matching identifier
        % write info to extraFiles file to keep track
        extraFiles{extraFileCtr,1} = thismetaFileID;
        extraFiles{extraFileCtr,2} = fileCtr;
        extraFiles{extraFileCtr,3} = metadata{fileCtr,7};
        extraFileCtr = extraFileCtr+1;
        % change flag for is_bad
        metadata{fileCtr,3} = true;
    end
end

%% add fileID and wormNum to metadata file
metadata = [metadata,metaFileID];
metadata = [metadata,num2cell(wormNum)];
[~,~,metadataLabel] = xlsread('/Volumes/behavgenom$/Serena/aggregationScreeningDocs/aggregation_metadata.xlsx','A1:J1');
metadataLabel{11} = 'metafileID';
metadataLabel{12} = 'wormNum';
metadata = vertcat(metadataLabel,metadata);

%% export updated metadata file and save extraFiles as a text document
%dlmcell('metadata.txt',metadata)
%dlmcell('extraFilesList.txt',extraFiles)