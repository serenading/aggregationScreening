clear 
close all

%% script generates list of files for each strain from metadata and checks it against metametadata
%% last ran script on 11 Apr 2018 and all entries match between meta and metameta files. voila!

load('strainsList/all.mat')
allLists = getFileList(strains);
[~,~,metametadata] = xlsread('/Volumes/behavgenom$/Serena/aggregationScreeningDocs/aggregation_metametadata.xlsx','A1:W2176');

for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    metanum_40 = length(allLists.([strain 'List_40']));
    metanum_5 = length(allLists.([strain 'List_5']));
    metametanum_40 = 0; % initialise metameta count
    metametanum_5 = 0; % initialise metameta count
    nameMatchInd = find(strcmp({metametadata{:,1}},strain)); % check for strain name match
    for nameMatchCtr = 1:length(nameMatchInd) % go through each name match
        fileCtr = nameMatchInd(nameMatchCtr);
        if isa(metametadata{fileCtr,17},'double') % check that it's a valid recording i.e. run number is a double rather than string 'x'.
            if metametadata{fileCtr,21} == 40 % check for wormNum
                metametanum_40 = metametanum_40+1;
            elseif metametadata{fileCtr,21} == 5 % check for wormNum
                metametanum_5 = metametanum_5+1;
            end
        end
    end
    % check to see if metadata count matches metametadata count
    if metanum_40 == metametanum_40 & metanum_5 == metametanum_5
        disp([num2str(strainCtr) ' strains ok'])
    elseif metanum_40 ~= metametanum_40
            warning([strain 'List_40 does not match. ' num2str(metanum_40) ' in meta file and ' num2str(metametanum_40) ' in metametafile'])
        if metanum_5 ~= metametanum_5
            warning([strain 'List_5 does not match. ' num2str(metanum_5) ' in meta file and ' num2str(metametanum_5) ' in metametafile'])
        end
    elseif metanum_5 ~= metametanum_5
        warning([strain 'List_5 does not match. ' num2str(metanum_5) ' in meta file and ' num2str(metametanum_5) ' in metametafile'])
    else
        warning(['something is not matching for strain ' strain])
    end
end
            
            