%% script extra names of all strains from the metadata file

[~,~,metadata] = xlsread('/Volumes/behavgenom$/Serena/aggregationScreeningDocs/aggregation_metadata.xlsx');
strainNames = {metadata{2:end,7}};
uniqueStrainNames = unique(strainNames);
strains = [{uniqueStrainNames{1:177}} {uniqueStrainNames{179:end}}]'; % remove "NONE"
save('strainsList/all.mat','strains');