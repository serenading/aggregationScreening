clear
close all

%% script combines individual .txt files within a mapping directory to generate a single feature matrix,
% which is then exported as a csv file and ready to be used for cegwas mapping using R

% set directory
directory = '/Users/sding/Documents/AggScreening/results/mapping/swDensityCompare/Tierpsy/Tierpsy_256/';

addpath('auxiliary/')

% list all .txt files for mapping
[fileList, ~] = dirSearch(directory, '.txt');
numFiles = numel(fileList);

% pre-allocate
featMat = cell(numFiles, 198);

% go through each file
for fileCtr = 1:numFiles
    filename = fileList{fileCtr};
    % get feature name
    filenameSplit = strsplit(filename,'/');
    featName = filenameSplit{end};
    featNameSplit = strsplit(featName,'.');
    featName = featNameSplit{1};
    % get file content and format
    fid = fopen(filename);
    rawContent = textscan(fid,'%s'); % grab content
    formatContent = rawContent{1};
    formatContent = reshape(formatContent,[numel(2,formatContent)/2]);
    % add feature name
    formatContent{2,1} = featName;
    % add feature to feature matrix
    featMat(fileCtr,:) = {formatContent{2,:}};
end

% add strain names to top row
featMat = vertcat({formatContent{1,:}},featMat)';
% save file
saveFileName = [directory 'cegwasMappingFeatures.csv'];
cell2csv(saveFileName,featMat);
% the resultant csv file is ready for cegwas mapping using the R package