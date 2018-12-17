clear
close all

%% script takes all the realhits csv files within the specified directory and generate a summary text file within the same directory listing hit traits and variance explained by the trait

addpath('../auxiliary/')
% set directory
directory = 'results/';
% list all .txt files for mapping
[fileList, ~] = dirSearch(directory, 'width.csv');
numFiles = numel(fileList);
% go through each file
for fileCtr = 1:numFiles
    filename = fileList{fileCtr};
    % get content
    fileContent = readtable(filename);
    % slice out columns with traits and var.exp
    fileContent = table2cell(fileContent(:,[5,12]));
    % separate trait column from variation explained column
    traitCol = {fileContent{:,1}};
    varExpCol = [fileContent{:,2}];
    % find unique traits
    [traits,traitRowIdx,~] = unique(traitCol);
    % pre-allocate array
    results = cell(numel(traits),2);
    % get trait name and variation explained for each trait
    for traitCtr = 1:numel(traits)
        results{traitCtr,1} = traits{traitCtr};
        results{traitCtr,2} = varExpCol(traitRowIdx(traitCtr));
    end
    % add top row
    rowTitle = {'trait','variance explained'};
    results = vertcat(rowTitle,results);
    % save results as text file
    savename = strrep(filename,'.csv','_summary.txt');
    dlmcell(savename,results)
end