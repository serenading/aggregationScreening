%% This is a one-off use script to divide the full aggregation dataset up in different ways
%% in order to accommodate copying onto separate hard drives
% author: serenading. May 2020

clear 
close all

%% Import metadata
metadataTable = readtable('/Users/sding/OneDrive - Imperial College London/aggScreening/source/metadata_aggregationScreening.csv','Delimiter',',');

%% Get five worm files
fiveWormLogInd = metadataTable.wormNum==5 & metadataTable.is_bad ==0;
dirnames = metadataTable.dirname(fiveWormLogInd);
basenames = metadataTable.basename(fiveWormLogInd);
filenames = cellfun(@(x,y) [x '/' y],dirnames,basenames,'un',0);
writecell(filenames,'/Users/sding/Desktop/AggregationScreening/divvy/fivewormfilenames.csv')

%% Get forty worm files
fortyWormLogInd = metadataTable.wormNum==40 & metadataTable.is_bad ==0;
dirnames = metadataTable.dirname(fortyWormLogInd);
basenames = metadataTable.basename(fortyWormLogInd);
filenames = cellfun(@(x,y) [x '/' y],dirnames,basenames,'un',0);
writecell(filenames,'/Users/sding/Desktop/AggregationScreening/divvy/fortywormfilenames.csv')

%% get control files logical index
N2LogInd = strcmp(metadataTable.strain_name,'N2');
CB4856LogInd = strcmp(metadataTable.strain_name,'CB4856');
DA609LogInd = strcmp(metadataTable.strain_name,'DA609');
controlLogInd = N2LogInd | CB4856LogInd | DA609LogInd;

%% get forty worm control files
fortywormControlLogInd = fortyWormLogInd & controlLogInd;
dirnames = metadataTable.dirname(fortywormControlLogInd);
basenames = metadataTable.basename(fortywormControlLogInd);
filenames = cellfun(@(x,y) [x '/' y],dirnames,basenames,'un',0);
writecell(filenames,'/Users/sding/Desktop/AggregationScreening/divvy/fortywormcontrolfilenames.csv')

%% get forty worm non-control files
fortywormNonControlLogInd = fortyWormLogInd & ~controlLogInd;
dirnames = metadataTable.dirname(fortywormNonControlLogInd);
basenames = metadataTable.basename(fortywormNonControlLogInd);
filenames = cellfun(@(x,y) [x '/' y],dirnames,basenames,'un',0);
writecell(filenames,'/Users/sding/Desktop/AggregationScreening/divvy/fortywormnoncontrolfilenames.csv')

%% get forty worm files but capping n = 15 control files
% choose 15 random N2 files to keep
N2LogIndPos = find(fortyWormLogInd & N2LogInd);
r = randperm(numel(N2LogIndPos),15);
N2Ind2Keep = N2LogIndPos(r);
N2LogInd2Keep = false(1,numel(N2LogInd));
N2LogInd2Keep(N2Ind2Keep) = true;
% choose 15 random CB4856 files to keep
CB4856LogIndPos = find(fortyWormLogInd & CB4856LogInd);
r = randperm(numel(CB4856LogIndPos),15);
CB4856Ind2Keep = CB4856LogIndPos(r);
CB4856LogInd2Keep = false(1,numel(CB4856LogInd));
CB4856LogInd2Keep(CB4856Ind2Keep) = true;
% choose 15 random DA609 files to keep
DA609LogIndPos = find(fortyWormLogInd & DA609LogInd);
r = randperm(numel(DA609LogIndPos),15);
DA609Ind2Keep = DA609LogIndPos(r);
DA609LogInd2Keep = false(1,numel(DA609LogInd));
DA609LogInd2Keep(DA609Ind2Keep) = true;
% get control logical index keeping 15 files from each of the three strains
controlSubsampledLogInd = N2LogInd2Keep | CB4856LogInd2Keep | DA609LogInd2Keep;
assert(nnz(controlSubsampledLogInd) == 15*3, 'total number of control files is not 45');
% get forty worm filenames keeping max 15 files from each of the three
% control strains
fortywormCappedControlLogInd = fortyWormLogInd & (~controlLogInd | controlSubsampledLogInd');
assert(nnz(fortywormCappedControlLogInd) == nnz(fortywormNonControlLogInd)+15*3)
% write filenames
dirnames = metadataTable.dirname(fortywormCappedControlLogInd);
basenames = metadataTable.basename(fortywormCappedControlLogInd);
filenames = cellfun(@(x,y) [x '/' y],dirnames,basenames,'un',0);
writecell(filenames,'/Users/sding/Desktop/AggregationScreening/divvy/fortywormcappedcontrolfilenames.csv')

%% get forty worm divergent files

% get divergent control files logical index for the 10 strains that are not CB4856 or N2
divergentLogInd = strcmp(metadataTable.strain_name,'CX11314') | strcmp(metadataTable.strain_name,'DL238')|...
strcmp(metadataTable.strain_name,'ED3017') | strcmp(metadataTable.strain_name,'EG4725')|...
strcmp(metadataTable.strain_name,'JT11398') | strcmp(metadataTable.strain_name,'JU258')|...
strcmp(metadataTable.strain_name,'JU775') | strcmp(metadataTable.strain_name,'LKC34')|...
strcmp(metadataTable.strain_name,'MY16') | strcmp(metadataTable.strain_name,'MY23');
% get forty worm divergent files
fortywormDivergentControlLogInd = fortyWormLogInd & (controlLogInd|divergentLogInd);
dirnames = metadataTable.dirname(fortywormDivergentControlLogInd);
basenames = metadataTable.basename(fortywormDivergentControlLogInd);
filenames = cellfun(@(x,y) [x '/' y],dirnames,basenames,'un',0);
writecell(filenames,'/Users/sding/Desktop/AggregationScreening/divvy/fortywormdivergentfilenames.csv')

%% get forty worm non-divergent files
fortywormNonDivergentLogInd = fortyWormLogInd & ~fortywormDivergentControlLogInd;
dirnames = metadataTable.dirname(fortywormNonDivergentLogInd);
basenames = metadataTable.basename(fortywormNonDivergentLogInd);
filenames = cellfun(@(x,y) [x '/' y],dirnames,basenames,'un',0);
writecell(filenames,'/Users/sding/Desktop/AggregationScreening/divvy/fortywormnondivergentfilenames.csv')