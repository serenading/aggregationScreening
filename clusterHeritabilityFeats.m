clear
close all

%% Script probably obsolete.
% Script clusters Tierpsy feature values in strain x feature format, where
% each row is a unique strain and each column is one of the 641
% heritability features that Erik and Loraina sent over in Oct 2020.
% Author: serenading. Oct 2020

%% Load data
ftPath = '/Users/sding/OneDrive - Imperial College London/HFSP/5 worm locomotion dataset/fiveWormMappingFile_20200519_153722.txt'; % this is the file sent to Erik but instead of .tsv it's an identical .txt file so it works with readtable
featureTable = readtable(ftPath,'Delimiter','\t');
h2Path = '/Users/sding/OneDrive - Imperial College London/aggScreening/source/heritability_overlap_20201001.tsv';
h2Feats = cellstr(tdfread(h2Path).trait);

%% Trim features down to the 641 heritability feats
strain_name = featureTable.strain;
featureMat = table2array(featureTable(:,h2Feats));
featureTable.strain_name = strain_name;

%% Cluster
cgObj = clustergram(featureMat,'RowLabels',strain_name,'ColumnLabels',h2Feats,...
    'Colormap',redbluecmap,'ShowDendrogram','off','OptimalLeafOrder',true)

%% view group info
% right click on node and export group info as GroupInfo
GroupInfo.ColumnNodeNames
GroupInfo.RowNodeNames