clear
close all

% Script plots heritability scores that Erik and Loraina sent in Oct 2020
% and finds overlapping features with Tierpsy256.
% Author: serenading. Oct 2020. 

%% Set parameters
% what's the path to the heritability tsv file?
fpath = '/Users/sding/OneDrive - Imperial College London/aggScreening/source/AndersenLab/heritability_overlap_20201001.tsv';
% display top features or not?
dispFeats = true;
% how many top features to display?
n_feats_display = 30;

%% Load scores
scores = tdfread(fpath);

%% Plot
figure;

% h2
subplot(1,2,1)
histogram(scores.H2,'Normalization','count','BinWidth',0.05)
xlabel('H2')
ylabel('count')
xlim([0.25 1])
xticks([0.25:0.25:1])
ylim([0 250])
title('broad sense heritability')
% format
set(gca,'FontSize',15)

% H2
subplot(1,2,2)
histogram(scores.h2,'Normalization','count','BinWidth',0.05)
xlabel('h2')
ylabel('count')
xlim([0.25 1])
xticks([0.25:0.25:1])
ylim([0 250])
title('narrow sense heritability')
% format
set(gca,'FontSize',15)

%% Get top features
% features already sorted according to h2 value so no need to sort first
% get h2 features
feats_h2 = scores.trait;
if dispFeats
    % display top h2 features
    disp(['Top ' num2str(n_feats_display) ' h2 features are:'])
    feats_h2(1:n_feats_display,:)
    % display top 15 h2 scores
    disp(['Top ' num2str(n_feats_display) ' h2 scores are:'])
    scores.h2(1:n_feats_display)
end

% sort H2 to get top features
[sortH2, sortH2I] = sort(scores.H2,'descend');
% get H2 features
feats_H2 = scores.trait(sortH2I,:);
if dispFeats
    % display top 15 H2 features
    disp(['Top ' num2str(n_feats_display) ' H2 features are:'])
    feats_H2(1:n_feats_display,:)
    % display top 15 H2 scores
    disp(['Top ' num2str(n_feats_display) ' H2 scores are:'])
    sortH2(1:n_feats_display)
end

%% Compare top features
% Load in Tierpsy 256
feats_256Tierpsy = table2cell(readtable('strainsList/Tierpsy_256_short.csv','PreserveVariableNames',true,'ReadVariableNames',false));
% Load in mapped Tierpsy 256
feats_256TierpsyMapped = table2cell(readtable('strainsList/Tierpsy_256_cegwas2mapped.csv','ReadVariableNames',false));
feats_256TierpsyMapped = feats_256TierpsyMapped(:,1);
% Find features that are in both heritability table and in Tierpsy256
sharedFeats = intersect(cellstr(scores.trait),feats_256Tierpsy');
% Find features that are in both heritability table and in mapped Tierpsy256
sharedMappedFeats = intersect(cellstr(scores.trait),feats_256TierpsyMapped');

%% See whether a particular feature is contained in a particular set
feat = 'd_ang_vel_midbody_w_backward_abs_50th' %'d_ang_vel_head_base_w_forward_abs_90th'
ismember(feat,cellstr(scores.trait))
ismember(feat,feats_256TierpsyMapped)
ismember(feat,sharedFeats)
ismember(feat,sharedMappedFeats)