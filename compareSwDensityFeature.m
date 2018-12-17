clear
close all

%% Script reads saved feature data (Tierpsy sw features or other calculated features)
% to assess density-dependent (40 vs. 5 worms) effects on single worm features.
% Script exports relevant feature for mapping.


%% set parameters
strainSet = 'all'; % 'controls','divergent','all'
feature = 'Tierpsy_4548'; % 'speed_nonTierpsy','food','blob','eigen','width','length','area','axis','speed','velocity','curvature','Tierpsy_256','Tierpsy_4548'
applyBonferroniCorrection = true;
saveResults = true;

%% default parameters
useTierpsyFeat = true;
if strcmp(feature, 'speed_nonTierpsy')
    useTierpsyFeat = false;
end

if useTierpsyFeat
    metricRange = 4; % use 4 (default) for median values of Tierpsy features; use 5 for mean
else
    metricRange = [1:9];
end

if applyBonferroniCorrection
    mappingStrainNumThreshold = 20; % threshold only used for Tierpsy features - num of strains showing density-dependent difference 
else
    mappingStrainNumThreshold = 40;
end

%% prep work
addpath('auxiliary/')
% load the strain names included in the specified strainSet
load(['strainsList/' strainSet '.mat'])
% get list of file names for each strain
[strainFileList,~] = getFileList(strains);

%% calculate featureMetric values
%% initialise
if useTierpsyFeat
    % load feature values
    load(['results/TierpsyFeat_' feature '_all.mat'])
    sw40_feature = allFeats;
    clear('allFeats')
    load(['results/5worm_TierpsyFeat_' feature '_all.mat'])
    sw5_feature = allFeats;
    clear('allFeats')
    % get feature list
    featList = sw40_feature.(strains{1})(:,1);
else
    if strcmp(feature,'speed_nonTierpsy')
        % load feature values
        load('results/speed_all_smooth3s.mat')
        sw40_feature = sw_feature;
        clear('sw_feature','mw_feature','cluster_feature')
        load('results/5worm_speed_all_smooth3s.mat')
        sw5_feature = sw_feature;
        clear('sw_feature','mw_feature','cluster_feature')
    end
    % get feature list
    featList = {feature};
end
numFeats = numel(featList);

%% generate feature matrix csv for mapping all 40 and 5 worm Tierpsy_256 data median values
if strcmp(feature,'Tierpsy_256')
    Tierpsy_256_40worms = cell(numel(strains),257);
    Tierpsy_256_5worms = cell(numel(strains),257);
    Tierpsy_256_40worms(:,1) = strains;
    Tierpsy_256_5worms(:,1) = strains;
    for strainCtr = 1:numel(strains)
        strain = strains{strainCtr};
        strainFeats_40worms = median(cell2mat(sw40_feature.(strain)(:,2:end)),2)'; format long g
        strainFeats_5worms = median(cell2mat(sw5_feature.(strain)(:,2:end)),2)'; format long g
        Tierpsy_256_40worms(strainCtr,2:end) = num2cell(strainFeats_40worms);
        Tierpsy_256_5worms(strainCtr,2:end) = num2cell(strainFeats_5worms);
    end
    % add feature name title
    featnames = horzcat({'strain'}, sw40_feature.(strain)(:,1)');
    Tierpsy_256_40worms = vertcat(featnames,Tierpsy_256_40worms);
    Tierpsy_256_5worms = vertcat(featnames,Tierpsy_256_5worms);
    if saveResults
        saveDir = '/Users/sding/Documents/AggScreening/results/mapping/swDensityCompare/Tierpsy/Tierpsy_256/';
        % for some reason cannot direct export via cell2csv otherwise feat names truncated, must go through .txt as intermediate step
        % first export as .txt
        dir40 = [saveDir 'cegwasMappingFeatures40worms.txt'];
        dir5 = [saveDir 'cegwasMappingFeatures5worms.txt'];
        dlmcell(dir40,Tierpsy_256_40worms);
        dlmcell(dir5,Tierpsy_256_5worms);
        % then read in .txt
        fid40 = fopen(dir40);
        fid5 = fopen(dir5);
        rawContent40 = textscan(fid40,'%s'); % grab content
        rawContent5 = textscan(fid5,'%s'); % grab content
        formatContent40 = rawContent40{1};
        formatContent5 = rawContent5{1};
        formatContent40 = reshape(formatContent40,[257,numel(formatContent40)/257])';
        formatContent5 = reshape(formatContent5,[257,numel(formatContent5)/257])';
        % then export as csv
        dir40 = strrep(dir40,'.txt','.csv');
        dir5 = strrep(dir5,'.txt','.csv');
        cell2csv(dir40,formatContent40);
        cell2csv(dir5,formatContent5);
    end
end

%% go through each feature
for featCtr = 1:numFeats
    % through each strain
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        num40Files = size(sw40_feature.(strain),2); % -1 because first column contains feature name
        num5Files = size(sw5_feature.(strain),2);
        if useTierpsyFeat
            num40Files = num40Files-1; % -1 because first column contains feature name
            num5Files = num5Files-1;
        end
        % initialise 3D matrix to hold feature metric values
        featVals_40.(strain) = NaN(num40Files,9);
        featVals_5.(strain) = NaN(num5Files,9);
        % populate matrix with variables
        % (different metrics in columns will look the same for Tierpsy features as only plate average is available)
        for repCtr = 1:num40Files
            if useTierpsyFeat
                featThisRep_40 = sw40_feature.(strain){featCtr,repCtr+1}; % +1 because first column contains feature name
            else
                featThisRep_40 = sw40_feature.(strain){featCtr,repCtr};
            end
            if ~isempty(featThisRep_40)
                featVals_40.(strain)(repCtr,1) = min(featThisRep_40);
                featVals_40.(strain)(repCtr,2) = prctile(featThisRep_40,10);
                featVals_40.(strain)(repCtr,3) = prctile(featThisRep_40,25);
                featVals_40.(strain)(repCtr,4) = nanmedian(featThisRep_40);
                featVals_40.(strain)(repCtr,5) = nanmean(featThisRep_40);
                featVals_40.(strain)(repCtr,6) = prctile(featThisRep_40,75);
                featVals_40.(strain)(repCtr,7) = prctile(featThisRep_40,90);
                featVals_40.(strain)(repCtr,8) = max(featThisRep_40);
                featVals_40.(strain)(repCtr,9) = iqr(featThisRep_40);
            end
        end
        for repCtr = 1:num5Files
            if useTierpsyFeat
                featThisRep_5 = sw5_feature.(strain){featCtr,repCtr+1};
            else
                featThisRep_5 = sw5_feature.(strain){featCtr,repCtr};
            end
            if ~isempty(featThisRep_5)
                featVals_5.(strain)(repCtr,1) = min(featThisRep_5);
                featVals_5.(strain)(repCtr,2) = prctile(featThisRep_5,10);
                featVals_5.(strain)(repCtr,3) = prctile(featThisRep_5,25);
                featVals_5.(strain)(repCtr,4) = nanmedian(featThisRep_5);
                featVals_5.(strain)(repCtr,5) = nanmean(featThisRep_5);
                featVals_5.(strain)(repCtr,6) = prctile(featThisRep_5,75);
                featVals_5.(strain)(repCtr,7) = prctile(featThisRep_5,90);
                featVals_5.(strain)(repCtr,8) = max(featThisRep_5);
                featVals_5.(strain)(repCtr,9) = iqr(featThisRep_5);
            end
        end
    end
    %% save variable
    if ~useTierpsyFeat
        if saveResults
            save(strcat('results/featMetricVals/swDensityCompare_',featList{featCtr},'_all.mat'),'featVals_40','featVals_5');
        end
    end
   
    %% compare featVals between 40 and 5 worm recordings
    % initialise
    pVals = cell(length(strains),10);
    medianDiff = cell(length(strains),10);
    significantStrainCtr = 0; % counter only useful for Tierpsy features - to be added to saved output file name to help distinguish useful features
    % go through each strain, conduct t test for each metric between 40 and 5 worm datasets
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        pVals{strainCtr,1} = strain;
        medianDiff{strainCtr,1} = strain;
        for metricCtr = metricRange
            [~,p]=ttest2(featVals_40.(strain)(:,metricCtr),featVals_5.(strain)(:,metricCtr));
            if applyBonferroniCorrection
                if p<0.05/9 % correct for multiple comparison
                    pVals{strainCtr,1+metricCtr} = p;
                    significantStrainCtr = significantStrainCtr+1;
                else
                    pVals{strainCtr,1+metricCtr} = 'ns';
                end
            else
                if p<0.05
                    pVals{strainCtr,1+metricCtr} = p;
                    significantStrainCtr = significantStrainCtr+1;
                else
                    pVals{strainCtr,1+metricCtr} = 'ns';
                end
            end
            medianDiff{strainCtr,1+metricCtr} = nanmedian(featVals_40.(strain)(:,metricCtr))-nanmedian((featVals_5.(strain)(:,metricCtr)));
        end
    end
    % add heading
    headingText = {'strain', 'min','10prc','25prc','median','mean','75prc','90prc','max','iqr'};
    pVals = vertcat(headingText,pVals);
    medianDiff = vertcat(headingText,medianDiff);
    % trim down for Tierpsy features (keep only the median)
    if useTierpsyFeat
        assert(numel(metricRange)==1)
        pVals = pVals(:,[1,metricRange+1]);
        medianDiff = medianDiff(:,[1,metricRange+1]);
    end
    % save and export
    if saveResults
        dirName = 'results/swDensityCompare/';
        featName = featList{featCtr};
        if useTierpsyFeat
            dirName = [dirName 'Tierpsy/' feature '/'];
            featName = char(featName);
        end
        saveFileName = strcat(dirName, featName);
        if applyBonferroniCorrection
            saveFileName = strcat(saveFileName, '_corrected');
        end
        if useTierpsyFeat
            saveFileName = strcat(saveFileName, '_sigStr', num2str(significantStrainCtr));
        end
        save(strcat(saveFileName, '.mat'),'pVals');
        % 
        if useTierpsyFeat
            featpValExportName = strcat('results/mapping/swDensityCompare/Tierpsy/', feature, '/', featName, '_medianDiff.txt');
            % remove DA609 from mapping
            removeIdx = find(strcmp(medianDiff(:,1),'DA609')); % should be 18
            medianDiff = vertcat(medianDiff(1:removeIdx-1,:),medianDiff(removeIdx+1:end,:));
            % export
            if applyBonferroniCorrection & significantStrainCtr >= mappingStrainNumThreshold 
                dlmcell(featpValExportName,medianDiff);
            elseif ~applyBonferroniCorrection & significantStrainCtr >= mappingStrainNumThreshold
                dlmcell(featpValExportName,medianDiff);
            end
        else
            % remove DA609 from mapping
            removeIdx = find(strcmp(medianDiff(:,1),'DA609')); % should be 18
            medianDiff = vertcat(medianDiff(1:removeIdx-1,:),medianDiff(removeIdx+1:end,:));
            % export
            featpValExportName = strcat('results/mapping/swDensityCompare/',featName, '_medianDiff.txt');
            dlmcell(featpValExportName,medianDiff);
        end
    end
end