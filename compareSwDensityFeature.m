clear
close all

strainSet = 'all'; % 'controls','divergent','all'
feature = 'food'; % 'speed','food'
applyBonferroniCorrection = false;
saveResults = true;

if strcmp(feature, 'speed')
    useTierpsyFeat = false;
elseif strcmp(feature, 'food')
    useTierpsyFeat = true;
end

if useTierpsyFeat
    metricRange = 4;
else
    metricRange = [1:9];
end

%% prep work
addpath('auxiliary/')
% load the strain names included in the specified strainSet
load(['strainsList/' strainSet '.mat'])
% get list of file names for each strain
[strainFileList,~,~] = getFileList(strains);

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
    if strcmp(feature,'speed')
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

%% go through each feature
for featCtr = 1:numFeats
    % through each strain
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        if useTierpsyFeat
            sw40_feature.(strain) = (sw40_feature.(strain)(:,2:end));  % remove first column containing feature names
            sw5_feature.(strain) = (sw5_feature.(strain)(:,2:end));
        end
        num40Files = size(sw40_feature.(strain),2);
        num5Files = size(sw5_feature.(strain),2);
        % initialise 3D matrix to hold feature metric values
        featVals_40.(strain) = NaN(num40Files,9);
        featVals_5.(strain) = NaN(num5Files,9);
        % populate matrix with variables
        % (different metrics in columns will look the same for Tierpsy features as only plate average is available)
        for repCtr = 1:num40Files
            featThisRep_40 = sw40_feature.(strain){featCtr,repCtr};
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
            featThisRep_5 = sw5_feature.(strain){featCtr,repCtr};
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
    % trim down for Tierpsy features
    if useTierpsyFeat
        assert(numel(metricRange)==1)
        pVals = pVals(:,[1,metricRange+1]);
        medianDiff = medianDiff(:,[1,metricRange+1]);
    end
    % save and export
    if saveResults
        dirName = 'results/swDensityCompare/';
        if useTierpsyFeat
            dirName = [dirName 'Tierpsy/'];
        end
        saveFileName = strcat(dirName, featList{featCtr}, '_all');
        if applyBonferroniCorrection
            saveFileName = strcat(saveFileName, '_corrected');
        end
        if useTierpsyFeat
            saveFileName = strcat(saveFileName, '_sigStr', num2str(significantStrainCtr));
        end
        save(strcat(saveFileName, '.mat'),'pVals');
        % 
        if useTierpsyFeat
            featpValExportName = strcat('results/mapping/swDensityCompare/Tierpsy/', featList{featCtr}, '_medianDiff_all.txt');
            if applyBonferroniCorrection & significantStrainCtr >= 5 
                dlmcell(featpValExportName,medianDiff);
            elseif ~applyBonferroniCorrection & significantStrainCtr >= 10
                dlmcell(featpValExportName,medianDiff);
            end
        else
            featpValExportName = strcat('results/mapping/swDensityCompare/', featList{featCtr}, '_medianDiff_all.txt');
            dlmcell(featpValExportName,medianDiff);
        end
    end
end