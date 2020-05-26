%% This script uses features extracted from a control dataset to determine which of the Tierpsy features have normal distribution.
% The control dataset to use are control strains (N2, CB4856, DA609) from the 2018 aggregation screening, at 5 worm density.
% Use both Kolmogorov-Smirnov and Shapiro-Wilk tests, as n is just over 50 for each strain.
% author: serenading. May 2020

clear
close all

addpath('auxiliary')

%% Things to consider:
% 1. effect size for determining what's non-Gaussian

%% Set analysis parameters
bonCorr = true;
whichTest = 'SW'; % 'KS' or 'SW' % 'SW' seems to work better, otherwise most features are not Gaussian when correcting for multiple comparison

%% Import features matrices and load strains list
% set which feature extraction timestamp to use
featExtractTimestamp = '20191024_122847';
extractStamp = featExtractTimestamp;

featureTable = readtable(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/fiveWorm/fiveWormFeaturesTable_' extractStamp '.csv'],'Delimiter',',','preserveVariableNames',true);
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features
n_feats = size(featureTable,2)-n_nonFeatVar;
featNames = featureTable.Properties.VariableNames(n_nonFeatVar+1:end);

%% Get control strain logical indices
N2LogInd = strcmp(featureTable.strain_name,'N2');
CB4856LogInd = strcmp(featureTable.strain_name,'CB4856');
DA609LogInd = strcmp(featureTable.strain_name,'DA609');
controlLogInd = N2LogInd | CB4856LogInd | DA609LogInd;
    
%% Go through each feature, test for normal distribution
testValFilename = ['/Users/sding/OneDrive - Imperial College London/aggScreening/results/featuresDistribution/normalitytest_pvalues_' extractStamp '.mat'];
if ~exist(testValFilename)

    % Initialise matrix to hold p values
    KSTestP = NaN(n_feats,4);
    SWTestP = NaN(n_feats,4);
    % convert feature table to matrix
    featureMat = table2array(featureTable(:,n_nonFeatVar+1:end));
    
    % Go through each feature to get KS test p-values on each control strain set
    for featCtr = 1:n_feats
        % Conduct one-sample Kolmogorov-Smirnov test on each of the control strains separately, record p value
        [~,KSTestP(featCtr,1)] = kstest(featureMat(N2LogInd,featCtr));
        [~,KSTestP(featCtr,2)] = kstest(featureMat(CB4856LogInd,featCtr));
        [~,KSTestP(featCtr,3)] = kstest(featureMat(DA609LogInd,featCtr));
        [~,KSTestP(featCtr,4)] = kstest(featureMat(controlLogInd,featCtr));
        % Conduct Shapiro-Wilk test on each of the control strains separately, record p value
        if std(featureMat(controlLogInd,featCtr))~=0 % SW test does not work when standard deviation is zero
            [~,SWTestP(featCtr,1)] = swtest(featureMat(N2LogInd,featCtr));
            [~,SWTestP(featCtr,2)] = swtest(featureMat(CB4856LogInd,featCtr));
            [~,SWTestP(featCtr,3)] = swtest(featureMat(DA609LogInd,featCtr));
            [~,SWTestP(featCtr,4)] = swtest(featureMat(controlLogInd,featCtr));
        end
        % Display progress
        disp([ num2str(featCtr/n_feats*100) '% of all features tested'])
    end
    
    %% Save p value table
    save(testValFilename,'SWTestP','featNames');
else
    load(testValFilename,'KSTestP','SWTestP', 'featNames');
end

%% Look into which features have normal distribution
% Apply Bonferroni correction for multiple comparisons if specified
if bonCorr
    pThreshold = 0.05/n_feats;
else
    pThreshold = 0.05;
end
% Check how many tests show normal distribution using either method
% set non-significant p-values to NaN
KSNormalLogInd = KSTestP > pThreshold; % gives 15 nnz /100 MC
SWNormalLogInd = SWTestP > pThreshold; % gives 7338 nnz /13294 MC
if strcmp(whichTest,'KS')
    % Check how many features have normal distribution for which strain set
    KSNormalLogInd_sum1 = sum(KSNormalLogInd,1); % shows about half of the features have normal distribution for one control strain
    KSNormalLogInd_sum2 = sum(KSNormalLogInd,2); % shows for how many control strain sets the feature is normally distributed
    % Keep record of which features have normal distribution within at least one control strain set
    KSNormalLogInd_sum2LogInd = KSNormalLogInd_sum2>=1;
    normalFeatNames = featNames(KSNormalLogInd_sum2LogInd)';
    nonNormalFeatNames = featNames(~KSNormalLogInd_sum2LogInd)';
elseif strcmp(whichTest,'SW')
    % Check how many features have normal distribution for which strain set
    SWNormalLogInd_sum1 = sum(SWNormalLogInd,1); % shows about half of the features have normal distribution for one control strain
    SWNormalLogInd_sum2 = sum(SWNormalLogInd,2); % shows for how many control strain sets the feature is normally distributed
    % Keep record of which features have normal distribution within at least one control strain set
    SWNormalLogInd_sum2LogInd = SWNormalLogInd_sum2>=1;
    normalFeatNames = featNames(SWNormalLogInd_sum2LogInd)';
    nonNormalFeatNames = featNames(~SWNormalLogInd_sum2LogInd)';
else
    error('Must specify whichTest as either KS or SW')
end
% save
save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/featuresDistribution/whichFeatNormal' whichTest 'test_' extractStamp '.mat'],'normalFeatNames','nonNormalFeatNames');
% generate short version of feat names for older versions of the features
if strcmp(featExtractTimestamp,'20191024_122847')
    nonNormalFeatNames = cellfun(@(x) strrep(strrep(strrep(x,'angular','ang'),'relative','rel'),'velocity','vel'), nonNormalFeatNames, 'UniformOutput', false);
    normalFeatNames = cellfun(@(x) strrep(strrep(strrep(x,'angular','ang'),'relative','rel'),'velocity','vel'), normalFeatNames, 'UniformOutput', false);
    save(['/Users/sding/OneDrive - Imperial College London/aggScreening/results/featuresDistribution/whichFeatNormal' whichTest 'test_short_' extractStamp '.mat'],'normalFeatNames','nonNormalFeatNames');
end

%% Using SW test correcting for multiple comparison, 103 out of 4548 features are non-Gaussian. 
% Most of these are path curvature features, some are blob features.