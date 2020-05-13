%% This script uses features extracted from a control dataset to determine
% which of the Tierpsy features have normal distribution.
% The control dataset to use are control strains (N2, CB4856, DA609) from
% the 2018 aggregation screening, at 5 worm density.
% Use both Kolmogorov-Smirnov and Shapiro-Wilk tests, as n is just over 50
% for each strain.
% author: serenading. May 2020

clear
close all

addpath('auxiliary')

%% Import features matrices and load strains list
featExtractTimestamp = '20191024_122847'; %'20191024_122847' or '20181203_141111'
featureTable = readtable(['/Users/sding/Dropbox/aggScreening/results/fiveWorm/fiveWormFeaturesTable_' featExtractTimestamp '.csv'],'Delimiter',',','preserveVariableNames',true);
n_nonFeatVar = 17; % the first n columns of the feature table that do not contain features
n_feats = size(featureTable,2)-n_nonFeatVar;
featNames = featureTable.Properties.VariableNames(n_nonFeatVar+1:end);

%% Get control strain logical indices
N2LogInd = strcmp(featureTable.strain_name,'N2');
CB4856LogInd = strcmp(featureTable.strain_name,'CB4856');
DA609LogInd = strcmp(featureTable.strain_name,'DA609');
controlLogInd = N2LogInd | CB4856LogInd | DA609LogInd;

    
%% Go through each feature, test for normal distribution
if ~exist('/Users/sding/Dropbox/aggScreening/results/featuresDistribution/normalitytest_pvalues.mat')

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
    save('/Users/sding/Dropbox/aggScreening/results/featuresDistribution/normalitytest_pvalues.mat','KSTestP','SWTestP','featNames');
else
    load('/Users/sding/Dropbox/aggScreening/results/featuresDistribution/normalitytest_pvalues.mat','KSTestP','SWTestP', 'featNames');
end

%% Look into which features have normal distribution (p>0.05)
% Check how many tests show normal distribution using either method
KSNormalLogInd = KSTestP > 0.05; % gives 15 nnz
SWNormalLogInd = SWTestP > 0.05; % gives 7338 nnz
% Check how many features have normal distribution for which strain set
SWNormalLogInd_sum1 = sum(SWNormalLogInd,1); % shows about half of the features have normal distribution for one control strain
SWNormalLogInd_sum2 = sum(SWNormalLogInd,2); % shows for how many control strain sets the feature is normally distributed
% Keep record of which features have normal distribution within at least one control strain set
SWNormalLogInd_sum2LogInd = SWNormalLogInd_sum2>=1;
normalFeatNames = featNames(SWNormalLogInd_sum2LogInd)';
nonNormalFeatNames = featNames(~SWNormalLogInd_sum2LogInd)';
% save
save('/Users/sding/Dropbox/aggScreening/results/featuresDistribution/whichFeatNormal.mat','normalFeatNames','nonNormalFeatNames');