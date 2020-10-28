%[sweptStrains, notSweptStrains] = getSweptStrains(method)

%% Script makes a list of swept vs. non-swept strains from the imaging dataset in order to e.g. run classification separately on those.
% method: 'conservative' or 'liberal'. 'conservative' classes a strain as
% swept if all four chromosomes (I, IV, V or X) are swept (>=30% of the
% same haplotype). 'liberal' classes a strain as swept if any of the four
% chromosomes are swept.
% Author: @serenading. Oct 2020.

%% Set parameter
method = 'conservative'; % 'conservative' or 'liberal'

%% Load haplotype data
fpath = '/Users/sding/OneDrive - Imperial College London/aggScreening/source/AndersenLab/sweep_summary_20200815.tsv';
hapshares = tdfread(fpath);

%% Get swept strains master list
if strcmp(method,'conservative')
    sweptMasterLogInd = hapshares.swept_chroms==4;
elseif strcmp(method,'liberal')
    sweptMasterLogInd = hapshares.swept_chroms>=1;
else
    warning('Please specify method as conservative or liberal')
end
% find all swept strains from the Andersen lab list
sweptStrains_master = cellstr(hapshares.isotype(sweptMasterLogInd,:));

%% Get list of swept vs. non-swept strains that were imaged in the dataset
load('strainsList/all.mat')
sweptLogInd = false(numel(strains),1); % to keep track of swept strains inside the sweptMaster list
strainLogInd = false(numel(strains),1); % to keep track of strains inside the master list (for example DA609 will not be found in this list)
for strainCtr = 1:numel(strains)
    strain = strains(strainCtr);
    if ismember(strain,sweptStrains_master)
        sweptLogInd(strainCtr) = true;
    end
    if ismember(strain,hapshares.isotype)
        strainLogInd(strainCtr) = true;
    end
end
sweptStrains = strains(sweptLogInd);
nonSweptStrains = strains(~sweptLogInd & strainLogInd);

%% Save
save(['strainsList/swept_strains_' method '.mat'],'sweptStrains','nonSweptStrains')