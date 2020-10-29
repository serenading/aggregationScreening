function featureTable = balanceStrainReps(featureTable, featExtractTimestamp, crossVal_k, strains)

%% Function subsamples replicates to ensure that all strains are balanced in the number of replicates represented in a classification task. 
% Because there are more replicates of the divergent panel and of DA609 and
% this will necessarily impact classification accuracy when they are
% overrepresented. 

% Author: @serenading. Oct 2020.

dropInd = [];
if strcmp(featExtractTimestamp,'20201021_114135') | strcmp(featExtractTimestamp,'20200630_171156')
    n_reps2keep = 3*5; % 3 actual replicates x 5 fold augmentation for augmented data
else
    n_reps2keep = crossVal_k+1;
end

for strainCtr = 1:numel(strains)
    allInd = find(strcmp(strains{strainCtr},featureTable.strain_name));
    if numel(allInd)>n_reps2keep
        % Ensure that when augmented trajectories are used, whole replicates are subsampled together
        if strcmp(featExtractTimestamp,'20201021_114135') | strcmp(featExtractTimestamp,'20200630_171156')
            firstInd = allInd(1:5:numel(allInd)); % get the first indices of each of the augmented set of five
            dropInd = vertcat(dropInd,datasample(firstInd,numel(firstInd)-3,'Replace',false));
        else
            dropInd = vertcat(dropInd,datasample(allInd,numel(allInd)-n_reps2keep,'Replace',false));
        end
    end
end

% Ensure that when augmented trajectories are used, dropInd are expanded to include all five augmented rows
if strcmp(featExtractTimestamp,'20201021_114135') | strcmp(featExtractTimestamp,'20200630_171156')
    dropInd = [dropInd, dropInd+1, dropInd+2, dropInd+3, dropInd+4];
    dropInd = sort(dropInd(:));
end

% Drop the extra replicates
dropRowLogInd = false(1,numel(featureTable.strain_name));
dropRowLogInd(dropInd) = true;
featureTable = featureTable(~dropRowLogInd,:);
disp(['Random sampling applied so that each of the ' num2str(numel(strains)) ' strains has ' num2str(n_reps2keep) ' replicates for classification.'])