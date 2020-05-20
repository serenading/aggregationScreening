%% Function shortens Tierpsy feature names that are above 63 characters, the apparently the limit for Matlab

% Input & Output: cell array containing Tierpsy feature names

function featureNames = shortenFeatNames(featureNames)

for featCtr = 1:numel(featureNames)
    featureName = featureNames{featCtr};
    
    % replace long names with short
    featureName = strrep(featureName,'velocity','vel');
    featureName = strrep(featureName,'relative','rel');
    featureName = strrep(featureName, 'angular', 'ang');
    
%     % further shorten if necessary (normally this will not be necessary)
%     featNameLength = numel(featureName);
%     if featNameLength > 63
%         featureName = featureName(1:63);
%     end
    featureNames{featCtr} = featureName;
end

end