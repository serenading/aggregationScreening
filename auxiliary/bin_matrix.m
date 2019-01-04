function [ binnedIM, binning_map, binnedIMstd ] = bin_matrix( IM, bsz, flag_std, hfun )
%UNTITLED Summary of this function goes here
%   quite slow using accumarray though :(

%% input check

if nargin < 3 || isempty(flag_std)
    flag_std = false;
end

% allow user to ask for a function to accumarray
if nargin < 4 || isempty(hfun)
    flag_onlysum = true;
else
    flag_onlysum = false;
end

if numel(bsz) == 1
    bsz = bsz(1,[1 1]);
end

%% useful numbers

sz = size(IM);

if mod(sz(1),bsz(1)) || mod(sz(2),bsz(2))
    IM = IM(1:end-mod(sz(1),bsz(1)),1:end-mod(sz(2),bsz(2)));   %this prevents errors later in case one try to bin with a boxsize that doesn't fit in the image an integer number of times
end
sz = size(IM);  %maybe got chopped by line before
newsz = floor(sz./bsz); %should now be of integers, but never too sure


%% creation of binning_map

id = 1:prod(newsz);   %creates indices for pixels in the input image
id = reshape(id,newsz);         %makes it a row


ii = ceil((1:bsz(1)*newsz(1))/bsz(1)); %creates blowed up indices
jj = ceil((1:bsz(2)*newsz(2))/bsz(2));

binning_map = id(ii,jj);

%% actual averaging over the boxsize

if flag_onlysum % just do average
    
    norm = 1/prod(bsz);
    binnedIM = accumarray(binning_map(:),IM(:)) * norm ; %average
    binnedIM = reshape(binnedIM, newsz);
    
else
    binnedIM = accumarray(binning_map(:),IM(:),[],hfun) ; 
    binnedIM = reshape(binnedIM, newsz);
end

if flag_std && nargout == 3 
    binnedIMstd = accumarray(binning_map(:),IM(:),[],@std) ; %average
    binnedIMstd = reshape(binnedIMstd, newsz);

end

