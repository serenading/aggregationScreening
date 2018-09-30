function [ ] = superimpose_a_mask( h, mask, RGB, alpha )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if max(RGB)>1 || numel(RGB)~=3 || ~isvector(RGB)
    error('RGB has to be a rgb vector matlab-style')
end

alpha = min(alpha,1); %%this to check that alpha is maximum 1

sz = size(mask);

RGB(1,1,1:3) = RGB(:);

mask_RGB = mask(:,:,[1,1,1]).* RGB(ones(sz(1),1),ones(sz(2),1),:);

if isgraphics(h,'figure')
    ha = findall(h,'type','axes'); %find axes in figure handle
    ha = ha(1); %only one is needed
elseif isgraphics(h,'axes')
    ha = h;
else
    error('handle required');
end

hold(ha, 'on');
image(mask_RGB,'alphadata',alpha.*mask);
hold(ha,'off');

end

