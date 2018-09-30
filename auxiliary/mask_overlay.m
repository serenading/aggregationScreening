function [ hf ] = mask_overlay( I, BW, RGB, alpha )
%mask_overlay Plots a mask overlaid to a greyscale image
%   I has to be greyscale, BW has to be logical, colour has to be a RGB
%   triplet matlab-style and alpha a double between 0 and 1


if max(RGB)>1 || numel(RGB)~=3 || ~isvector(RGB)
    error('RGB has to be a rgb vector matlab-style')
end

if size(I)~=size(BW)
    error('I and BW have to be of same size');
else 
    sz = size(I);
end

alpha = min(alpha,1); %%this to check that alpha is maximum 1

RGB(1,1,1:3) = RGB(:);

mask_RGB = BW(:,:,[1,1,1]).* RGB(ones(sz(1),1),ones(sz(2),1),:);

hf = figure;
imshow(I,[]);
hold on
image(mask_RGB,'alphadata',alpha.*BW)


end

