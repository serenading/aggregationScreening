function [ corr_coeff ] = image_corr( a, b )
%same as xcorr2 but specifying valid when calling conv2

if nargin == 1
	b = a;
end

c = conv2(a, rot90(conj(b),2), 'valid');
end