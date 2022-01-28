function [spectrum2D] = get2Dfrom4D(spectrum4D)
%GET2DFROM4D Converts 4D data to 2D intensity map.
%   Takes a 4D double (x-angle, y-angle, energy bins, 2) and
%   returns a 2D double (x-angle, y-angle) with each element being the
%   integrated intensity at each pixel.
%
%   
%   Future plans:
%
xdim = size(spectrum4D, 1);
ydim = size(spectrum4D, 2);
% spectrumSummed = squeeze(sum(spectrum4D, 3));
% spectrum2D = reshape(spectrumSummed(:, :, 2), xdim, ydim).';

spectrum2D = zeros(xdim, ydim);

for i = 1:ydim
    for j = 1:xdim
        spectrum2D(j,i) = trapz(squeeze(spectrum4D(i,j,:,1)), squeeze(spectrum4D(i,j,:,2)));
    end
end

