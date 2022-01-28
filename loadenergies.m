function [energies] = loadenergies(file)
%LOADENERGIES Load data from energy/x-angle/y-angle/counts text file.
%
%   S = loadenergies(file) returns a 4-D double of dimensions (xAngleBins,
%   yAngleBins, energyBins, (energy, counts)) where (:,:,:,1) denotes the
%   energy at a given coordinate, and (:,:,:,2) denotes the number of
%   photons at a given coordinate.
% 
%   Note: This function is designed to process output data from the 
%   isv10T.nb code, and is valid for any arbitrary number of x-angle, 
%   y-angle, and energy bins as long as the coordinate nesting is as 
%   follows in the input file:
%   ENERGY > X-ANGLE > Y-ANGLE
%   
%   Future plans:
%       - Allow input of 1-quadrant of code & extrapolate other quadrants.
%       - Allow robustness of coordinate nesting.


raw = load(file);

indexAngleCombinations = size(find(raw(:,1) == raw(1,1)),1);
energyBins = size(raw,1) / indexAngleCombinations;
yAngleBins = size(find(raw(:,2) == (raw(1,2)) & raw(:,1) == (raw(1,1))),1);
xAngleBins = indexAngleCombinations / yAngleBins;

reshapedRaw = reshape(raw, xAngleBins, yAngleBins, energyBins, 4);
% raw = reshape(raw, 101, 101, [], 4);
shape = size(reshapedRaw);

energies = zeros(shape(1), shape(2), shape(3), 2);
energies(:,:,:,1) = reshapedRaw(:,:,:,1);
energies(:,:,:,2) = reshapedRaw(:,:,:,4);
