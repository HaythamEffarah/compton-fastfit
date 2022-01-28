function data4D = get4Dfromcoeff(coeffData)
%GET4DFROMCOEFF Takes array of polynomial coefficients and returns 4D data.
%

for i = 1:m
    for j = 1:n
        coeffData(i,j,tmp,2) = polyval(squeeze(spectrumPolyCoefficients(i, j,:)), truncatedEnergies, [], squeeze(meanErrorMatrix(i,j,:)));
    end
end

