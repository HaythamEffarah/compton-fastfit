function [spectrumPolyCoefficients, eBounds, meanErrorMatrix, spectrumFitted] = polynomialspectrumfit(data4D, polyOrder)
%POLYNOMIALSPECTRUMFIT Fits polynomials pixel-by-pixel to x-ray distribution.
%
%   polynomialspectrumfit.m
%   Takes an LCS energy vs. x-angle vs. y-angle vs. counts array and
%   returns a pixel-by-pixel polynomial fit of order n.
%
%   Use:
%
%
%   Inputs:
%       data4D: data imported from CSV file from direct simulation code
%       output with ENERGY > X-ANGLE > Y-ANGLE nesting
%
%       polyOrder: int of polynomial order we seek to fit
%
%   Required Codes:
%       loadenergies.m
%
%   LINES FOR USE AS A SCRIPT
%   polyOrder = 4;
%   path = "/Users/haytham/Documents/GitHub/lcs-fxn-fit/";
%   file = "/xraydist/35keV_6um_xyn.csv";
%   spectrumData = loadenergies(path + file);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% tic
% fprintf('\nLoading spectrum ' + file + '.\n');
% fprintf('This step will take approximately 35 seconds.\n');
% spectrumData = loadenergies(file);
% toc


tic
% fprintf('Starting fit.\n');

% --- Begin Trevor edit
% This is hardcoded in!!
% You may want to make this a function later for changing if needed
% Or maybe set it as some percentage region around the middle so it is not
%   dependent on the size of the input
% Must define the "peaky_zone" because we only want to fit the peak in this zone
peaky_zone = [40:60];
options = optimoptions('fmincon', 'Display', 'off', ...
    'FunctionTolerance', 1e-10, 'StepTolerance', 1e-12);
% --- End Trevor edit

[m, n, o, ~] = size(data4D);
spectrumPolyCoefficients = zeros(m, n, (polyOrder + 1));
meanErrorMatrix = zeros(m, n, 2);
spectrumFitted = zeros(m, n, o, 2);
scalefactor = zeros(m, n);
eBounds = zeros(m,n,2);
spectrumFitted(:,:,:,1) = data4D(:,:,:,1);

for i = 1:m
    for j = 1:n
        energy = squeeze(data4D(i,j,:,1));
        dEnergy = energy(2) - energy(1);
        intensity = squeeze(data4D(i,j,:,2));
        power = trapz(energy,intensity);

        nonZeroLocations = find(intensity);

        truncatedIntensities = intensity(nonZeroLocations);
        truncatedEnergies = energy(nonZeroLocations);

        %Produce polynomial coefficients to fit to the data at pixel (i,j)
        %given an x-rescaling function mu to produce the correct x-values
        [p,~,mu] = polyfit(truncatedEnergies, truncatedIntensities, polyOrder);

        % --- Begin Trevor edit
        if ismember(i, peaky_zone) && ismember(j, peaky_zone)
            tmp = (truncatedEnergies - mu(1)) / mu(2);
            fun = @(coeff) objective(tmp, truncatedIntensities, coeff);
            nonlcon = @(coeff) nonlin_constraints(tmp, truncatedIntensities, coeff);
            p = fmincon(fun, p, [], [], [], [], [], [], nonlcon, options);
        end

        %Correct range of energies (x) to include only non-negative values in
        %the evaluated polynomial distribution
        y = polyval(p,truncatedEnergies,[],mu);

        spectrumPolyCoefficients(i, j,:) = p;
        meanErrorMatrix(i,j,:) = mu;

        spectrumFitted(i,j,nonZeroLocations,2) = polyval(squeeze(spectrumPolyCoefficients(i, j,:)), energy(nonZeroLocations), [], squeeze(meanErrorMatrix(i,j,:)));

        %Clean up stray values outside of single mode distribution
        [~, maxLoc] = max(spectrumFitted(i,j,:,2));
        scanUp = maxLoc;
        scanDown = maxLoc;
        while spectrumFitted(i,j,scanUp,2) >= (0.5 * spectrumFitted(i,j,maxLoc,2))
            scanUp = scanUp + 1;
        end

        while (spectrumFitted(i,j,scanUp,2) - spectrumFitted(i,j,scanUp+1,2)) > 0 ...
                && spectrumFitted(i,j,scanUp,2) > 0
            scanUp = scanUp + 1;
        end
        spectrumFitted(i,j,scanUp:end,2) = 0;


        while spectrumFitted(i,j,scanDown,2) >= (0.5 * spectrumFitted(i,j,maxLoc,2))
            scanDown = scanDown - 1;
        end

        while (spectrumFitted(i,j,scanDown,2) - spectrumFitted(i,j,scanDown-1,2)) > 0 ...
                && spectrumFitted(i,j,scanDown,2) > 0
            scanDown = scanDown - 1;
        end
        spectrumFitted(i,j,1:scanDown,2) = 0;

        % Now redefine nonZeroLocations to include those locations after
        % cleaning up the polynomial
        nonZeroLocations = find(spectrumFitted(i,j,:,2)>0);
        start = nonZeroLocations(1);
        finish = nonZeroLocations(end);

        fittedPower = trapz(squeeze(spectrumFitted(i,j,:,1)), squeeze(spectrumFitted(i,j,:,2)));
        scalefactor = power / fittedPower;

        % --- Begin Trevor edit
        % No longer want to renormalize since it fits to the integral
        % I commented out the line below
        spectrumPolyCoefficients(i, j,:) = scalefactor * spectrumPolyCoefficients(i, j,:);
        % --- End Trevor edit

        spectrumFitted(i,j,nonZeroLocations,2) = polyval(squeeze(spectrumPolyCoefficients(i, j,:)), energy(nonZeroLocations), [], squeeze(meanErrorMatrix(i,j,:)));

        eBounds(i,j,1) = energy(start) - dEnergy;
        eBounds(i,j,2) = energy(finish) + dEnergy;
    end
end
toc

% --- Begin Trevor edit
% Functions to call for fitting
function [y] = model(x, coeff)
  y = polyval(coeff, x);
end

function [c, ceq] = nonlin_constraints(x, ydata, coeff)
  c = [];
  yfit = model(x, coeff);
  trap = trapz(x, yfit) - trapz(x, ydata);

  [~, argmax] = max(ydata);
  ydif = diff(yfit); % setting slope to 0 seems to not work as well and apparently is not needed

  peak = yfit(argmax) - ydata(argmax);

  %ceq = [trap, ydif(argmax), peak];
  ceq = [trap, peak];
end

function [f] = objective(x, ydata, coeff)
  yfit = model(x, coeff);
  f = sum((yfit - ydata) .^ 2);
end
end % extra end to close out function
% --- End Trevor edit
