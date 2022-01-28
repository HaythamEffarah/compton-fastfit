%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

    This script is intended for use in the MATLAB GUI. Please see
    README.md for code functionality.

    Required files and directories:
        /fwhm/
            fwhm.m
            fwonem.m
        get2Dfrom4D.m
        loadenergies.m
        polynomialspectrumfit.m
        /Perceptually uniform colormaps/

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SEC 1. Load anchor files into MATLAB workspace (code tested for three)

% Values taken from leading number in anchors_100k file names
energies = [20 60 100];

% Load x-ray output spectra into data cell
if ~exist('data', 'var')
    
    data{length(energies)} = [];
    for i = 1:length(energies)
        data{i} = loadenergies(sprintf('anchors_100k/%ikeV_6um_xyn.csv', energies(i)));
    end
    
end
%% SEC 2. Convert data into smooth fitted polynomials

polyOrder = 8; %CHOOSE POLYNOMIAL ORDER TO FIT TO ENERGY DISTRIBUTIONS

% Initialize variables and iterator bounds
[m, n, ~, ~] = size(data{1});
[~, o] = size(data);
dataFitted{length(energies)} = [];
muFit = zeros(m,n,length(energies),2);
fitCoeffs = zeros(m,n, polyOrder + 1, o);
testEBound = zeros(m,n,length(energies),2);

% Fit polynomials to anchor distributions
for k = 1:length(energies)
    [fitCoeffs(:,:,:,k), testEBound(:,:,k,:), muFit(:,:,k,:),...
        dataFitted{k}] = polynomialspectrumfit(data{k}, polyOrder);
end

%--------------------------------------------------------------------------
% DESIGN NOTE: polynomial spectrum fit is only really used for arrays, so 
% perhaps changing the function to accept an array rather than iterating 
% out here would make things cleaner and more intuitive
%--------------------------------------------------------------------------

%% SEC 3. Create interpolation map

% Initialize variables and iterator bounds, partially redundant to last section
[m, n, ~, ~] = size(data{1});
[~, o] = size(data);
fwonemArray = zeros(m, n, o);
eBound = zeros(m,n,o,2);
meanEnergyArray = zeros(m, n, o);
mu = zeros(m,n,3,2);
muPoly = zeros(m,n,polyOrder+1,2);
predictEBounds = zeros(m,n,3,2);
muBounds = zeros(m,n,2,2);
muPredict = zeros(size(muFit));
predictCoeffs = zeros(m, n, 3, 3); % 3 for # of coeffs, 3 for # of parameters to fit
polyOnPoly = zeros(m, n, polyOrder + 1, 3); %polyfit the fit coefficients


% Create normalized 2D data files for plotting
dataImageRaw = cellfun(@get2Dfrom4D,dataFitted,'UniformOutput',false);
dataImage = zeros(m, n, o);
dataImageNorm = zeros(m, n, o);

% Create normalized image files to see intensity distributions
for k = 1:o
    dataImage(:,:,k) = dataImageRaw{k};
    dataImageNorm(:,:,k) = dataImageRaw{k} / max(dataImageRaw{k}, [], 'all');
end


% Pixel by pixel calculation of fitting parameters
for i = 1:m
    for j = 1:n
        
        for k = 1:o
            % Get FW1%M for all files at current pixel
            [fwonemArray(i,j,k), eBound(i,j,k,1), eBound(i,j,k,2)] ...
                = fwonem(dataFitted{k}(i,j,:,1),dataFitted{k}(i,j,:,2));

            % Get mean energy for all files at current pixel
            meanEnergyArray(i,j,k) = sum(dataFitted{k}(i,j,:,1) ...
                .* dataFitted{k}(i,j,:,2)) ./ sum(dataFitted{k}(i,j,:,2));
            
        end
        
        % Solve for fitting coefficients at current pixel
        [predictCoeffs(i,j,:,1), ~, mu(i,j,1,:)] =  polyfit(energies, fwonemArray(i,j,:), 2);
        [predictEBounds(i,j,:,1), ~, muBounds(i,j,1,:)] =  polyfit(energies, testEBound(i,j,:,1), 2);
        [predictEBounds(i,j,:,2), ~, muBounds(i,j,2,:)] =  polyfit(energies, testEBound(i,j,:,2), 2);
        [predictCoeffs(i,j,:,2), ~, mu(i,j,2,:)] = polyfit(energies, meanEnergyArray(i,j,:), 2);
        muPredict(i,j,:,1) = polyfit(energies, muFit(i,j,:,1), 2);
        muPredict(i,j,:,2) = polyfit(energies, muFit(i,j,:,2), 2);
        [predictCoeffs(i,j,:,3), ~, mu(i,j,3,:)] = polyfit(energies, dataImageNorm(i,j,:), 2);
        
        % Fit polynomials to polynomial fitting coefficients
        for l = 1:polyOrder + 1
            polyOnPoly(i,j,l,:) = polyfit(energies, fitCoeffs(i,j,l,:), 2);
        end
        
    end
end

%--------------------------------------------------------------------------
% NOTE: At this point, all the information necessary to produce
% interpolated spectra is here. polyOnPoly is the matrix containing the
% quadratic fit of all polynomial coefficients.
%--------------------------------------------------------------------------
%% SEC 4. Predict spectrum at new energy using just polynomials

predictedEnergy = 20; % CHANGE THIS TO CHANGE PREDICTION
en = linspace(.625 * predictedEnergy, predictedEnergy + 2, 1000); % CHANGE BASED ON predictedEnergy; energy bins

% Initialize variables
newSpec = zeros(m, n, length(en), 2);
newSpecCoeffs = zeros(m, n, polyOrder+1);
newFWONEM = zeros(m, n);
newMeanEnergy = zeros(m, n);
newPhotonDist = zeros(m, n);
newEBounds = zeros(m, n, 2);
newMu = zeros(m,n,2);


for i = 1:m
    for j = 1:n
        newSpec(i, j, :, 1) = en;

        newFWONEM(i,j) = polyval(squeeze(predictCoeffs(i,j,:,1)), ...
            predictedEnergy, [], mu(i,j,1,:));
        newMu(i,j,1) = polyval(squeeze(muPredict(i,j,:,1)), ...
            predictedEnergy);
        newMu(i,j,2) = polyval(squeeze(muPredict(i,j,:,2)), ...
            predictedEnergy);
        newMeanEnergy(i,j) = polyval(squeeze(predictCoeffs(i,j,:,2)), ...
            predictedEnergy, [], mu(i,j,2,:));
        
        newPhotonDist(i,j) = polyval(squeeze(predictCoeffs(i,j,:,3)), ...
            predictedEnergy, [], mu(i,j,3,:));
        newEBounds(i,j,1) = polyval(squeeze(predictEBounds(i,j,:,1)), ...
            predictedEnergy, [], muBounds(i,j,1,:));
        newEBounds(i,j,2) = polyval(squeeze(predictEBounds(i,j,:,2)), ...
            predictedEnergy, [], muBounds(i,j,2,:));
         
        for l = 1:polyOrder + 1
            newSpecCoeffs(i,j,l) = ...
                polyval(squeeze(polyOnPoly(i,j,l,:)), predictedEnergy);
        end
         

        tmp = en >= newEBounds(i,j,1) & en <= newEBounds(i,j,2);
        newEn = en(tmp);
        
        muNew = [mean(newEn); std(newEn)]; %% useless

        newSpec(i, j, tmp, 2) = polyval(squeeze(newSpecCoeffs(i,j,:)), ...
            newEn, [], newMu(i,j,:));

        %Clean up stray values outside of single mode distribution
        [~, maxLoc] = max(newSpec(i,j,:,2));
        scanUp = maxLoc;
        scanDown = maxLoc;
        while newSpec(i,j,scanUp,2) >= (0.5 * newSpec(i,j,maxLoc,2))
            scanUp = scanUp + 1;
        end

        while (newSpec(i,j,scanUp,2) - newSpec(i,j,scanUp+1,2)) > 0 ...
                && newSpec(i,j,scanUp,2) > 0
            scanUp = scanUp + 1;
        end
        newSpec(i,j,scanUp:end,2) = 0;


        while newSpec(i,j,scanDown,2) >= (0.5 * newSpec(i,j,maxLoc,2))
            scanDown = scanDown - 1;
        end

        while (newSpec(i,j,scanDown,2) - newSpec(i,j,scanDown-1,2)) > 0 ...
                && newSpec(i,j,scanDown,2) > 0
            scanDown = scanDown - 1;
        end
        newSpec(i,j,1:scanDown,2) = 0;
    end
end

newSpec(newSpec < 0) = 0;
scaling = zeros(m,n);

for i = 1:m
    for j = 1:n
        power = trapz(squeeze(newSpec(i,j,:,1)), squeeze(newSpec(i,j,:,2)));
        scalefactor = newPhotonDist(j,i) / power;
        scaling(i,j) = scalefactor;
        newSpec(i,j,:,2) = newSpec(i,j,:,2) * scalefactor;
    end
end

% New spectrum generated as newSpec!


%% SEC 4A. Produce a set of LCS spectra

% Initialize variables
newSpec = zeros(m, n, 1000, 2);
newSpecCoeffs = zeros(m, n, polyOrder+1);
newFWONEM = zeros(m, n);
newMeanEnergy = zeros(m, n);
newPhotonDist = zeros(m, n);
newEBounds = zeros(m, n, 2);
newMu = zeros(m,n,2);
table1Energies = [25 40 55 65 80 95]; % energies chosen for Table 1
fastfitArray{size(table1Energies,2)} = [];
magicI = 1; % very sloppy, forgive me

for pp = table1Energies
    predictedEnergy = pp;
    en = linspace((pp * .625), pp + 2, 1000);

    
    for i = 1:m
        for j = 1:n
            newSpec(i, j, :, 1) = en;

            newFWONEM(i,j) = polyval(squeeze(predictCoeffs(i,j,:,1)), ...
                predictedEnergy, [], mu(i,j,1,:));
            newMu(i,j,1) = polyval(squeeze(muPredict(i,j,:,1)), ...
                predictedEnergy);
            newMu(i,j,2) = polyval(squeeze(muPredict(i,j,:,2)), ...
                predictedEnergy);
            newMeanEnergy(i,j) = polyval(squeeze(predictCoeffs(i,j,:,2)), ...
                predictedEnergy, [], mu(i,j,2,:));

            newPhotonDist(i,j) = polyval(squeeze(predictCoeffs(i,j,:,3)), ...
                predictedEnergy, [], mu(i,j,3,:));
            newEBounds(i,j,1) = polyval(squeeze(predictEBounds(i,j,:,1)), ...
                predictedEnergy, [], muBounds(i,j,1,:));
            newEBounds(i,j,2) = polyval(squeeze(predictEBounds(i,j,:,2)), ...
                predictedEnergy, [], muBounds(i,j,2,:));

            for l = 1:polyOrder + 1
                newSpecCoeffs(i,j,l) = ...
                    polyval(squeeze(polyOnPoly(i,j,l,:)), predictedEnergy);
            end


            tmp = en >= newEBounds(i,j,1) & en <= newEBounds(i,j,2);
            newEn = en(tmp);

            muNew = [mean(newEn); std(newEn)]; %% useless

            newSpec(i, j, tmp, 2) = polyval(squeeze(newSpecCoeffs(i,j,:)), ...
                newEn, [], newMu(i,j,:));

            %Clean up stray values outside of single mode distribution
            [~, maxLoc] = max(newSpec(i,j,:,2));
            scanUp = maxLoc;
            scanDown = maxLoc;
            while newSpec(i,j,scanUp,2) >= (0.5 * newSpec(i,j,maxLoc,2))
                scanUp = scanUp + 1;
            end

            while (newSpec(i,j,scanUp,2) - newSpec(i,j,scanUp+1,2)) > 0 ...
                    && newSpec(i,j,scanUp,2) > 0
                scanUp = scanUp + 1;
            end
            newSpec(i,j,scanUp:end,2) = 0;


            while newSpec(i,j,scanDown,2) >= (0.5 * newSpec(i,j,maxLoc,2))
                scanDown = scanDown - 1;
            end

            while (newSpec(i,j,scanDown,2) - newSpec(i,j,scanDown-1,2)) > 0 ...
                    && newSpec(i,j,scanDown,2) > 0
                scanDown = scanDown - 1;
            end
            newSpec(i,j,1:scanDown,2) = 0;
        end
    end

    newSpec(newSpec < 0) = 0;
    scaling = zeros(m,n);

    for i = 1:m
        for j = 1:n
            power = trapz(squeeze(newSpec(i,j,:,1)), squeeze(newSpec(i,j,:,2)));
            scalefactor = newPhotonDist(j,i) / power;
            scaling(i,j) = scalefactor;
            newSpec(i,j,:,2) = newSpec(i,j,:,2) * scalefactor;
        end
    end
    
    fastfitArray{magicI} = newSpec;
    magicI = magicI + 1;    

end

%% SEC 5. Find minimum energy bandwidth circular aperture
% Definitely takes longer than it needs to, as it's a very "dumb" scan, but
% it gets the job done in a few minutes.

% Initialize variables
en = zeros(1,1000);
newSpec = zeros(m, n, length(en), 2);
oldSpec = newSpec;
newSpecCoeffs = zeros(m, n, polyOrder+1);
newFWONEM = zeros(m, n);
newMeanEnergy = zeros(m, n);
newPhotonDist = zeros(m, n);
newEBounds = zeros(m, n, 2);
newMu = zeros(m,n,2);

circlePixels= zeros(m,n,50);

% Create circle mask
centerX = ceil(m / 2); 
centerY = ceil(m / 2);
% radius = 20;
for radius = 1:50
    [columnsInImage, rowsInImage] = meshgrid(1:m, 1:n);
    circlePixels(:,:,radius) = (rowsInImage - centerY).^2 ... 
        + (columnsInImage - centerX).^2 <= radius.^2;
end

% Define bandwidth of interest
E_center = 72;
E_halfWidth = E_center * .02;
E_max = E_center + E_halfWidth;
E_min = E_center - E_halfWidth;

% SET start energy and end energy such that eEND > eSTART
eSTART = 70;
eEND = 80;

% SET number of steps
steps = 1:100;

% Calculates step size in energy
deltaE = (eEND - eSTART) / length(steps);

fwhmTracker = zeros(length(steps),50);
integratedSpectrum = zeros(1000,length(steps),50);
newSpecGood = integratedSpectrum;
flux = fwhmTracker;
fluxNorm = fwhmTracker;

for pp = steps
    
    predictedEnergy = eSTART + (pp - 1) * deltaE;

    en = linspace((predictedEnergy * .75), predictedEnergy + 2, 1000);
    for i = 1:m
        for j = 1:n
            newSpec(i, j, :, 1) = en;
            oldSpec(i, j, :, 1) = en;


            newFWONEM(i,j) = polyval(squeeze(predictCoeffs(i,j,:,1)), predictedEnergy, [], mu(i,j,1,:));
            newMu(i,j,1) = polyval(squeeze(muPredict(i,j,:,1)), predictedEnergy);
            newMu(i,j,2) = polyval(squeeze(muPredict(i,j,:,2)), predictedEnergy);
            newMeanEnergy(i,j) = polyval(squeeze(predictCoeffs(i,j,:,2)), predictedEnergy, [], mu(i,j,2,:));


            newPhotonDist(i,j) = polyval(squeeze(predictCoeffs(i,j,:,3)), predictedEnergy, [], mu(i,j,3,:));
            newEBounds(i,j,1) = polyval(squeeze(predictEBounds(i,j,:,1)), predictedEnergy, [], muBounds(i,j,1,:));
            newEBounds(i,j,2) = polyval(squeeze(predictEBounds(i,j,:,2)), predictedEnergy, [], muBounds(i,j,2,:));

            for l = 1:polyOrder + 1
                newSpecCoeffs(i,j,l) = polyval(squeeze(polyOnPoly(i,j,l,:)), predictedEnergy);
            end


            tmp = en >= newEBounds(i,j,1) & en <= newEBounds(i,j,2);
            newEn = en(tmp);
            muNew = [mean(newEn); std(newEn)];
            oldSpec(i, j, tmp, 2) = polyval(squeeze(fitCoeffs(i,j,:,2)), newEn, [], muFit(i,j,2,:));
            newSpec(i, j, tmp, 2) = polyval(squeeze(newSpecCoeffs(i,j,:)), newEn, [], newMu(i,j,:));


        end
    end

    newSpec(newSpec < 0) = 0;

    for i = 1:m
        for j = 1:n
            power = trapz(squeeze(newSpec(i,j,:,1)), squeeze(newSpec(i,j,:,2)));
            scalefactor = newPhotonDist(j,i) / power;
            newSpec(i,j,:,2) = newSpec(i,j,:,2) * scalefactor;
        end
    end
    
    % Identify energy range for E_good photons
    logicalIndices = (en < E_max) & (en > E_min);
    
    for radius = 1:50
        % Get photons within radius of aperture
        newSpecAperture = squeeze(circlePixels(:,:,radius)) .* newSpec;
        
        % Get total integrated spectrum within aperture
        integratedSpectrum(:,pp,radius) = squeeze(sum(newSpecAperture(:,:,:,2), [1,2]));
        
        % Get integrated spectrum of photons that meet E_good requirements
        newSpecGood(logicalIndices,pp,radius) = squeeze(sum(newSpecAperture(:,:,logicalIndices,2), [1,2]));
    
        flux(pp,radius) = trapz(en, integratedSpectrum(:,pp,radius));
        [~,temp] = max(integratedSpectrum(:,pp,radius));
        fwhmTracker(pp,radius) = fwhm(squeeze(newSpec(1,1,:,1)), squeeze(integratedSpectrum(:,pp,radius)))./ en(temp) ;
    end
    
end

%Find total number of "good" photons at each aperture, at each energy
E_good = zeros(size(newSpecGood,2),radius);
E_bad = E_good;

for i = steps
    for radius = 1:50
        temporary = eSTART + (pp - 1) * deltaE;
        E_good(i,radius) = trapz(linspace(((i+69) * .75), (i+69) + 2, 1000), newSpecGood(:,i,radius));
        E_bad(i,radius) = trapz(linspace(((i+69) * .75), (i+69) + 2, 1000), integratedSpectrum(:,i,radius)) - E_good(i,radius);
    end
end


%% SEC 6. Find minimum energy bandwidth annular aperture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definitely takes longer than it needs to, as it's a very "dumb" scan, but
% it gets the job done.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize variables
en = zeros(1,1000);
newSpec = zeros(m, n, length(en), 2);
newSpecCoeffs = zeros(m, n, polyOrder+1);
newFWONEM = zeros(m, n);
newMeanEnergy = zeros(m, n);
newPhotonDist = zeros(m, n);
newEBounds = zeros(m, n, 2);
newMu = zeros(m,n,2);

circlePixels= zeros(m,n,51); circleBlocks = zeros(m,n,51);

% Create circle mask
centerX = ceil(m / 2); 
centerY = ceil(m / 2);

for radius = 1:51
    [columnsInImage, rowsInImage] = meshgrid(1:m, 1:n);
    circlePixels(:,:,radius) = (rowsInImage - centerY).^2 ... 
        + (columnsInImage - centerX).^2 <= (radius-1).^2;
    
end

circleBlock = ~circlePixels;

% Define bandwidth of interest
E_center = 72;
E_halfWidth = E_center * .02;
E_max = E_center + E_halfWidth;
E_min = E_center - E_halfWidth;

% SET start energy and end energy such that eEND > eSTART
eSTART = 70;
eEND = 80;

% SET number of steps
steps = 1:20;

% Calculates step size in energy
deltaE = (eEND - eSTART) / length(steps);

fwhmTracker = zeros(length(steps),51,51);
integratedSpectrum = zeros(1000,length(steps),51,51);
newSpecGood = integratedSpectrum;
E_good = zeros(size(newSpecGood,2),51,51);
E_bad = E_good;

for pp = steps
    
    predictedEnergy = eSTART + (pp - 1) * deltaE;

    en = linspace((predictedEnergy * .75), predictedEnergy + 2, 1000);
    for i = 1:m
        for j = 1:n
            newSpec(i, j, :, 1) = en;
            oldSpec(i, j, :, 1) = en;


            newFWONEM(i,j) = polyval(squeeze(predictCoeffs(i,j,:,1)), predictedEnergy, [], mu(i,j,1,:));
            newMu(i,j,1) = polyval(squeeze(muPredict(i,j,:,1)), predictedEnergy);
            newMu(i,j,2) = polyval(squeeze(muPredict(i,j,:,2)), predictedEnergy);
            newMeanEnergy(i,j) = polyval(squeeze(predictCoeffs(i,j,:,2)), predictedEnergy, [], mu(i,j,2,:));


            newPhotonDist(i,j) = polyval(squeeze(predictCoeffs(i,j,:,3)), predictedEnergy, [], mu(i,j,3,:));
            newEBounds(i,j,1) = polyval(squeeze(predictEBounds(i,j,:,1)), predictedEnergy, [], muBounds(i,j,1,:));
            newEBounds(i,j,2) = polyval(squeeze(predictEBounds(i,j,:,2)), predictedEnergy, [], muBounds(i,j,2,:));

            for l = 1:polyOrder + 1
                newSpecCoeffs(i,j,l) = polyval(squeeze(polyOnPoly(i,j,l,:)), predictedEnergy);
            end


            tmp = en >= newEBounds(i,j,1) & en <= newEBounds(i,j,2);
            newEn = en(tmp);
            muNew = [mean(newEn); std(newEn)];
            newSpec(i, j, tmp, 2) = polyval(squeeze(newSpecCoeffs(i,j,:)), newEn, [], newMu(i,j,:));


        end
    end

    newSpec(newSpec < 0) = 0;

    for i = 1:m
        for j = 1:n
            power = trapz(squeeze(newSpec(i,j,:,1)), squeeze(newSpec(i,j,:,2)));
            scalefactor = newPhotonDist(j,i) / power;
            newSpec(i,j,:,2) = newSpec(i,j,:,2) * scalefactor;
        end
    end
    
    % Identify energy range for E_good photons
    logicalIndices = (en <= E_max) & (en >= E_min);
    
    %Find max pixel radius that no longer contains any E_good photons
    goodMap = get2Dfrom4D(newSpec(:,:,logicalIndices,:));
    [row, col] = find(goodMap);
    if isempty(find(goodMap,1)) == 1
        continue
    end
    radiusMax = max([max(abs(51-col)), max(abs(51-row))]);
    radiusMin = floor(min(sqrt( (row - 51).^2 + (col - 51).^2)));
    
        for radius = radiusMin+1:radiusMax+1
            for radius2 = radiusMin+1:radiusMax+1
                if radius2 > radius
                    continue
                end
                % Get photons within radius of aperture
                newSpecAperture = squeeze(circlePixels(:,:,radius)) .* ...
                    squeeze(circleBlock(:,:,radius2)) .* newSpec;
                if sum(newSpecAperture(:,:,:,2)) == 0
                    continue
                end
                % Get total integrated spectrum within aperture
                %%% NOTE: This is only valid when mrad spacing is conserved for
                %%% both anchor and generated distributions
                integratedSpectrum(:,pp,radius,radius2) = squeeze(sum(newSpecAperture(:,:,:,2), [1,2]));

                % Get integrated spectrum of photons that meet E_good requirements
                newSpecGood(logicalIndices,pp,radius,radius2) = squeeze(sum(newSpecAperture(:,:,logicalIndices,2), [1,2]));

                % Define E_good and E_bad
                E_good(pp,radius,radius2) = trapz(en, newSpecGood(:,pp,radius,radius2).');
                E_bad(pp,radius,radius2) = trapz(en, ...
                    integratedSpectrum(:,pp,radius,radius2)) - E_good(pp,radius,radius2);
            end
        end
end

