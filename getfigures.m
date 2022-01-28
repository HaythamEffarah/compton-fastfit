% generate figures for manuscript

%% USAGE NOTE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each section in getfigures.m is titled as "Figure #." followed by a brief
% description. Under the title of each section are instructions on which
% parts of the "comptonfastfit.m" script need to be run to produce the
% figure in question. It's admittedly sloppy, but the code blocks in each
% section of this script require the workspace variables that are generated
% in the required sections (SEC.).
% 
% For best results, ONLY run the comptonfastfit.m sections that are listed
% in each of the following Figure code blocks. The robustness of the code
% in maintaining the workspace variables when running different
% combinations of sections has not been fully tested.
%   
%   Code Section Requirement Summary:
%       Figure 1. SEC. 1
%       Figure 2. SECS. 1-2
%       Figure 3. NONE
%       Figure 4. SECS. 1-2
%       Figure 5. SECS. 1-4
%       Figure 6. SECS. 1-4
%       Figure 7. SECS. 1-3, 5
%       Figure 8. SECS. 1-3, 5
%       Figure 9. SECS. 1-3, 6
%
% Have fun!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMKollkNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWk;,,lXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMNXWMMMMMMMMMMMMMMWk;,,lXMMMMMMMMMMMMMMWXXWMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMWOlcdXMMMMMMMMMMMMMWk;,,lXMMMMMMMMMMMMMNkc:xXMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMW0l;;oKMMMMMMMMMMMMWk;,,lXMMMMMMMMMMMMNd;,:kWMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMKo;;l0WMMMMMWWWMMWk;,,lXMMMWWMMMMMMXo;,cOWMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMMXd;,cOWMMWOolkWMWk;,,lXMWKdokXMMWKl;;l0WMMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMMMNx;,:kNMWk;,c0MWk;,,lXMNo;,oXMW0c;;oKWMMMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMMMMNk:,:xNMXl,;dNWk;,,lXWO:,:OWWOc;;dXMMMMMMMMMMMMMMMMMMMMMMM
% MMMMWNXWMMMMMMMMMMN00NMWOc;;dXWO:,:OWk;,,lXXo,;oXNk:,:xNMWXOKWMMMMMMMMMWXXWMMMMM
% MMMW0c:lx0XWMMMMMXo;:o0NW0c;;oKXo,,oXk;,,lKO:,:OXx;,:kNWKxc;c0WMMMMMNKkoc:xNMMMM
% MMMWKxl;,;cdkKWMMW0d:,;o0NKl;;l0O:,:xo;,,cxl,;d0d;,cONKx:;;lONMMWX0xl:,;coONMMMM
% MMMMMMNKko:;,:lx0NWW0o:;;oO0o;,cl:,;;;;;;;;;;;cl;,l0Kx:;;lONWNKkoc;,:lx0XWMMMMMM
% MMMMMMMMMWXOdl;,;cdOKN0o:,;oo:;;::::::::::::::::;;ldc,;lOXX0xl:;;cokKNMMMMMMMMMM
% MMMMMMMMMMMMWN0ko:;;:oxOkl;;;,cxOOOOOOOOOOOOOO0Oo;,,,cxOkoc;,:lxOXWMMMMMMMMMMMMM
% MMMMMMMMMMMNKXNWWXOdc;,;cc;;;,:;;;;;;;;;;;;;;;;;;;;;;:l:,;cokKNWWXKXWMMMMMMMMMMM
% MMMMMMMMMMXo;:lodxO0Okdc;;;;;;:dddddddddddddddddl;;;;;;:okO0Okdolc;c0MMMMMMMMMMM
% MMMMMMMMMMNOdol:;,,;:cl:;;;;,;o0KKKKKKKKKKKKKKKKk:,;;;,:ll:;,,;:codxKMMMMMMMMMMM
% MMMMMMMMMMMMMWNXKOkxdl:;;;;;;;:looood0KKKKkoooooc;;;;;;,;codxO0XNWMMMMMMMMMMMMMM
% kxxxxxxxxxxxxxxxxxxxxdc;;;;;;;;,,,,,cOKKKKd;,;,,;;;;;;;,;oxxxxxxxxxxxxxxxxxxxxx0
% :;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;,cOKKKKd;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;o
% KKKKKKKKKKKKKKKKK00kxdc;;;;;;;;;;;;,cOKKKKd;;;;;;;;;;;;,;lxkO0KKKKKKKKKKKKKKKKKX
% MMMMMMMMMMMNXKOkxdlc:;;;;;;;;;;;;;;,cOKKKKd;;;;;;;;;;;;;,,;:lodxO0KNWMMMMMMMMMMM
% MMMMMMMMMMXd:;;,;:cldxxc;;;;;;;;;;;,cOKKKKd;;;;;;;;;;;,:dxdol:;;,;:l0MMMMMMMMMMM
% MMMMMMMMMMNOodkO0KXKOdl;,,;;;;;;;;;,cOKKKKd;,;;;;;;;;;,;cok0XXKOkxoxXMMMMMMMMMMM
% MMMMMMMMMMMMWMMNKko:;,:lxo;;;;;;;;;,cOKKKKd;,;;;,,;;;cxoc;,:lx0XWMWMMMMMMMMMMMMM
% MMMMMMMMMMMWXOdl:,;cdk00d:;;:;;;;;;,cOKKKKd;,;;;,;::;;oOK0xl:,;cokKNMMMMMMMMMMMM
% MMMMMMMMNKkoc;,:lx0NWKd:,;lkx:,;:;;;:oddddc;;;;;;;oOd:,;oONNKkoc;,:lx0XWMMMMMMMM
% MMMMWXOxl:,;cdkKNMWKd:;;lOXk:,:xkc,;lc;,,:l:,;dkc,;dK0d:,;oONMWXOxl:,;cdkKWMMMMM
% MMMW0c;;:lx0XWMMMXx:;;lONNx:,:kXx;,lKk;,,lKx;,lK0c,;oKWKd:,;oKMMMMNKkoc;,:xNMMMM
% MMMMXkdkKWMMMMMMMNkloONMXd;,cOWKc,;kWk;,,lK0c,;xNKl;;lKWWKdldKMMMMMMMWXOxd0WMMMM
% MMMMMMMMMMMMMMMMMMMWWMMKo;,l0WNd;,lKWk;,,lKWx;,cKMXo;;c0WMWNWMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMMMWKl;;lKWM0c,;kWWk;,,lXMKl,;xNMXd;,cOWMMMMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMMWOc;;oXMMWx;,oXMWk;,,lXMWk;;lKMMNx:,:xNMMMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMWk:,;dNMMMWX0OKWMWk;,,lXMMN00XWMMMWk:,;dNMMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMNx:,:xNMMMMMMMMMMMWk;,,lXMMMMMMMMMMMWOc;;oXMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMWk:,:kWMMMMMMMMMMMMWk;,,lXMMMMMMMMMMMMW0l,;oXMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMWXOx0WMMMMMMMMMMMMMWk;,,lXMMMMMMMMMMMMMMXxx0WMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWk;,,lXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
% MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWk;,,lXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

%% Figure 1. Angle-correlated spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% SEC. 1 in comptonfastfit.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Illustrate flux density at 60 keV and then create lineouts of four
% different regions (the whole thing, the center, off-center, and a thin
% band

% Create plotting meshgrid
x = linspace(-5,5,101);
y = x;
[X, Y] = meshgrid(x,y);

%create graphics matrix
%In order: flux, a, b, c, d, overlapping
H = gobjects(6,3);

for k = 1:size(H,1)
    H(k,1) = figure;
    H(k,2) = axes('Parent',H(k,1));
end

%%% H(1,2) 60 keV photon distribution plot

fig1_sub1 = get2Dfrom4D(data{2});
fig1_sub1 = fig1_sub1 / max(fig1_sub1, [], 'all');

% Create ax(1) plot
surf(H(1,2), X, Y, fig1_sub1, 'EdgeColor', 'none');

% Set plot details
view(H(1,2),[0 90]);
colorbar(H(1,2))
h1 = colorbar(H(1,2));
set(get(h1,'label'),'string','Relative flux density');
xlabel(H(1,2), '\theta_x (mrad)')
ylabel(H(1,2), '\theta_y (mrad)')
H(1,2).XTick = -4:2:4;
H(1,2).YTick = -4:2:4;
H(1,2).FontSize = 24;
pbaspect(H(1,2),[1 1 1])

colormap(H(1,2),'jet')
shading(H(1,2), 'interp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% H(2,2) Lineout of entire integrated spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1_x = squeeze(data{2}(1,1,:,1));
spectrum_a = squeeze(sum(data{2}(:,:,:,2), [1,2]));
spectrum_a_norm = spectrum_a ./ max(spectrum_a);

patch(H(2,2), [fig1_x fliplr(fig1_x)], [spectrum_a_norm ...
    fliplr(spectrum_a_norm)], [94 94 94]./256,'FaceAlpha',0.5)

xlabel(H(2,2), 'Energy (keV)')
ylabel(H(2,2), 'Relative flux')
H(2,2).FontSize = 24;
xlim(H(2,2),[45 62]);
ylim(H(2,2), [0 1.1]);
pbaspect(H(2,2),[1 1 1]);
box(H(2,2), 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% H(3,2) Lineout of annulus with 2.0 < r < 2.3 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create circle mask
[m, n, ~, ~] = size(data{2});
circlePixels= zeros(m,n,51); circleBlocks = zeros(m,n,51);
centerX = ceil(m / 2); 
centerY = ceil(m / 2);
for radius = 1:51
    [columnsInImage, rowsInImage] = meshgrid(1:m, 1:n);
    circlePixels(:,:,radius) = (rowsInImage - centerY).^2 ... 
        + (columnsInImage - centerX).^2 <= (radius-1).^2;
    
end
circleBlock = ~circlePixels;

spectrum_b = squeeze(circlePixels(:,:,25)) .* ...
    squeeze(circleBlock(:,:,19)) .* data{2};

spectrum_b = squeeze(sum(spectrum_b(:,:,:,2), [1,2]));
spectrum_b_norm = spectrum_b ./ max(spectrum_a);

patch(H(3,2), [fig1_x fliplr(fig1_x)], [spectrum_b_norm ...
    fliplr(spectrum_b_norm)], [4 51 255]./256,'FaceAlpha',0.5)

xlabel(H(3,2), 'Energy (keV)')
ylabel(H(3,2), 'Relative flux')
H(3,2).FontSize = 24;
xlim(H(3,2),[45 62]);
ylim(H(3,2), [0 1.1]);
pbaspect(H(3,2),[1 1 1]);
box(H(3,2), 'on');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% H(4,2) Lineout of circlle with 0 < r < 0.7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spectrum_c = squeeze(circlePixels(:,:,10)) .* data{2};

spectrum_c = squeeze(sum(spectrum_c(:,:,:,2), [1,2]));
spectrum_c_norm = spectrum_c ./ max(spectrum_a);

patch(H(4,2), [fig1_x fliplr(fig1_x)], [spectrum_c_norm ...
    fliplr(spectrum_c_norm)], [0 250 146]./256,'FaceAlpha',0.5)

xlabel(H(4,2), 'Energy (keV)')
ylabel(H(4,2), 'Relative flux')
H(4,2).FontSize = 24;
xlim(H(4,2),[45 62]);
ylim(H(4,2), [0 1.1]);
pbaspect(H(4,2),[1 1 1]);
box(H(4,2), 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% H(5,2) Lineout of offset circle (-3.5,3) with r < 0.7 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create off-center circle mask
[m, n, ~, ~] = size(data{2});
offsetCirclePixels= zeros(m,n,51);
centerX = 81;
centerY = 86;
for radius = 1:51
    [columnsInImage, rowsInImage] = meshgrid(1:m, 1:n);
    offsetCirclePixels(:,:,radius) = (rowsInImage - centerY).^2 ... 
        + (columnsInImage - centerX).^2 <= (radius-1).^2;
    
end

spectrum_d = squeeze(offsetCirclePixels(:,:,10)) .* data{2};
spectrum_d = squeeze(sum(spectrum_d(:,:,:,2), [1,2]));
spectrum_d_norm = spectrum_d ./ max(spectrum_a);

patch(H(5,2), [fig1_x fliplr(fig1_x)], [spectrum_d_norm ...
    fliplr(spectrum_d_norm)], [255 38 0]./256,'FaceAlpha',0.5)

xlabel(H(5,2), 'Energy (keV)')
ylabel(H(5,2), 'Relative flux')
H(5,2).FontSize = 24;
xlim(H(5,2),[45 62]);
ylim(H(5,2), [0 1.1]);
pbaspect(H(5,2),[1 1 1]);
box(H(5,2), 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% H(6,2) All lineouts overlapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

patch(H(6,2), [fig1_x fliplr(fig1_x)], [spectrum_a ./ max(spectrum_a) fliplr(spectrum_a ./ max(spectrum_a))], ...
    [94 94 94]./256,'FaceAlpha',0.5)
hold on
patch(H(6,2), [fig1_x fliplr(fig1_x)], [spectrum_b ./ max(spectrum_a) fliplr(spectrum_b ./ max(spectrum_a))], ...
    [4 51 255]./256,'FaceAlpha',0.5)
patch(H(6,2), [fig1_x fliplr(fig1_x)], [spectrum_c ./ max(spectrum_a) fliplr(spectrum_c ./ max(spectrum_a))], ...
    [0 250 146]./256,'FaceAlpha',0.5)
patch(H(6,2), [fig1_x fliplr(fig1_x)], [spectrum_d ./ max(spectrum_a) fliplr(spectrum_d ./ max(spectrum_a))], ...
    [255 38 0]./256,'FaceAlpha',0.5)

xlabel(H(6,2), 'Energy (keV)')
ylabel(H(6,2), 'Relative flux')
H(6,2).FontSize = 24;
xlim(H(6,2),[45 62]);
pbaspect(H(6,2),[1 1 1]);
box(H(6,2), 'on');

%% Figure 2. Compare local distribution over a range of observation angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% SECS. 1-2 in comptonfastfit.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create figure and axes
f3 = figure;
ax3 = axes('Parent',f3);

%Plot energy distribution near center
% plot(ax(1), squeeze(data2(51,51,:,1)), squeeze(data2(51,51,:,2)) );
hold(ax3,'off')

x = squeeze(dataFitted{2}(1,1,:,1));
plot(ax3, x, ...
    squeeze(dataFitted{2}(51,58,:,2)) ./ ...
    max(squeeze(dataFitted{2}(51,58,:,2))), ...
    'LineWidth',2);

hold(ax3,'on')


plot(ax3, x, ...
    squeeze(dataFitted{2}(51,63,:,2)) ./ ...
    max(squeeze(dataFitted{2}(51,58,:,2))), ...
    'LineWidth',2);

plot(ax3, x, ...
    squeeze(dataFitted{2}(51,73,:,2)) ./ ...
    max(squeeze(dataFitted{2}(51,58,:,2))), ...
    'LineWidth',2);

hold(ax3,'on')

xlabel(ax3, 'Energy (keV)')
ylabel(ax3, 'X-ray flux (arb. units)')
title(ax3, 'Local X-ray variation over \theta_x (\theta_y = 0)');
ax3.FontSize = 24;
legend(ax3,'0.7 mrad', '1.2 mrad', '2.2 mrad', ...
    'Location','northwest')
ylim(ax3,[0 1.1]);
xlim(ax3,[57 60]);
pbaspect(ax3,[1 1 1]);
% f3.Position = [100 100 700 700];

%% Figure 3. Energy and relative flux density scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% NONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create array of directly simulated spectra from 20keV - 100keV
if ~exist('spectrumArray', 'var')
    fileNames = dir('anchors_100k/*.csv');
    num = size(fileNames,1);
    files = strings(num, 1);
    fileNumbers = zeros(num,1);

    for i = 1:num
        files(i) = convertCharsToStrings((fileNames(i).name));
        fileNumbers(i) = sscanf(files(i),'%i');
    end

    [~,sorter] = sort(fileNumbers);
    files = files(sorter(1:5));
    
    spectrumArray{num} = [];

    for i = 1:5
        spectrumArray{i} = loadenergies(files(i));
    end
end

%Create arrays of normalized flux density images
[m, n, ~, ~] = size(spectrumArray{num});
spectrumImageArray = zeros(m,n,num);

for i = 1:num
    spectrumImageArray(:,:,i) = get2Dfrom4D(spectrumArray{i});
end

normImageArray = zeros(size(spectrumImageArray));
for i = 1:num
    normImageArray(:,:,i) = spectrumImageArray(:,:,i) / sum(spectrumImageArray(:,:,i), 'all');
end

%Plot relative flux density, y = 0 mrad, x = (-5) -> 0 mrad
% Note that MATLAB indexes are (row, column), so scanning over x would be
% scanning over column while holding row constant
x = sort(fileNumbers);
normImageArrayToPlot = squeeze(normImageArray(50,1:50,:));
figure1 = figure;
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'off');
h = plot(x,normImageArrayToPlot ./ normImageArrayToPlot(50,1),'LineWidth',2,'Parent',axes1);
set(h, {'color'}, num2cell(jet(50),2));
ylabel({'Relative X-ray flux density'});
xlabel({'E_{CE} (keV)'});
set(axes1,'FontSize',24);
pbaspect(axes1,[1 1 1])
figure1.Position = [100 100 500 500];

%Plot relative flux density, x = 0 mrad, y = (-5) -> 0 mrad
% Note that MATLAB indexes are (row, column), so scanning over y would be
% scanning over row while holding column constant
x = sort(fileNumbers);
normImageArrayToPlot = squeeze(normImageArray(1:50,50,:));

%Plot relative flux density
figure1 = figure;
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'off');
h = plot(x,normImageArrayToPlot ./ normImageArrayToPlot(50,1),'LineWidth',2,'Parent',axes1);
colors = num2cell(jet(129),2);
set(h, {'color'}, colors(129-49:end));
ylabel({'Relative X-ray flux density'});
xlabel({'E_{CE} (keV)'});
set(axes1,'FontSize',24);
pbaspect(axes1,[1 1 1])
figure1.Position = [100 100 500 500];


%Create arrays of mean energies
[m, n, ~, ~] = size(spectrumArray{num});
modeEnergyArray = zeros(m,n,num);
meanEnergyArray = zeros(m,n,num);

for i = 1:m
    for j = 1:n
        for k = 1:num
        [~,modeLoc] = max(spectrumArray{k}(i,j,:,2));
        modeEnergyArray(i,j,k) = spectrumArray{k}(i,j,modeLoc,1);
        meanEnergyArray(i,j,k) = sum(spectrumArray{k}(i,j,:,1) .* spectrumArray{k}(i,j,:,2)) ./ sum(spectrumArray{k}(i,j,:,2));
        end
    end
end


%Plot mean energies, x = 0 mrad, y = (-5) -> 0 mrad
% Note that MATLAB indexes are (row, column), so scanning over y would be
% scanning over row while holding column constant
x = sort(fileNumbers);
meanEnergyToPlot = squeeze(meanEnergyArray(1:50,50,:));

figure1 = figure;
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'off');
h = plot(x,meanEnergyToPlot,'LineWidth',2,'Parent',axes1);
colors = num2cell(jet(129),2);
set(h, {'color'}, colors(129-49:end));
ylabel({'Local mean energy (keV)'});
xlabel({'E_{CE} (keV)'});
set(axes1,'FontSize',24);
pbaspect(axes1,[1 1 1])
figure1.Position = [100 100 500 500];

%Plot mean energies, y = 0 mrad, x = (-5) -> 0 mrad
% Note that MATLAB indexes are (row, column), so scanning over x would be
% scanning over column while holding row constant
x = sort(fileNumbers);
meanEnergyToPlot = squeeze(meanEnergyArray(50,1:50,:));
figure1 = figure;
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'off');
h = plot(x,meanEnergyToPlot,'LineWidth',2,'Parent',axes1);
set(h, {'color'}, num2cell(jet(50),2));
ylabel({'Local mean energy (keV)'});
xlabel({'E_{CE} (keV)'});
set(axes1,'FontSize',24);
pbaspect(axes1,[1 1 1])
figure1.Position = [100 100 500 500];

%% Figure 4. Image subtraction of anchor distribution (60 keV).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% SECS. 1-2 in comptonfastfit.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create plotting meshgrid
x = linspace(-5,5,101);
y = x;
[X, Y] = meshgrid(x,y);

%create figure and axes
f1 = figure;
f1.Position = [100 100 2000 700];

subplot(1,2,2)
ax = axes('Parent',f1);
for k = 1:2
    ax(k) = subplot(1,2,k);

end

% Redundantly re-define data files
data_compare{1} = data{2};
data_compare{2} = dataFitted{2};

% Create normalized 2D data files for plotting
directSimImageSet = cellfun(@get2Dfrom4D,data_compare,'UniformOutput',false);
for i = 1:size(directSimImageSet,2)
    dataImageSetRaw{i} = directSimImageSet{i};
    directSimImageSet{i} = directSimImageSet{i} / max(directSimImageSet{i}, [], 'all');
    numPixels(i) = numel(directSimImageSet{i});
end

%%% AX(1) INTENSITY SUBTRACTION PLOT
fluxDiff = directSimImageSet{1} - directSimImageSet{2};
fluxAPE = (abs(fluxDiff) ./ directSimImageSet{1}) .* 100;

% mean absolute percentage error
fluxMAPE = sum(fluxAPE , 'all') / numPixels(1) ;

% Create ax(1) plot
surf(ax(1), X, Y, fluxAPE, 'EdgeColor', 'none');

% Set plot details
view(ax(1),[0 90]);
title(ax(1), 'Flux');
colorbar(ax(1))
h1 = colorbar(ax(1));
set(get(h1,'label'),'string','APE (%)');
% caxis(ax(1), [0, 1]);
xlabel(ax(1), '\theta_x (mrad)')
ylabel(ax(1), '\theta_y (mrad)')
ax(1).XTick = -4:2:4;
ax(1).YTick = -4:2:4;
ax(1).FontSize = 24;
pbaspect(ax(1),[1 1 1])


%%% AX(2) MEAN ENERGY SUBTRACTION PLOT

%Plot mean energy pixel-by-pixel error 2D map
[m, n, ~, ~] = size(data_compare{1});
meanEnergyArray = zeros(m,n,2);

for i = 1:m
    for j = 1:n
        for k = 1:2
            [~,modeLoc] = max(data_compare{k}(i,j,:,2));
            meanEnergyArray(i,j,k) = sum(data_compare{k}(i,j,:,1) .* data_compare{k}(i,j,:,2)) ./ sum(data_compare{k}(i,j,:,2));
        end
    end
end

energyDiff = abs(squeeze((meanEnergyArray(:,:,1) - meanEnergyArray(:,:,2))));
energyAPE = (abs(energyDiff) ./ squeeze(meanEnergyArray(:,:,1))) .* 100;

% mean absolute percentage error
energyMAPE = sum(energyAPE , 'all') / numPixels(1) ;

surf(ax(2), X, Y, energyAPE, 'EdgeColor', 'none');
view(ax(2),[0 90]);
title(ax(2), 'Mean energy');
colorbar(ax(2))
h2 = colorbar(ax(2));
% caxis(ax(2), [0, 1]);
set(get(h2,'label'),'string','APE (%)');
xlabel(ax(2), '\theta_x (mrad)')

ax(2).XTick = -4:2:4;
ax(2).YTick = -4:2:4;
ax(2).FontSize = 24;
pbaspect(ax(2),[1 1 1])

%Set colormaps
colormap(ax(1),'inferno')
colormap(ax(2),'inferno')


%Shading
shading(ax(1), 'interp')
shading(ax(2), 'interp')

sprintf('Relative Intensity MAPE: %0.5e', fluxMAPE)
sprintf('Mean Energy MAPE: %0.5e', energyMAPE)

%% Figure 5. Image subtraction of interpolation (80 keV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% SECS. 1-4 in comptonfastfit.m
%
% NOTE: To reproduce Fig. 5 exactly, please make sure SEC. 4 is used to
% generate a spectrum at 80 keV, which is written as default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save time for loading in direct simulation test data
if ~exist('data_tester', 'var')
%     data_tester = loadenergies('table1/25keV_6um_xyn.csv');
%     data_tester = loadenergies('anchors_100k/80keV_6um_xyn.csv');
    data_tester = loadenergies('anchors_100k/20keV_6um_xyn.csv');
end

% Create plotting meshgrid
x = linspace(-5,5,101);
y = x;
[X, Y] = meshgrid(x,y);

%create figure and axes
f1 = figure;
f1.Position = [100 100 2000 700];

% ax = tight_subplot(1,3,[.03],[.1 .01],[.01 .01]);

subplot(1,2,2)
ax = axes('Parent',f1);
for k = 1:2
    ax(k) = subplot(1,2,k);

end

%%% AX(1) INTENSITY SUBTRACTION PLOT

data_compare{1} = data_tester;
data_compare{2} = newSpec;

% Create normalized 2D data files for plotting
directSimImageSet = cellfun(@get2Dfrom4D,data_compare,'UniformOutput',false);
for i = 1:size(directSimImageSet,2)
    dataImageSetRaw{i} = directSimImageSet{i};
    directSimImageSet{i} = directSimImageSet{i} / max(directSimImageSet{i}, [], 'all');
    numPixels(i) = numel(directSimImageSet{i});
end

fluxDiff = directSimImageSet{1} - directSimImageSet{2};
fluxAPE = (abs(fluxDiff) ./ directSimImageSet{1}) .* 100;

% mean absolute percentage error
fluxMAPE = sum(fluxAPE , 'all') / numPixels(1) ;

% Create ax(1) plot
surf(ax(1), X, Y, fluxAPE, 'EdgeColor', 'none');

% Set plot details
view(ax(1),[0 90]);
title(ax(1), 'Flux');
colorbar(ax(1))
h1 = colorbar(ax(1));
set(get(h1,'label'),'string','APE (%)');
xlabel(ax(1), '\theta_x (mrad)')
ylabel(ax(1), '\theta_y (mrad)')
ax(1).XTick = -4:2:4;
ax(1).YTick = -4:2:4;
ax(1).FontSize = 24;
pbaspect(ax(1),[1 1 1])


%%% AX(2) MEAN ENERGY SUBTRACTION PLOT

[m, n, ~, ~] = size(data_compare{1});
meanEnergyArray = zeros(m,n,2);

for i = 1:m
    for j = 1:n
        for k = 1:2
            [~,modeLoc] = max(data_compare{k}(i,j,:,2));
            meanEnergyArray(i,j,k) = sum(data_compare{k}(i,j,:,1) .* data_compare{k}(i,j,:,2)) ./ sum(data_compare{k}(i,j,:,2));
        end
    end
end

energyDiff = abs(squeeze((meanEnergyArray(:,:,1) - meanEnergyArray(:,:,2))));
energyAPE = (abs(energyDiff) ./ squeeze(meanEnergyArray(:,:,1))) .* 100;

% mean absolute percentage error
energyMAPE = sum(energyAPE , 'all') / numPixels(1) ;

surf(ax(2), X, Y, energyAPE, 'EdgeColor', 'none');
view(ax(2),[0 90]);
title(ax(2), 'Mean energy');
colorbar(ax(2))
h2 = colorbar(ax(2));
set(get(h2,'label'),'string','APE (%)');
xlabel(ax(2), '\theta_x (mrad)')

ax(2).XTick = -4:2:4;
ax(2).YTick = -4:2:4;
ax(2).FontSize = 24;
pbaspect(ax(2),[1 1 1])

%Set colormaps
colormap(ax(1),'inferno')
colormap(ax(2),'inferno')

%Shading
shading(ax(1), 'interp')
shading(ax(2), 'interp')

sprintf('Relative Intensity MAPE: %0.5e', fluxMAPE)
sprintf('Mean Energy MAPE: %0.5e', energyMAPE)

%% Table 1. MAPE analysis over several energies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% SECS. 1-3, 4A in comptonfastfit.m 
%
% NOTE: Things might get funky if you run SEC. 4 and SEC. 4A so beware
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save time for loading in direct simulation test data
% Create array of directly simulated spectra from 20keV - 100keV
if ~exist('directSimArray', 'var')
    fileNames = dir('table1/*.csv');
    num = size(fileNames,1);
    files = strings(num, 1);
    fileNumbers = zeros(num,1);

    for i = 1:num
        files(i) = convertCharsToStrings((fileNames(i).name));
        fileNumbers(i) = sscanf(files(i),'%i');
    end

    [~,sorter] = sort(fileNumbers);
    files = files(sorter(1:6));
    
    directSimArray{num} = [];

    for i = 1:num
        directSimArray{i} = loadenergies(files(i));
    end
end

MAPEtable = zeros(6,2);

% Create 2D images that only consider total X-ray photons per pixel
directSimImageSet = cellfun(@get2Dfrom4D,directSimArray,'UniformOutput',false);
fastfitImageSet = cellfun(@get2Dfrom4D,fastfitArray,'UniformOutput',false);

% Calculate flux MAPE for each pair of spectra
for i = 1:6
    directSimImageSet{i} = directSimImageSet{i} / max(directSimImageSet{i}, [], 'all');
    fastfitImageSet{i} = fastfitImageSet{i} / max(fastfitImageSet{i}, [], 'all');
    
    numPixels(i) = numel(directSimImageSet{i});
    fluxDiff = directSimImageSet{i} - fastfitImageSet{i};
    fluxAPE = (abs(fluxDiff) ./ directSimImageSet{i}) .* 100;
    fluxMAPE = sum(fluxAPE , 'all') / numPixels(i);
    MAPEtable(i,1) = fluxMAPE;
end


% Do the same as above but with mean energies
[m, n, ~, ~] = size(directSimArray{1});
fastfitMeanEnergyArray = zeros(m,n,6);
directSimMeanEnergyArray = zeros(m,n,6);

for i = 1:m
    for j = 1:n
        for k = 1:6
            fastfitMeanEnergyArray(i,j,k) = ...
                sum(fastfitArray{k}(i,j,:,1) .* ...
                fastfitArray{k}(i,j,:,2)) ./ sum(fastfitArray{k}(i,j,:,2));
            
            directSimMeanEnergyArray(i,j,k) = ...
                sum(directSimArray{k}(i,j,:,1) .* ...
                directSimArray{k}(i,j,:,2)) ./ sum(directSimArray{k}(i,j,:,2));
            
            energyDiff = abs(squeeze((directSimMeanEnergyArray(:,:,k) - fastfitMeanEnergyArray(:,:,k))));
            energyAPE = (abs(energyDiff) ./ squeeze(directSimMeanEnergyArray(:,:,k))) .* 100;
            energyMAPE = sum(energyAPE , 'all') / numPixels(1);
            MAPEtable(k,2) = energyMAPE;
        end
    end
end

disp(MAPEtable)

%% Figure 6. Zoom in lineout plots for interpolation at 80 keV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% SEC. 1-4 in comptonfastfit.m
%
% NOTE: To reproduce Fig. 6 exactly, please make sure SEC. 4 is used to
% generate a spectrum at 80 keV, which is the default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save time for loading in 80keV data
if ~exist('data_tester', 'var')
    data_tester = loadenergies('anchors_100k/80keV_6um_xyn.csv');
end

%create figure and axes
f3 = figure;
ax3 = axes('Parent',f3);

%Plot energy distribution near center
% plot(ax(1), squeeze(data2(51,51,:,1)), squeeze(data2(51,51,:,2)) );
hold(ax3,'off')

plot(ax3, squeeze(data_tester(51,53,:,1)), squeeze(data_tester(51,53,:,2))...
    / trapz(squeeze(data_tester(51,53,:,1)), squeeze(data_tester(51,53,:,2))),...
    'LineStyle','None','Marker','*','MarkerEdgeColor','red');
hold(ax3,'on')
plot(ax3, squeeze(newSpec(51,53,:,1)), squeeze(newSpec(51,53,:,2))...
    / trapz(squeeze(newSpec(51,53,:,1)), squeeze(newSpec(51,53,:,2))),...
    'LineWidth',1.5,'Color','red');

plot(ax3, squeeze(data_tester(51,62,:,1)), squeeze(data_tester(51,62,:,2))...
    / trapz(squeeze(data_tester(51,62,:,1)), squeeze(data_tester(51,62,:,2))),...
    'LineStyle','None','Marker','*', 'MarkerEdgeColor','blue');
plot(ax3, squeeze(newSpec(51,62,:,1)), squeeze(newSpec(51,62,:,2))...
    / trapz(squeeze(newSpec(51,62,:,1)), squeeze(newSpec(51,62,:,2))),...
    'LineWidth',1.5,'Color','blue');

xlabel(ax3, 'Energy (keV)')
ylabel(ax3, 'X-ray flux (arb. units)')
title(ax3, 'Local energy spectra for E_{CE} = 80 keV');
ax3.FontSize = 16;
legend(ax3,'LCS Code','FastFit',...
    'LCS Code','FastFit','Location','northwest')
xlim(ax3,[78 80]);
pbaspect(ax3,[1 1 1]);

%% Figure 7. Plot E_good - E_bad vs. aperture vs. E_gamma as a surface plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% SECS. 1-3 and SEC. 5 in comptonfastfit.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure;
ax = axes('Parent',f);

x = .1:.1:5;
y = eSTART:deltaE:eEND-deltaE;
[X, Y] = meshgrid(x,y);


sc = surfc(ax, X, Y, ((E_good - E_bad) ./ (E_good + E_bad)), 'EdgeColor', 'none');
ax.ZLim(2) = 15;
sc(2).ZLocation = 'zmax';
sc(2).LineColor = 'black';
sc(2).ShowText = 'on';
clabel([],sc(2),'FontSize',13, 'FontWeight', 'bold', 'Color', 'black')
sc(2).LineWidth = 2;
sc(2).LevelList = [0.8,0.8];
xl = xline(ax, 2.8,'--','2.8 mrad', 'LineWidth', 2.5, 'Color', 'white');
xl.FontSize = 16;
% xlim = 0.1:5;
pbaspect(ax,[1 1 1])
shading interp
colormap jet

view(ax,[0 90]);
title(ax, '');
colorbar(ax)
h1 = colorbar(ax);
set(get(h1,'label'),'string','( N_{good} - N_{bad} ) / N_{total}');
xlabel(ax, 'aperture radius (mrad)')
ylabel(ax, 'E_{CE} (keV)')
ax.FontSize = 20;

%% Figure 8. CIRCULAR aperture, compare Direct Simulation to FastFit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% SECS. 1-3 and SEC. 5 in comptonfastfit.m
%
% NOTE: Running SEC. 6 will overwrite some variables defined in SEC. 5 and
% vice versa.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('filesCircular', 'var')
    filesCircular = dir('circular_aperture_check/*.csv'); 
    circularData{length(filesCircular)} = [];
    for i=1:length(filesCircular)
        circularData{i} = loadenergies(append('circular_aperture_check/',filesCircular(i).name));
    end
end

dataAperture{length(filesCircular)} = [];
goodIndices = dataAperture;
dataIntegrated = dataAperture;
dataGood = dataAperture;
E_good_DATA = zeros(length(filesCircular),1);
E_bad_DATA = E_good_DATA;

for i = 1:length(filesCircular)
    dataAperture{i} = squeeze(circlePixels(:,:,28)) .* circularData{i};
    dataAperture{i}(:,:,:,1) = circularData{i}(:,:,:,1);
    
    goodIndices{i} = (squeeze(circularData{i}(1,1,:,1)) < E_max) & (squeeze(circularData{i}(1,1,:,1)) > E_min);
    
    dataIntegrated{i} = squeeze(sum(dataAperture{i}(:,:,:,2), [1,2]));
    
    dataGood{i} = zeros(size(circularData{i}(1,1,:,1),3),1);
    dataGood{i}(goodIndices{i}) = squeeze(sum(dataAperture{i}(:,:,goodIndices{i},2), [1,2]));
    E_good_DATA(i) = trapz(squeeze(circularData{i}(1,1,:,1)), dataGood{i});
    E_bad_DATA(i) = trapz(squeeze(circularData{i}(1,1,:,1)), dataIntegrated{i} - dataGood{i});
    
end

metric = (E_good - E_bad) ./ (E_good + E_bad);

metricDATA = (E_good_DATA - E_bad_DATA) ./ (E_good_DATA + E_bad_DATA);


% Minimum found at metric(19,28) which corresponds to 73.6 keV and 2.8 mrad
% aperture. 

f5 = figure;
ax = axes('Parent',f5);

% We want to plot 73.0 keV to 74.2 keV
ySTART = 72.8; 
yEND = 74.4;
range = int32(((ySTART - eSTART) / deltaE ) + 1 : ...
    ((yEND - eSTART) / deltaE ) + 1);

ySTART_DATA = 73.0; 
yEND_DATA = 74.2;

plot(ax, ySTART:deltaE:yEND, metric(range,28), ...
    'LineWidth',1.5,'Color','red')
hold on
plot(ax, ySTART_DATA:deltaE*2:yEND_DATA, metricDATA, ...
    'LineStyle','None','Marker','*', 'MarkerEdgeColor','blue',...
    'MarkerSize', 12, 'LineWidth', 1.5)

xlabel(ax, 'E_{CE} (keV)')
ylabel(ax, '( N_{good} - N_{bad} ) / N_{total}')
legend(ax,'FastFit','LCS Code','Location','northwest')
ax.FontSize = 16;
xlim(ax,[ySTART yEND]);
ylim(ax,[0.4, 1.0]);
pbaspect(ax,[1 1 1]);

%% Figure 9. ANNULUS: Find viable aperture radii and associated flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires:
% SECS. 1-3 and SEC. 6 in comptonfastfit.m
%
% NOTE: Running SEC. 6 will overwrite some variables defined in SEC. 5 and
% vice versa.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxAperture = zeros(pp,2);

% Identify apertures that contain metric value over 0.8
metric = (E_good - E_bad) ./ (E_good + E_bad);
flux = (E_good + E_bad);

viable = find(metric > 0.8);
[a, b, c] = ind2sub(size(metric),viable);
viableIndices = zeros(length(viable),3);
for i = 1:size(viable)
    viableIndices(i,:) = [a(i),b(i),c(i)];
end


viableArray = metric > 0.8;
fluxCut = flux .* viableArray;
maxFlux = zeros(size(steps));
maxIndexLinear = zeros(size(steps));
apertures = zeros(size(steps,1),2);

for test = steps(1):steps(end)
    [maxFlux(test),maxIndexLinear(test)] = max(fluxCut(test,:,:),[],'all','linear');
    [~,t2,t3] = ind2sub(size(flux(1,:,:)),maxIndexLinear(test));
    apertures(test,:) = [t2-1,t3-1];
end


f = figure;
ax = axes('Parent',f);
x = eSTART:deltaE:eEND-deltaE;
ax.FontSize = 20;
yyaxis(ax,'right');
xlabel(ax, 'E_{CE}')
plot(ax, x(3:17), apertures(3:17,1).'./10, 'LineWidth',2);
hold on;
plot(ax, x(3:17), apertures(3:17,2).'./10, 'LineWidth',2);
patch([x(3:17) fliplr(x(3:17))], [apertures(3:17,1).'./10 fliplr(apertures(3:17,2).'./10)],[0.9290 0.6940 0.1250],'FaceAlpha',0.5)
ylabel(ax, 'aperture (mrad)')

yyaxis(ax,'left');
plot(ax, x(3:17), maxFlux(3:17), 'LineWidth',2);
ylabel(ax, 'Total X-ray flux (arb. units)')
hold off
pbaspect(ax,[1 1 1])
xlim(ax,[70 80]);
xl1 = xline(ax, 78,'-','no valid apertures', 'LineWidth', 2.5, 'Color', 'black', 'FontSize', 14);
xl2 = xline(ax, 71,'-','no valid apertures', 'LineWidth', 2.5, 'Color', 'black', 'FontSize', 14, 'LabelHorizontalAlignment', 'left');


