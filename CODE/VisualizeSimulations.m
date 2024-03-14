% Script to plot Camino simulation signal results.
% Simulation geometry: square-packed spheres with separation distance
% described by separation factor (sepf) times the sphere radius (r).

% MiJo 20.1.2023. Original version. Written with Matlab R2021a.
% MiJo 3.2.2024. Modified for publishing. 

%% User input
clear all, clc, close all

% Choose data to plot
path = 'My Documents\GitHub\Jokivuolle_MedPhys_2024\DATA';
Gmax = 45; %in units of mT/m. Maximum gradient strength.
b = 2500; %in units of s/mm^2. b-value.
sepf = 1; % Separation factor. Available values in the simulation library: 1, 1.1, 1.7
k = 0; %in units of s^-1. Exchange rate.
D0 = 2; %in units of µm^2/ms. Intrinsic diffusivity.

stype = 3; % Switch: Plot either (1) total signal, (2) intracellular signal, or (3) extracellular signal.
Nreps = 100; % Control the number of diffusing particles: one repetition is a collection of 10 000 diffusing particles. 

% Choose cell sizes to plot: (give in µm)
r_plot = []; % Available values in the simulation library: [2:2:20]. If empty, plot all.

fix_colorscale = 0; %1 = use same colorscale for all plots. 0 = individual colorscales for each plot.

%% Data to load
gr_fn = sprintf('Gmax%d_b%d', Gmax, b);
gr_folder = fullfile(path,'Gradients', gr_fn);
signal_fn = ['sphere_', gr_fn, sprintf('_sepf%d_exch%d_D%d', sepf, k, D0)];
signal_folder = fullfile(path,'Signals',signal_fn);

if isequal(stype, 1)
    main_title = sprintf('Simulated signal in spheres. Total signal.\nSeparation %d. Exchange %d s^{-1}. Gmax = %d mT/m. b = %d s/mm^2. D0 = %d s/mm^2. Nwalkers = %d.', sepf, k, Gmax, b, D0, Nreps*10000); 
elseif isequal(stype, 2)
    main_title = sprintf('Simulated signal in spheres. Intracellular signal.\nSeparation %d. Exchange %d s^{-1}. Gmax = %d mT/m. b = %d s/mm^2. D0 = %d s/mm^2. Nwalkers = %d.', sepf, k, Gmax, b, D0, Nreps*10000); 
elseif isequal(stype, 3)
    main_title = sprintf('Simulated signal in spheres. Extracellular signal.\nSeparation %d. Exchange %d s^{-1}. Gmax = %d mT/m. b = %d s/mm^2. D0 = %d s/mm^2. Nwalkers = %d.', sepf, k, Gmax, b, D0, Nreps*10000); 
else
    disp('Parameter stype needs to be 1 (total signal), 2 (intracellular) or 3 (extracellular).')
end

%% Read data 
load(fullfile(gr_folder, 'Mat-files','big_DELTA.mat'))
load(fullfile(gr_folder, 'Mat-files','small_delta.mat'))
load(fullfile(gr_folder, 'Mat-files','bvalues.mat'))
load(fullfile(signal_folder, 'Mat-files','R'))

Sarray = CreateSarray(stype, signal_folder, Nreps);

%% Pick only certain cell sizes to plot

if isempty(r_plot)
    r_plot = r;
    indices = 1:length(r_plot);
else
    indices = length(r_plot);
    for i = 1:length(r_plot)
        indices(i) = find(r == r_plot(i));
    end
end

Nslices = length(r_plot);
%% Plot
[rows, cols] = NoRowsCols(Nslices);

D_ms = DELTA*1e+3;
d_ms = delta*1e+3;

if isequal(fix_colorscale, 1)
    temp = Sarray(:);
    global_max = max(temp);
    global_min = min(temp(temp>0));
    climits = [global_min, global_max];
else
    climits = [];
end

% Shared plotting preparations for all cell sizes:
% create grid for contour function. Must match the matrix imagesc-function
% visualizes:
[X, Y] = meshgrid([1:length(d_ms)], [1:length(D_ms)]);

% The figure
fig = figure(1); clf
t = tiledlayout(rows, cols, 'TileSpacing','tight','Padding','tight');
for ri = indices
    tile_id = find(indices == ri); % which tile number we are plotting at
    tiletitle = sprintf('%d µm', r(ri));
    ax = WplotTile(t, tile_id, Sarray(:,:,ri), X, Y, d_ms, D_ms, climits, tiletitle);
end

xlabel(t, 'delta (ms)')
ylabel(t, 'DELTA (ms)')
sgtitle(main_title)

set(fig, 'Position', get(0, 'Screensize'))

%% Help functions

function Sarray = CreateSarray(stype, signal_folder, Nreps)
% Read either (1) total signal, (2) intracellular signal, or (3) extracellular signal:
    if isequal(stype, 1)
        Sarray = 0;
        Nwalkers = 0;
        for i = 1:Nreps
            load(fullfile(signal_folder,'Intra_raw',num2str(i),'Sarray_in.mat'))
            load(fullfile(signal_folder, 'Intra_raw',num2str(i),'Nw_in.mat'))
            load(fullfile(signal_folder,'Extra_raw',num2str(i),'Sarray_ex.mat'))
            load(fullfile(signal_folder, 'Extra_raw',num2str(i),'Nw_ex.mat'))
            Sarray = Sarray + Sarray_in + Sarray_ex;
            Nwalkers = Nwalkers + Nw_in + Nw_ex;
            clear Sarray_in Sarray_ex Nw_in Nw_ex 
        end
        Sarray = bsxfun(@rdivide, Sarray, permute(Nwalkers', [1 3 2])); 
        Sarray = permute(Sarray, [2 1 3]);
        load(fullfile(signal_folder, 'Mat-files','sep_factor')) 
        Sarray_mask = Sarray;
    elseif isequal(stype, 2)
        Sarray = 0; %initialise for summing over repetitions
        Nwalkers = 0;
        for i = 1:Nreps
            load(fullfile(signal_folder,'Intra_raw',num2str(i),'Sarray_in.mat'))
            load(fullfile(signal_folder, 'Intra_raw',num2str(i),'Nw_in.mat'))
            Sarray = Sarray + Sarray_in;
            Nwalkers = Nwalkers + Nw_in;
            clear Sarray_in Nw_in
        end
        Sarray = bsxfun(@rdivide, Sarray, permute(Nwalkers', [1 3 2])); % permute: makes Nwalkers into a 3D-vector of 1 row and col. Nwalkers needs to be a ROW vector. bsxfun "broadcasts" the 3D vector so that it becomes size of Sarray. rdivide is element-wise (right) division.
        Sarray = permute(Sarray, [2 1 3]); % Sarray size: N_delta x N_DELTA x N_r -> permute to correspond the cylinder case
        Sarray_mask = Sarray;
    elseif isequal(stype, 3)
        Sarray = 0; %initialise for summing over repetitions
        Nwalkers = 0;
        for i = 1:Nreps
            load(fullfile(signal_folder,'Extra_raw',num2str(i),'Sarray_ex.mat'))
            load(fullfile(signal_folder, 'Extra_raw',num2str(i),'Nw_ex.mat'))
            Sarray = Sarray + Sarray_ex;
            Nwalkers = Nwalkers + Nw_ex;
            clear Sarray_ex Nw_ex
        end
        Sarray = bsxfun(@rdivide, Sarray, permute(Nwalkers', [1 3 2])); % permute: makes Nwalkers into a 3D-vector of 1 row and col. Nwalkers needs to be a ROW vector. bsxfun "broadcasts" the 3D vector so that it becomes size of Sarray. rdivide is element-wise (right) division.
        Sarray = permute(Sarray, [2 1 3]); % Sarray size: N_delta x N_DELTA x N_r -> permute to correspond the cylinder case
        Sarray_mask = Sarray;
    else
        disp('Parameter stype needs to be 1 (total signal), 2 (intracellular) or 3 (extracellular).')
    end

    % Signal equal to 1 is likely to originate from unfeasible gradient
    % combinations (i.e. simulations without a gradient).
    Sarray(Sarray_mask >= 1) =  0;
end

function [rows, cols] = NoRowsCols(Nplots)
% Function to determine adequate number of rows and cols for a
% tiled/subplot layout based on the number of required plots

% MiJo 29.5.2023

    if Nplots > 20
        cols = 6;
        rows = ceil(Nplots/cols);
    elseif 10 <= Nplots && Nplots <= 20
        cols = 5;
        rows = ceil(Nplots/cols);
    elseif 6 <= Nplots && Nplots < 10
        cols = 5;
        rows = ceil(Nplots/cols);
    elseif 1 < Nplots && Nplots <= 5
        cols = Nplots;
        rows = 1;
    elseif Nplots == 1
        cols = 1;
        rows = 1;
    end
end

function [ax1, ax2] = WplotTile(t, i, array2D, X, Y, d_ms, D_ms, climits, tiletitle)

% For visualisation purposes, find the min and max signal:
    % find the minimum and maximum signal points:
    s_max = max(array2D(:));
    s_min = min(array2D(array2D > 0)); % find nonzero minimum
    index_max = find(array2D == s_max, 1); % N.B.! Sometimes there are several indeces. Now you only find the first one.
    index_min = find(array2D == s_min, 1);
    % Calculate the contrast (signal difference) per cell size:
    diff = s_max - s_min;

% PLOT
    ax1 = nexttile(t);
    s = imagesc(array2D); %show 2D data
    hold on

% Modify axes
    ax1.YDir = 'normal'; %imagesc plot 'reversed', need to set it back
    % You need to correct for the mismatch between the actual delta & DELTA
    % values and the matrix indexing in Matlab:
    xticks(linspace(1, length(d_ms), 4));
    yticks(linspace(1, length(D_ms), 7));
    xticklabels(num2cell(linspace(min(d_ms), max(d_ms), 4)));
    yticklabels(num2cell(linspace(min(D_ms), max(D_ms), 7)));

% Colormap/caxis modifications
    cmap1 = colormap(ax1, 'parula');
    c1 = colorbar(ax1);
    if isempty(climits)
        caxis([s_min s_max])
    else
        caxis(climits)
    end
    % Add isocontours. First, create a mask of the area containing values
    % to prevent contours appear at the 'sides' of the triangular plane
    binary = array2D;
    binary(binary > 0) = 1;
    edg = edge(binary);
    temp = array2D;
    temp(edg) = NaN;
    levels = c1.Ticks;
    [Con, h] = contour(X, Y, temp, levels, 'k');
    clabel(Con, h, 'LabelSpacing', 200) 
    colorbar('off');

% Mark certain parameter pairs. Using plot3 as the index is the linear
% index for 2D matrix -> works for X and Y
    plot3(X(index_max), Y(index_max), s_max, 'k.', 'MarkerSize', 20);
    plot3(X(index_min), Y(index_min), s_max, 'k.', 'MarkerSize', 20);
    text(X(index_max+5), Y(index_max-5), s_max, sprintf('%.3f', s_max), 'FontWeight', 'bold')
    text(X(index_min+5), Y(index_min-5), s_max, sprintf('%.3f', s_min), 'FontWeight', 'bold')

% Title & final edits to axes
    title(tiletitle)
    axis equal tight, grid on

% Contrast scale:
    ax2 = axes(t);
    ax2.Layout.Tile = i;
    ax2.Visible = 'off';
    cmap2 = colormap(ax2, white);
    c2 = colorbar(ax2);
    % specifying ticklabels
    c2.Label.String = 'Contrast';
    caxis([0 1])
    c2.Ticks = [0, diff, 1];
    ticklabels = {'0',sprintf('%.2f', diff),'1'};
    c2.TickLabels = ticklabels;
    
end % function