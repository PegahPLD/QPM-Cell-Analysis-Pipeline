%% Step 1: QPM Phase Unwrapping and Thickness Calculation
clear all; close all;

% Parameters
imagenum = 10;               % Total number of QPM image stacks
dish = 221;                  % Dish/sample ID
nframe = 30;                 % Number of phase-shifting frames per image
paraboloid_coeff = -19;      % Paraboloid correction coefficient
lambda = 0;                  % SPUD threshold (not used here)

% Loop through all image stacks
for i = 1:imagenum
    fileorig = strcat(num2str(i), '.tif');
    
    % Load and stack all frames
    imtemp = imread(fileorig, 1);
    [xSize, ySize] = size(imtemp);
    inputImg = zeros(xSize, ySize, nframe);
    
    for ii = 1:nframe
        inputImg(:,:,ii) = double(imread(fileorig, ii));
    end
    
    % Compute complex wavefront from fringes
    comp = fringes2phase(inputImg);
    
    % Paraboloid correction
    [xMat, yMat] = meshgrid(1:ySize, 1:xSize);
    xMat = xMat - ySize / 2;
    yMat = yMat - xSize / 2;
    comp_correct = comp .* exp(1i * paraboloid_coeff * (xMat.^2 + yMat.^2) / 1e6);

    % Oversample and shift Fourier domain
    comp_small = comp_correct(1:4:end, 1:4:end);
    comp_oversample = zeros(xSize, ySize);
    x_start = round(xSize * 3/8);
    y_start = round(ySize * 3/8);
    comp_oversample(x_start:x_start+size(comp_small,1)-1, y_start:y_start+size(comp_small,2)-1) = comp_small;
    
    fft3 = abs(fftshift(fft2(comp_oversample)));
    [xtmp, ytmp] = find(fft3 == max(fft3(:)));
    
    xDiff = (xtmp - (xSize / 2 + 1)) / 4;
    yDiff = (ytmp - (ySize / 2 + 1)) / 4;

    % Apply shift correction
    comp_correct2 = comp_correct .* exp(-1i * 2 * pi * (xDiff * xMat / xSize + yDiff * yMat / ySize));
    comp_correct2 = comp_correct2 ./ sum(comp_correct2(:));

    % Get wrapped phase
    psi = angle(comp_correct2);
    psi = 2 * pi ./ (max(psi(:)) - min(psi(:))) .* (psi - min(psi(:))) - pi;
    imagesc(psi); title('Wrapped Phase');

    %% Phase Unwrapping using Goldstein Algorithm
    im_mag = abs(comp_correct2);        % Magnitude image
    im_phase = angle(comp_correct2);    % Phase image
    im_mask = ones(size(comp_correct2));% Full mask
    
    % Compute residues and branch cuts
    residue_charge = PhaseResidues_r1(im_phase, im_mask);
    max_box_radius = floor(min(size(im_phase)) / 2);
    branch_cuts = BranchCuts_r1(residue_charge, max_box_radius, im_mask);
    
    im_mask(branch_cuts) = 0;
    im_mag1 = im_mag .* im_mask;
    
    % Select starting point manually
    im_phase_quality = im_mag1;
    figure; imagesc(im_phase_quality); colormap(gray); title('Select Reference Point');
    [xpoint, ypoint] = ginput(1);
    colref = round(xpoint);
    rowref = round(ypoint);
    close;

    % Flood-fill phase unwrapping
    unwrapped = FloodFill_r1(im_phase, im_mag, branch_cuts, im_mask, colref, rowref);

    % Replace NaNs with surrounding means
    nan_indices = find(isnan(unwrapped));
    for k = 1:length(nan_indices)
        [row, col] = ind2sub(size(unwrapped), nan_indices(k));
        neighbors = unwrapped(max(1,row-1):min(end,row+1), max(1,col-1):min(end,col+1));
        unwrapped(row, col) = mean(neighbors(~isnan(neighbors)));
    end

    im_unwrapped(:,:,i) = unwrapped;

    %% Thickness & OPD Calculations
    % result: thickness in microns
    resultss = unwrapped * 633 / (4 * pi * 0.033 * 1000); % 633nm wavelength, 0.033 refractive index diff
    results(:,:,i) = resultss;

    % OPD in meters
    OPD(:,:,i) = resultss * 0.033 * 1e-6;
    Phase(:,:,i) = unwrapped;

    % Save images
    figure('WindowState','maximized'); 
    imagesc(resultss); colormap(jet); colorbar; axis square; caxis([-1 30]); 
    title('Thickness (μm)');
    saveas(gcf, ['thickness_' num2str(i) '.png']);

    % Normalize and save grayscale
    I2 = uint8(255 * mat2gray(resultss));
    imwrite(I2, ['myImage_' num2str(i) '.tif']);
    close all;
end

% Save all results
save(sprintf('results%d.mat', dish), 'results', '-mat');
save(sprintf('phase%d.mat', dish), 'im_unwrapped', '-mat');
save(sprintf('OPD%d.mat', dish), 'OPD', '-mat');
