%% Cell Crop Step
clear all; close all;

% Load pre-processed result image and center image
load results.mat
LIScent = imread("LIScent.tiff");
save('LIScent.mat','LIScent','-mat')

% Set scale and dish identifier
pixelsize = 0.35;  % micrometers per pixel
dish = 221;

%% Load and parse live cell data
filename = 'Results221.csv';
data = readtable(filename);
sliceNumbers = unique(data.Slice);

% Preallocate matrix for cropped live cells
Imcropmatrix = zeros(height(data)/numel(sliceNumbers), 4, numel(sliceNumbers));

% Organize bounding boxes per slice
for i = 1:numel(sliceNumbers)
    sliceData = data(data.Slice == sliceNumbers(i), {'BX', 'BY', 'Width', 'Height'});
    Imcropmatrix(1:height(sliceData), :, i) = table2array(sliceData);
end

%% Load and parse dead cell data
filenamed = 'dead221.csv';
datad = readtable(filenamed);
sliceNumbersd = unique(datad.Slice);

deadcropmatrix = zeros(height(datad)/numel(sliceNumbersd), 4, numel(sliceNumbersd));
for i = 1:numel(sliceNumbersd)
    sliceDatad = datad(datad.Slice == sliceNumbersd(i), {'BX', 'BY', 'Width', 'Height'});
    deadcropmatrix(1:height(sliceDatad), :, i) = table2array(sliceDatad);
end

% Visualize and capture user input point
imshow(LIScent)
[xpoint, ypoint] = ginput(1); 

% Initialize distance matrices
weighteddistance = zeros(size(Imcropmatrix,1), size(Imcropmatrix,3));
celldist = zeros(size(Imcropmatrix,1), 1);

% Process each image slice
for k = 1:size(Imcropmatrix,1)
    for im = 1:numel(sliceNumbers)
        % Compute distance from click point to cell center
        celldist(k) = sqrt((Imcropmatrix(k,1,im)+Imcropmatrix(k,3,im)/2 - xpoint)^2 + ...
                           (Imcropmatrix(k,2,im)+Imcropmatrix(k,4,im)/2 - ypoint)^2) * pixelsize;

        % Crop original image region
        my_im = sprintf('s%d_%d', im, k);
        D = results(:,:,im);
        orig3Dim.(my_im) = imcrop(D, Imcropmatrix(k,:,im));
        
        % Normalize and save
        I = uint8(255 * mat2gray(orig3Dim.(my_im)));
        imwrite(I, ['grayim' num2str(im) '_' num2str(k) '.tif'])

        % Process distances to dead cells
        for m = 1:size(deadcropmatrix,1)
            dist = sqrt((Imcropmatrix(k,1,im)+Imcropmatrix(k,3,im)/2 - ...
                         (deadcropmatrix(m,1,im)+deadcropmatrix(m,3,im)/2))^2 + ...
                        (Imcropmatrix(k,2,im)+Imcropmatrix(k,4,im)/2 - ...
                         (deadcropmatrix(m,2,im)+deadcropmatrix(m,4,im)/2))^2) * pixelsize;

            weighteddistance(k,im) = weighteddistance(k,im) + 1 / dist^2;

            % Crop and save dead region
            dead_im = sprintf('s%d_%d', im, m);
            dead3D.(dead_im) = imcrop(D, deadcropmatrix(m,:,im));
            Id = uint8(255 * mat2gray(dead3D.(dead_im)));
            imwrite(Id, ['deadim' num2str(im) '_' num2str(m) '.tif']);
        end
    end
end

% Save all output matrices
save(sprintf('orig3Dim%d.mat', dish), 'orig3Dim', '-mat');
save(sprintf('Imcropmatrix%d.mat', dish), 'Imcropmatrix', '-mat');
save(sprintf('celldist%d.mat', dish), 'celldist', '-mat');
save(sprintf('deadcropmatrix%d.mat', dish), 'deadcropmatrix', '-mat');
