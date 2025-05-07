%% Cell Process Step
clear all; close all;

% Load previous outputs
load orig3Dim221.mat
load ImcropMatrix221.mat

% Constants
dish = 221;
pixel_size = 0.35;
alpha = 0.18;

% Feature extraction loop
for k = 1:size(Imcropmatrix,1)
    for im = 1:size(Imcropmatrix,3)
        my_im = sprintf('s%d_%d', im, k);
        mask1 = imread(['grayim' num2str(im) '_' num2str(k) '.tif']);
        mask1 = imbinarize(mask1);

        % Background correction
        meanbg1(k) = mean(orig3Dim.(my_im)(mask1 == 0));
        bgc3D.(my_im) = orig3Dim.(my_im) - meanbg1(k);

        % Region with largest area
        stats1 = regionprops3(mask1, orig3Dim.(my_im), 'SurfaceArea');
        [~, idx] = max(stats1.SurfaceArea);
        cc1 = bwconncomp(mask1);
        BW21 = ismember(labelmatrix(cc1), idx);
        finalmask.(my_im) = BW21;
        imwrite(BW21, ['finalmaskim' num2str(im) '_' num2str(k) '.tif']);

        % Background-subtracted image
        Z = zeros(size(bgc3D.(my_im)));
        Z(finalmask.(my_im) == 1) = bgc3D.(my_im)(finalmask.(my_im) == 1);

        % Feature calculations
        volume(k,im) = sum(Z(:)) * pixel_size^2;
        Areaa(k,im) = sum(finalmask.(my_im)(:)) * pixel_size^2;
        voldivarea(k,im) = volume(k,im) / Areaa(k,im);

        props = regionprops(finalmask.(my_im), Z, ...
            'Area', 'Centroid','WeightedCentroid','MaxFeretProperties','MinFeretProperties','Circularity');
        
        flatness(k,im) = props.MaxFeretDiameter / props.MinFeretDiameter;
        circ(k,im) = props.Circularity;

        drymass(k,im) = (1 / alpha) * sum(Z(:)) * pixel_size^2 / 1000;
        [delta_hx, delta_hy] = gradient(double(Z));
        integrand = sqrt(1 + (delta_hx * pixel_size).^2 + (delta_hy * pixel_size).^2);
        SA(k,im) = sum(integrand(:)) * pixel_size^2;

        roundness(k,im) = (4 * pi * volume(k,im)) / SA(k,im)^2;
        SAV(k,im) = SA(k,im) / volume(k,im);
        SDM(k,im) = SA(k,im) / drymass(k,im);
        PAV(k,im) = Areaa(k,im) / volume(k,im);

        nc = nnz(finalmask.(my_im));
        Hmean(k,im) = mean(Z(finalmask.(my_im) == 1));
        Spherecity(k,im) = pi^(1/3)*(6*Hmean(k,im))^(2/3)/SA(k,im);
        HV(k,im) = var(Z(finalmask.(my_im) == 1));
        KT(k,im) = kurtosis(Z(finalmask.(my_im) == 1));
        SK(k,im) = skewness(Z(finalmask.(my_im) == 1));
        Ecc(k,im) = (props.MaxFeretDiameter - props.MinFeretDiameter) / ...
                    (props.MaxFeretDiameter + props.MinFeretDiameter);
    end
end

% Save all features
variables = {'finalmask','bgc3D','circ','voldivarea','flatness','roundness','volume','Areaa','drymass',...
             'SA','SAV','SDM','PAV','Spherecity','HV','KT','SK','Ecc'};
for v = variables
    save(sprintf('%s%d.mat', v{1}, dish), v{1}, '-mat')
end

%% Feature Visualization and Statistical Analysis

% Define features and units
features = {'SA', 'SAV', 'HV', 'KT', 'SK', 'drymass', 'PAR', 'SDM', 'Areaa', 'volume', ...
            'circ','flatness','roundness', 'voldivarea','Ecc'};
Units = {'μm²', '1/μm', 'μm²', '1/μm⁴', '1/μm³', 'pg', '1/μm', 'μm/pg', ...
         'μm²', 'μm³', '', '', '', 'μm', ''};

% ANOVA & boxplot
for i = 1:length(features)
    feature = features{i};
    aov = anova(eval(feature));
    figure('WindowState','maximized'); boxchart(aov);
    ylabel([feature ' [' Units{i} ']']);
    fontsize(gca,25,'pixels');
    saveas(gcf, [feature '.bmp']);
end

%% Sign Rank Test
for i = 1:size(Imcropmatrix,3)
    for j = 1:length(features)
        feature = features{j};
        dataT = eval(feature);
        [p, h, stats] = signrank(dataT(:,2), dataT(:,i));

        fprintf('p-value (%s): %.4f\n', feature, p);
        if h
            fprintf('Significant difference in %s between img2 and img%d\n', feature, i);
        else
            fprintf('No significant difference in %s between img2 and img%d\n', feature, i);
        end
        fprintf('Test statistic: %d\n\n', stats.signedrank);
    end
end
