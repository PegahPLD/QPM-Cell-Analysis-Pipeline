%% Step 5: Feature Aggregation and Statistical Analysis
close all;
clear all;

% Total number of data sets (dishes)
datanum = 222;

% Initialize all feature arrays
circT = []; flatnessT = []; celldistT = []; roundnessT = []; voldivareaT = [];
volumeT = []; PAT = []; drymassT = []; SAT = []; SAVT = []; HVT = [];
SKT = []; KTT = []; PAVT = []; SDMT = []; SphericityT = []; EccT = [];
DishNT = []; perimT = []; perimdivareaT = []; complexityscoreT = [];

% Loop through all datasets
for i = 12:datanum
    % Only process if all expected files exist
    if isfile(sprintf('circ%d.mat', i))
        load(sprintf('circ%d.mat', i));                circT = [circ(:, 1:10); circT];
        load(sprintf('flatness%d.mat', i));            flatnessT = [flatness(:, 1:10); flatnessT];
        load(sprintf('celldist%d.mat', i));            celldistT = [celldist'; celldistT];
        load(sprintf('roundness%d.mat', i));           roundnessT = [roundness(:, 1:10); roundnessT];
        load(sprintf('voldivarea%d.mat', i));          voldivareaT = [voldivarea(:, 1:10); voldivareaT];
        load(sprintf('volume%d.mat', i));              volumeT = [volume(:, 1:10); volumeT];
        load(sprintf('PA%d.mat', i));                  PAT = [PA(:, 1:10); PAT];
        load(sprintf('drymass%d.mat', i));             drymassT = [drymass(:, 1:10); drymassT];
        load(sprintf('SA%d.mat', i));                  SAT = [SA(:, 1:10); SAT];
        load(sprintf('SAV%d.mat', i));                 SAVT = [SAV(:, 1:10); SAVT];
        load(sprintf('HV%d.mat', i));                  HVT = [HV(:, 1:10); HVT];
        load(sprintf('SK%d.mat', i));                  SKT = [SK(:, 1:10); SKT];
        load(sprintf('KT%d.mat', i));                  KTT = [KT(:, 1:10); KTT];
        load(sprintf('PAV%d.mat', i));                 PAVT = [PAV(:, 1:10); PAVT];
        load(sprintf('SDM%d.mat', i));                 SDMT = [SDM(:, 1:10); SDMT];
        load(sprintf('Spherecity%d.mat', i));          SphericityT = [Spherecity(:, 1:10); SphericityT];
        load(sprintf('Ecc%d.mat', i));                 EccT = [Ecc(:, 1:10); EccT];
        load(sprintf('perim%d.mat', i));               perimT = [perim(:, 1:10); perimT];
        load(sprintf('perimdivarea%d.mat', i));        perimdivareaT = [perimdivarea(:, 1:10); perimdivareaT];
        load(sprintf('complexityscore%d.mat', i));     complexityscoreT = [complexityScore(:, 1:10); complexityscoreT];

        DishN = ones(size(celldist)) * i;
        DishNT = [DishN'; DishNT];
    end
end

%% Grouping Data Based on Distance Threshold
threshold = 300;

% Separate features: Above threshold (far from LIS)
volumeTu       = volumeT(celldistT >= threshold, :);
voldivareaTu   = voldivareaT(celldistT >= threshold, :);
roundnessTu    = roundnessT(celldistT >= threshold, :);
flatnessTu     = flatnessT(celldistT >= threshold, :);
circTu         = circT(celldistT >= threshold, :);
SATu           = SAT(celldistT >= threshold, :);
SAVTu          = SAVT(celldistT >= threshold, :);
HVTu           = HVT(celldistT >= threshold, :);
KTTu           = KTT(celldistT >= threshold, :);
SKTu           = SKT(celldistT >= threshold, :);
drymassTu      = drymassT(celldistT >= threshold, :);
PAVTu          = PAVT(celldistT >= threshold, :);
SDMTu          = SDMT(celldistT >= threshold, :);
SphericityTu   = SphericityT(celldistT >= threshold, :);
AreaTu         = PAT(celldistT >= threshold, :);
EccTu          = EccT(celldistT >= threshold, :);
perimTu        = perimT(celldistT >= threshold, :);
perimdivareaTu = perimdivareaT(celldistT >= threshold, :);

% Below threshold (close to LIS)
volumeTb       = volumeT(celldistT < threshold, :);
voldivareaTb   = voldivareaT(celldistT < threshold, :);
roundnessTb    = roundnessT(celldistT < threshold, :);
flatnessTb     = flatnessT(celldistT < threshold, :);
circTb         = circT(celldistT < threshold, :);
SATb           = SAT(celldistT < threshold, :);
SAVTb          = SAVT(celldistT < threshold, :);
HVTb           = HVT(celldistT < threshold, :);
KTTb           = KTT(celldistT < threshold, :);
SKTb           = SKT(celldistT < threshold, :);
drymassTb      = drymassT(celldistT < threshold, :);
PAVTb          = PAVT(celldistT < threshold, :);
SDMTb          = SDMT(celldistT < threshold, :);
SphericityTb   = SphericityT(celldistT < threshold, :);
AreaTb         = PAT(celldistT < threshold, :);
EccTb          = EccT(celldistT < threshold, :);
perimdivareaTb = perimdivareaT(celldistT < threshold, :);

%% Create Summary Table
Totaldata = vertcat(SAT', SAVT', HVT', KTT', SKT', drymassT', PAVT', SDMT', PAT', volumeT', circT', ...
                    flatnessT', roundnessT', voldivareaT', perimT', perimdivareaT', celldistT', DishNT');

T = table(SAT, SAVT, HVT, KTT, SKT, drymassT, PAVT, SDMT, PAT, volumeT, circT, ...
          flatnessT, roundnessT, voldivareaT, perimT, celldistT, perimdivareaT, DishNT, ...
          'VariableNames', {'SAT', 'SAVT', 'HVT', 'KTT', 'SKT', 'drymassT', 'PAVT', 'SDMT', ...
                            'PAT', 'volumeT', 'circT', 'flatnessT', 'roundnessT', 'voldivareaT', ...
                            'perimT', 'celldistT', 'perimdivareaT', 'DishNT'});

%% Sign-Rank Statistical Tests
features = {'circ', 'voldivarea', 'flatness', 'roundness', 'volume', 'PA', 'drymass', ...
            'SA', 'SAV', 'SDM', 'PAV', 'Spherecity', 'HV', 'KT', 'SK', 'Ecc', 'perim', 'perimdivarea'};

% Compare: 1 min before vs 10 min after LIS
disp('--- 1 min before vs 10 min after LIS ---');
for featureIdx = 1:length(features)
    feature = features{featureIdx};
    dataT = eval([feature 'T']);
    [p, h, stats] = signrank(dataT(:,2), dataT(:,10));
    fprintf('p-value (%s): %.4f\n', feature, p);
    if h
        fprintf('Significant difference detected.\n');
    else
        fprintf('No significant difference detected.\n');
    end
    fprintf('Test statistic: %d\n\n', stats.signedrank);
end

% Compare: 10 min after vs 2 hr after LIS
disp('--- 10 min vs 2 hr after LIS ---');
for featureIdx = 1:length(features)
    feature = features{featureIdx};
    dataT = eval([feature 'T']);
    [p, h, stats] = signrank(dataT(:,10), dataT(:,9));
    fprintf('p-value (%s): %.4f\n', feature, p);
    if h
        fprintf('Significant difference detected.\n');
    else
        fprintf('No significant difference detected.\n');
    end
    fprintf('Test statistic: %d\n\n', stats.signedrank);
end

% Compare: 1 min before vs 2 hr after LIS
disp('--- 1 min before vs 2 hr after LIS ---');
for featureIdx = 1:length(features)
    feature = features{featureIdx};
    dataT = eval([feature 'T']);
    [p, h, stats] = signrank(dataT(:,2), dataT(:,9));
    fprintf('p-value (%s): %.4f\n', feature, p);
    if h
        fprintf('Significant difference detected.\n');
    else
        fprintf('No significant difference detected.\n');
    end
    fprintf('Test statistic: %d\n\n', stats.signedrank);
end

%% ANOVA Plot Generation
Units = { 'μm²', '1/μm',  'μm²', '1/μm⁴', '1/μm³', 'pg', '1/μm', 'μm²/pg', 'μm²', ...
          'μm³', '1', '1', '1', 'μm', 'μm', '1/μm', '1', '1' };

% Plot full data per feature
for featureIdx = 1:length(features)
    feature = features{featureIdx};
    aov = anova(eval([feature 'T']));
    figure('WindowState', 'maximized'); boxchart(aov);
    fontsize(gca, 25, 'pixels');
    ylabel([feature ' [' Units{featureIdx} ']']);
    saveas(gcf, [feature 'T.bmp'], 'bmp');
end

% Grouped boxplots by distance
for suffix = ["Tu", "Tb"]
    label = ternary(suffix == "Tu", 'Distances > 300 µm', 'Distances < 300 µm');
    for featureIdx = 1:length(features)
        feature = features{featureIdx};
        aov = anova(eval([feature + suffix]));
        figure('WindowState', 'maximized'); boxchart(aov);
        fontsize(gca, 25, 'pixels');
        xlabel(label);
        ylabel([feature ' [' Units{featureIdx} ']']);
        saveas(gcf, [feature + suffix + '.bmp'], 'bmp');
    end
end

%% Selected Comparison Plot (Columns 2, 10, 9)
columns_to_plot = [2, 10, 9];

for featureIdx = 1:length(features)
    feature = features{featureIdx};
    dataT = eval([feature 'T']);
    selected_data = dataT(:, columns_to_plot);
    aov = anova(selected_data);
    figure('WindowState', 'maximized'); boxchart(aov);
    fontsize(gca, 25, 'pixels');
    ylabel([feature ' [' Units{featureIdx} ']']);
    saveas(gcf, [feature 'T_Columns_2_10_9.bmp'], 'bmp');
end

%% Create Binary Significance Table
significance = zeros(length(features), 8);
for i = 3:10
    for featureIdx = 1:length(features)
        feature = features{featureIdx};
        dataT = eval([feature 'T']);
        [~, h] = signrank(dataT(:,2), dataT(:,i));
        significance(featureIdx, i-2) = h;
    end
end
