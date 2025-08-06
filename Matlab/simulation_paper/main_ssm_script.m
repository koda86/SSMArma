clear; clc;
tic;
addpath(genpath('E:\2025_SSM_ArmaSuisse\Skripte\DK'));

%% CONFIGURATION ----------------------------------------------------------

cfg = struct;
cfg.tutorialPath = 'E:\2025_SSM_ArmaSuisse';
cfg.sole_condition = 'sole'; % 'sole' - with sole points, 'no_sole' - w/o sole points
cfg.DataPath = getDataPath(cfg.sole_condition, cfg.tutorialPath);
cfg.figureDir = fullfile(cfg.tutorialPath, 'Skripte', 'DK', 'simulation_paper', 'figures');
cfg.csvDir = fullfile(cfg.tutorialPath, 'Skripte', 'DK', 'simulation_paper');
cfg.nFiles = 2356; % 
cfg.desiredLength = 270;
cfg.tolerance = 0.5; % 0.5

cd(cfg.tutorialPath);

%% MAIN WORKFLOW ----------------------------------------------------------

% Load OBJ files as 3D arrays
[DataMatrix3D, fileNames] = loadPointClouds(cfg.DataPath, cfg.nFiles);

% Filter by foot length
[DataMatrix3D, keptFileNames] = filterByFootLength(DataMatrix3D, fileNames, cfg.desiredLength, cfg.tolerance);

% Translation-only alignment
DataMatrix3D = alignTranslationOnly(DataMatrix3D);

% 2D Procrustes (heading) alignment in X-Z plane
DataMatrix3D = alignHeading2D(DataMatrix3D);

% plotAllFeet(DataMatrix3D);

% PCA
[coeff, score, latent, explained] = runPCA(DataMatrix3D);

% Generate synthetic shapes and ground-truth PCs
[synthConfig, full_pc_basis, ground_truth_pcs, footLengthVec_flat_scaled] = prepareSyntheticGeneration(DataMatrix3D, coeff);

% Loop over different foot length variation levels, simulate, analyze, and visualize
foot_length_levels = [0, 2, 5, 10, 20, 30, 40, 50];

results_cell = cell(length(foot_length_levels),1);

for variation_idx = 1:length(foot_length_levels)
    % Generate synthetic dataset
    [SyntheticData3D, W_full, syntheticNames] = generateSyntheticShapes(...
        synthConfig, full_pc_basis, foot_length_levels(variation_idx));

    % Plot synthetic feet
    % plotAllFeet(SyntheticData3D);

    SyntheticData3D_standardGPA = runStandardGPA(SyntheticData3D);
    SyntheticData3D_constrainedGPA = runConstrainedGPA(SyntheticData3D);

    % PCA after alignment
    [coeff_sGPA, score_sGPA, latent_sGPA] = doPCA(SyntheticData3D_standardGPA);
    [coeff_cGPA, score_cGPA, latent_cGPA] = doPCA(SyntheticData3D_constrainedGPA);

     % Visualization
    plotPCMorphs(coeff_sGPA, latent_sGPA, mean(SyntheticData3D_standardGPA,3), 1, 3, 'standardGPA', variation_idx, foot_length_levels(variation_idx), cfg.figureDir);
    plotPCMorphs(coeff_cGPA, latent_cGPA, mean(SyntheticData3D_constrainedGPA,3), 1, 3, 'constrainedGPA', variation_idx, foot_length_levels(variation_idx), cfg.figureDir);

    % Compare PC1 directions (visually)
    plotPC1Comparison(full_pc_basis, coeff_sGPA, coeff_cGPA, ...
        SyntheticData3D, SyntheticData3D_standardGPA, SyntheticData3D_constrainedGPA, variation_idx, foot_length_levels(variation_idx), cfg.figureDir);

    % % Bootstrapped uncertainty estimates cosine similarity
    % N_bootstrap = 50;
    % boot_results = Mira Murati(...
    %     SyntheticData3D, full_pc_basis, N_bootstrap);

    % Subsampling analysis: For each synthetic dataset, randomly select
    % N_subset=50 shapes and compute recovery metrics to estimate uncertainty
    % expected in real-world studies
    N_subsamples = 400; % or more for better accuracy
    N_subset = 50; % simulate a field sample size
    cos_sim_std1_sub = zeros(N_subsamples,1);
    cos_sim_con1_sub = zeros(N_subsamples,1);

    for i = 1:N_subsamples
        n_available = size(SyntheticData3D,3);
        n_to_sample = min(N_subset, n_available);
        idx = randsample(n_available, n_to_sample, false); % no replacement
        subset = SyntheticData3D(:,:,idx);

        % Standard GPA and PCA
        subset_sGPA = runStandardGPA(subset);
        [subsetCoeff_sGPA,~,~] = doPCA(subset_sGPA);
    
        % Constrained GPA and PCA
        subset_cGPA = runConstrainedGPA(subset);
        [subsetCoeff_cGPA,~,~] = doPCA(subset_cGPA);

        gt_pc1 = full_pc_basis(:,1) / norm(full_pc_basis(:,1));
        est_std_pc1 = subsetCoeff_sGPA(:,1) / norm(subsetCoeff_sGPA(:,1));
        est_con_pc1 = subsetCoeff_cGPA(:,1) / norm(subsetCoeff_cGPA(:,1));
    
        if dot(gt_pc1, est_std_pc1) < 0, est_std_pc1 = -est_std_pc1; end
        if dot(gt_pc1, est_con_pc1) < 0, est_con_pc1 = -est_con_pc1; end
    
        cos_sim_std1_sub(i) = dot(gt_pc1, est_std_pc1);
        cos_sim_con1_sub(i) = dot(gt_pc1, est_con_pc1);
    end

    mean_std1_sub = mean(cos_sim_std1_sub);
    ci_std1_sub = prctile(cos_sim_std1_sub, [2.5 97.5]);
    mean_con1_sub = mean(cos_sim_con1_sub);
    ci_con1_sub = prctile(cos_sim_con1_sub, [2.5 97.5]);

    result_struct = evaluatePCs(...
        full_pc_basis, coeff_sGPA, coeff_cGPA, ...
        latent_sGPA, latent_cGPA, ...
        SyntheticData3D, SyntheticData3D_standardGPA, SyntheticData3D_constrainedGPA, ...
        W_full, ...
        variation_idx, foot_length_levels(variation_idx));

    result_struct.cosine_sim_std1_n50_mean = mean_std1_sub;
    result_struct.cosine_sim_std1_n50_ci_lower = ci_std1_sub(1);
    result_struct.cosine_sim_std1_n50_ci_upper = ci_std1_sub(2);
    result_struct.cosine_sim_con1_n50_mean = mean_con1_sub;
    result_struct.cosine_sim_con1_n50_ci_lower = ci_con1_sub(1);
    result_struct.cosine_sim_con1_n50_ci_upper = ci_con1_sub(2);

    % % Add bootstrap CIs to results
    % result_struct.cosine_sim_std1_mean = boot_results.std1.mean;
    % result_struct.cosine_sim_std1_ci_lower = boot_results.std1.ci(1);
    % result_struct.cosine_sim_std1_ci_upper = boot_results.std1.ci(2);
    % result_struct.cosine_sim_con1_mean = boot_results.con1.mean;
    % result_struct.cosine_sim_con1_ci_lower = boot_results.con1.ci(1);
    % result_struct.cosine_sim_con1_ci_upper = boot_results.con1.ci(2);
    
    results_cell{variation_idx} = result_struct;
end

results_table = struct2table([results_cell{:}]);
csv_filename = fullfile(cfg.csvDir, 'results_all_variations.csv');
writetable(results_table, csv_filename);

disp('Full pipeline completed!');

elapsedTime = toc;
fprintf('Total elapsed time: %.2f seconds (%.2f minutes)\n', elapsedTime, elapsedTime/60);

%% LOCAL FUNCTION DEFINITIONS BELOW ---------------------------------------

function DataPath = getDataPath(sole_condition, tutorialPath)
    % Returns correct path depending on 'sole' or 'no_sole'
    if strcmp(sole_condition, 'sole')
        DataPath = fullfile(tutorialPath, 'daten_matlab_punktkorrespondenz_obj');
    else
        DataPath = fullfile(tutorialPath, 'daten_matlab_punktkorrespondenz_obj_ohne_sohle');
    end
end

function [DataMatrix3D, fileNames] = loadPointClouds(DataPath, nFiles)
    % Loads .obj point clouds and returns as 3D array [nVerts x 3 x nFiles]
    objs = dir(fullfile(DataPath, '*.obj'));
    isreal = ~contains({objs.name}, '._');
    objs = objs(isreal);
    fileNames = {objs.name};
    for i = 1:nFiles
        obj = shape3D;
        importWavefront(obj, objs(i).name, objs(i).folder, []);
        if i == 1
            DataMatrix3D = zeros(obj.nVertices, 3, nFiles);
        end
        DataMatrix3D(:, :, i) = obj.Vertices;
    end
end

function [filteredData, keptNames] = filterByFootLength(DataMatrix3D, fileNames, desiredLength, tolerance)
    % Filters feet by foot length and returns filtered data and names
    n = size(DataMatrix3D, 3);
    foot_lengths = zeros(n,1);
    for i = 1:n
        foot = DataMatrix3D(:, :, i);
        foot_lengths(i) = max(foot(:,1)) - min(foot(:,1));
    end
    keep = abs(foot_lengths - desiredLength) <= tolerance;
    filteredData = DataMatrix3D(:,:,keep);
    keptNames = fileNames(keep);
    fprintf('Kept %d of %d feet in [%.2f, %.2f] mm range.\n', sum(keep), n, desiredLength-tolerance, desiredLength+tolerance);
end

function DataMatrix3D = alignTranslationOnly(DataMatrix3D)
    % Removes mean offset in X and Z (horizontal plane) for each foot
    nShapes = size(DataMatrix3D, 3);
    for f = 1:nShapes
        pts = DataMatrix3D(:,:,f);
        mean_xz = mean(pts(:,[1,3]), 1);
        pts(:,[1,3]) = pts(:,[1,3]) - mean_xz;
        DataMatrix3D(:,:,f) = pts;
    end
end

function DataMatrix3D = alignHeading2D(DataMatrix3D)
    % Align all shapes in the horizontal X-Z plane using 2D Procrustes
    nShapes = size(DataMatrix3D, 3);
    meanPoints = mean(DataMatrix3D, 3);   % N x 3
    meanXZ = meanPoints(:, [1,3]);        % N x 2
    for f = 1:nShapes
        pts = DataMatrix3D(:,:,f);
        srcXZ = pts(:, [1,3]);
        [~, Z, ~] = procrustes(meanXZ, srcXZ, 'scaling', false, 'reflection', false);
        pts(:, [1,3]) = Z;
        DataMatrix3D(:, :, f) = pts;
    end
end

function plotAllFeet(DataMatrix3D)
    % Plots all feet in the given 3D data array
    n_feet = size(DataMatrix3D, 3);
    mFace = shape3D;
    v = viewer(mFace);
    for f = 1:n_feet
        shp = shape3D;
        shp.Vertices = DataMatrix3D(:,:,f);
        viewer(shp, v);
    end
end

function [coeff, score, latent, explained] = runPCA(DataMatrix3D)
    % Runs PCA on the data (each shape as a row vector)
    Y = arrayStructure2Vector(DataMatrix3D);
    [coeff, score, latent] = pca(Y);
    explained = latent / sum(latent) * 100;
end

function [synthConfig, full_pc_basis, ground_truth_pcs, footLengthVec_flat_scaled] = prepareSyntheticGeneration(DataMatrix3D, coeff)
    % Prepares PCA/PCs for synthetic shape generation and calibration
    Y = arrayToRowVector(DataMatrix3D);
    meanShapeVec = mean(Y, 1);
    nVertices = size(DataMatrix3D,1);
    meanShapeArray = rowVectorToArray(meanShapeVec);

    % Foot length mode
    x_coords = meanShapeArray(:,1);
    x_min = min(x_coords);
    x_max = max(x_coords);
    x_centered = 2 * (x_coords - x_min) / (x_max - x_min) - 1;
    footLengthVec = zeros(nVertices, 3);
    footLengthVec(:,1) = x_centered;

    % Normalize
    footLengthVec_flat = arrayToRowVector(footLengthVec);
    footLengthVec_flat = footLengthVec_flat / norm(footLengthVec_flat);

    % Calibration
    ground_truth_pcs_temp = [coeff(:,1), footLengthVec_flat'];
    shape_plus1  = meanShapeVec + [0, 1] * ground_truth_pcs_temp';
    shape_minus1 = meanShapeVec + [0, -1] * ground_truth_pcs_temp';
    foot_plus1  = rowVectorToArray(shape_plus1);
    foot_minus1 = rowVectorToArray(shape_minus1);
    len_plus1  = max(foot_plus1(:,1)) - min(foot_plus1(:,1));
    len_minus1 = max(foot_minus1(:,1)) - min(foot_minus1(:,1));
    delta_length_per_unit = (len_plus1 - len_minus1) / 2;
    fprintf('[PC2 calibration] Symmetric PC2 causes %.2f mm change per unit weight.\n', delta_length_per_unit);

    scaling_factor = 1 / delta_length_per_unit;
    footLengthVec_flat_scaled = footLengthVec_flat * scaling_factor;

    ground_truth_pcs = [coeff(:,1), footLengthVec_flat_scaled'];
    noise_pcs = coeff(:,3:5);
    full_pc_basis = [ground_truth_pcs, noise_pcs];

    synthConfig.meanShapeVec = meanShapeVec;
    synthConfig.nVertices = nVertices;
end

function [SyntheticData3D, W_full, syntheticNames] = generateSyntheticShapes(synthConfig, full_pc_basis, L_length)
    % Generates synthetic dataset for given PC2 variation (L_length)
    L_shape = 40;
    gridSize = ceil(sqrt(1000));
    range_shape  = linspace(-L_shape,  L_shape,  gridSize); 
    range_length = linspace(-L_length, L_length, gridSize); 
    [W1, W2] = meshgrid(range_shape, range_length);
    W = [W1(:), W2(:)]; 
    nSynthetic = size(W, 1);

    syntheticNames = arrayfun(@(i) sprintf('synthetic_%03d',i), 1:nSynthetic, 'UniformOutput',false);
    noise_sd = 15;
    rng(2, 'twister');
    noise_W = noise_sd * randn(nSynthetic, 3);  
    W_full = [ W, noise_W ];

    SyntheticData3D = zeros(synthConfig.nVertices, 3, nSynthetic);
    for i = 1:nSynthetic
        all_weights = W_full(i,:);
        shapeVec = synthConfig.meanShapeVec + all_weights * full_pc_basis';
        SyntheticData3D(:,:,i) = rowVectorToArray(shapeVec);
    end
end

function AlignedData = runStandardGPA(Data)
    nShapes = size(Data,3);
    maxIter = 10;
    tol = 1e-6;
    AlignedData = Data;
    meanPoints = mean(AlignedData,3);

    for iterGPA = 1:maxIter
        for f = 1:nShapes
            pts = AlignedData(:,:,f);
            [~, Z, ~] = procrustes(meanPoints, pts, 'scaling', true, 'reflection', false);
            AlignedData(:,:,f) = Z;
        end
        Mnew = mean(AlignedData,3);
        Mnew = Mnew - mean(Mnew,1);
        Mnew = Mnew / norm(Mnew,'fro');
        if norm(Mnew - meanPoints, 'fro') < tol, break; end
        meanPoints = Mnew;
    end
end

function AlignedData = runConstrainedGPA(Data)
    % Constrained GPA: Aligns shapes with only X/Z scaling and translation, and aligns soles
    % Data: [nVertices x 3 x nShapes] array

    SyntheticData3D_constrainedGPA = Data;
    nShapes = size(SyntheticData3D_constrainedGPA, 3);
    centroid_sizes = zeros(nShapes,1);

    % Center each shape (do not normalize size)
    for f = 1:nShapes
        pts = SyntheticData3D_constrainedGPA(:,:,f);
        c = mean(pts,1);
        centered = pts - c;
        centroid_sizes(f) = norm(centered,'fro');
        SyntheticData3D_constrainedGPA(:,:,f) = centered;
    end

    crudeMean = mean(SyntheticData3D_constrainedGPA, 3);

    % Choose the template closest to crude mean
    dists = zeros(nShapes,1);
    for f = 1:nShapes
        diff = SyntheticData3D_constrainedGPA(:,:,f) - crudeMean;
        dists(f) = norm(diff,'fro')^2;
    end
    [~, bestIdx] = min(dists);

    % Initialize meanPoints from that real-sized template
    meanPoints = SyntheticData3D_constrainedGPA(:,:,bestIdx);

    % Find sole indices once (on centered shapes)
    vertical_height_weights = 1;
    sole_indices = cell(1, nShapes);
    for f = 1:nShapes
        tmp = SyntheticData3D_constrainedGPA(:,:,f);
        sole_indices{f} = find(tmp(:,2) < vertical_height_weights);
    end

    % First iteration: zero out sole height
    for f = 1:nShapes
        verts = SyntheticData3D_constrainedGPA(:,:,f);
        sole_y = mean(verts(sole_indices{f},2));
        verts(:,2) = verts(:,2) - sole_y;
        SyntheticData3D_constrainedGPA(:,:,f) = verts;
    end
    meanPoints = mean(SyntheticData3D_constrainedGPA,3);

    % Main constrained-GPA loop
    iter = 10;
    scale = true; % horizontal scaling ON

    for i = 1:iter
        for f = 1:nShapes
            faceVerts = SyntheticData3D_constrainedGPA(:,:,f);

            % Only horizontal scale + rotation about y + horizontal translation
            T = computeConstrainedTransform(faceVerts, meanPoints, scale, sole_indices{f});
            faceVerts = applyTransformConstrained(faceVerts, T);

            % keep sole at y=0
            sole_y_after = mean(faceVerts(sole_indices{f},2));
            faceVerts(:,2) = faceVerts(:,2) - sole_y_after;

            SyntheticData3D_constrainedGPA(:,:,f) = faceVerts;
        end

        % update mean shape
        meanPoints = mean(SyntheticData3D_constrainedGPA,3);
    end

    % Output aligned data
    AlignedData = SyntheticData3D_constrainedGPA;
end

function [coeff, score, latent] = doPCA(Data)
    Y = arrayStructure2Vector(Data);
    [coeff, score, latent] = pca(Y);
end

function plotPCMorphs(coeff, latent, sharedMean, PCnum, k, GPA, variation_idx, variation_val, figureDir)
    % Visualize -3σ, mean, +3σ along a principal component
    Transform = arrayVector2Structure(coeff(:,PCnum)');
    sigma = sqrt(latent(PCnum));
    shapes = {sharedMean + Transform*sigma*-k, sharedMean, sharedMean + Transform*sigma*k};
    names = {'minus', 'mean', 'plus'};
    views = {[0,1,0],[0,0,1],[0,-1,0]};
    camups = {[0,0,0],[0,1,0],[0,0,0]};
    view_names = {'top', 'medial', 'bottom'};

    for iShape = 1:3
        m = shape3D; m.Vertices = shapes{iShape};
        v = viewer(m); set(gcf,'Color',[1 1 1]); set(gca,'Color',[1 1 1]);
        v.SceneLightLinked = true; v.SceneLightVisible = true;
        for j = 1:length(views)
            figHandle = v.Figure; figHandle.Color = [1 1 1];
            ax = findall(figHandle, 'Type', 'Axes'); set(ax, 'Color', [1 1 1]);
            view(views{j}); camup(camups{j});
            filename = sprintf('n%s_%d_PC%d_%s_%s_%s.png', GPA, num2str(variation_idx), PCnum, names{iShape}, GPA, view_names{j});
            exportgraphics(figHandle, fullfile(figureDir, filename), 'Resolution', 300, 'BackgroundColor', 'current');
        end
        close(v.Figure);
    end
end

function plotPC1Comparison(full_pc_basis, coeff_sGPA, coeff_cGPA, ...
    DataGT, DataStdGPA, DataConGPA, variation_idx, variation_val, figureDir)

    meanVerts = mean(DataGT,3);
    centroid = mean(meanVerts,1);
    Dgt  = reshape(full_pc_basis(:,1), [], 3);
    Dstd = reshape(coeff_sGPA(:,1), [], 3);
    Dcon = reshape(coeff_cGPA(:,1), [], 3);
    v_gt  = mean(Dgt,1);  v_gt  = v_gt  / norm(v_gt);
    v_std = mean(Dstd,1); v_std = v_std / norm(v_std);
    v_con = mean(Dcon,1); v_con = v_con / norm(v_con);
    if dot(v_con, v_gt) < 0, v_con = -v_con; end

    arrowLength = 50;
    isContact = meanVerts(:,2) < 1; % logical mask for contact region

    hFig = figure('Color', 'w'); hold on;
    scatter3(meanVerts(isContact,1), meanVerts(isContact,2), meanVerts(isContact,3), 5, [0.6 0.6 0.6], 'filled');
    scatter3(meanVerts(~isContact,1), meanVerts(~isContact,2), meanVerts(~isContact,3), 5, 'k', 'filled');
    quiver3(centroid(1), centroid(2), centroid(3), v_gt(1)*arrowLength, v_gt(2)*arrowLength, v_gt(3)*arrowLength, 0, 'r', 'LineWidth',2,'MaxHeadSize',0.8);
    quiver3(centroid(1), centroid(2), centroid(3), v_std(1)*arrowLength, v_std(2)*arrowLength, v_std(3)*arrowLength, 0, 'g', 'LineWidth',2,'MaxHeadSize',0.8);
    quiver3(centroid(1), centroid(2), centroid(3), v_con(1)*arrowLength, v_con(2)*arrowLength, v_con(3)*arrowLength, 0, 'b', 'LineWidth',2,'MaxHeadSize',0.8);
    view([0 -30 10]); camproj('orthographic');
    xlim([min(meanVerts(:,1)) max(meanVerts(:,1))]);
    ylim([min(meanVerts(:,2)) max(meanVerts(:,2))]);
    zlim([min(meanVerts(:,3)) max(meanVerts(:,3))]);
    set(gca, 'Position', [0 0 1 1]);
    axis equal off;
    set(gca, 'Visible', 'off', 'Box', 'off');
    figName = sprintf('meanShape_PC1_variation_%d.png', variation_val);
    saveas(hFig, fullfile(figureDir, figName));
    close(hFig);
end

function results = evaluatePCs(full_pc_basis, coeff_sGPA, coeff_cGPA, ...
    latent_sGPA, latent_cGPA, DataGT, DataStdGPA, DataConGPA, W_full, variation_idx, variation_val)
    % Evaluate recovery metrics for each PC

    MAX_NPCS = 5;  % Set this to your chosen maximum PCs
    results = struct();

    nPCsToUse = size(full_pc_basis,2);
    mean_ground_truth = mean(DataGT,3);
    mean_std = mean(DataStdGPA,3);
    mean_con = mean(DataConGPA,3);

    % ground-truth explained variance in shape space
    v_weights = var(W_full);
    v_pcnorms = vecnorm(full_pc_basis).^2;
    v_shape = v_weights .* v_pcnorms;
    v_shape_total = sum(v_shape);
    var_ground_truth_percent = v_shape / v_shape_total * 100;

    for k = 1:MAX_NPCS
        if k <= nPCsToUse
            ground_truth_k = full_pc_basis(:,k) / norm(full_pc_basis(:,k));
            est_std_k = coeff_sGPA(:,k) / norm(coeff_sGPA(:,k));
            est_con_k = coeff_cGPA(:,k) / norm(coeff_cGPA(:,k));
            if dot(ground_truth_k, est_std_k) < 0, est_std_k = -est_std_k; end
            if dot(ground_truth_k, est_con_k) < 0, est_con_k = -est_con_k; end
            cosine_similarity_std_k = dot(ground_truth_k, est_std_k);
            cosine_similarity_con_k = dot(ground_truth_k, est_con_k);

            if v_weights(k) > 0
                var_ground_truth_k = var_ground_truth_percent(k);
                var_std_k = latent_sGPA(k) / sum(latent_sGPA) * 100;
                var_con_k = latent_cGPA(k) / sum(latent_cGPA) * 100;
                delta_var_std = var_std_k - var_ground_truth_k;
                delta_var_con = var_con_k - var_ground_truth_k;
            else
                var_ground_truth_k = NaN; delta_var_std = NaN; delta_var_con = NaN;
            end
        else
            cosine_similarity_std_k = NaN;
            cosine_similarity_con_k = NaN;
            delta_var_std = NaN;
            delta_var_con = NaN;
        end
        results.(['cosine_similarity_std_k' num2str(k)]) = cosine_similarity_std_k;
        results.(['cosine_similarity_con_k' num2str(k)]) = cosine_similarity_con_k;
        results.(['delta_var_std' num2str(k)]) = delta_var_std;
        results.(['delta_var_con' num2str(k)]) = delta_var_con;
        if k==1
            if k <= nPCsToUse
                [d_std, ~, ~] = procrustes(mean_ground_truth, mean_std, 'Scaling', false, 'Reflection', false);
                [d_con, ~, ~] = procrustes(mean_ground_truth, mean_con, 'Scaling', false, 'Reflection', false);
                results.d_std = d_std;
                results.d_con = d_con;
            else
                results.d_std = NaN;
                results.d_con = NaN;
            end
        end
    end
    results.variation_idx = variation_idx;
    results.variation_val = variation_val;
end

function boot_results = bootstrapCosineSimilarities(SyntheticData3D, full_pc_basis, N_bootstrap)
    % Bootstraps cosine similarities between ground truth PC1 and estimated PC1s after GPA

    nShapes = size(SyntheticData3D,3);
    boot_cos_sim_std1 = zeros(N_bootstrap, 1);
    boot_cos_sim_con1 = zeros(N_bootstrap, 1);
    for b = 1:N_bootstrap
        idx = randsample(nShapes, nShapes, true);
        bootData = SyntheticData3D(:,:,idx);

        % Standard GPA and PCA
        bootData_sGPA = runStandardGPA(bootData);
        [bootCoeff_sGPA,~,~] = doPCA(bootData_sGPA);

        % Constrained GPA and PCA
        bootData_cGPA = runConstrainedGPA(bootData);
        [bootCoeff_cGPA,~,~] = doPCA(bootData_cGPA);

        gt_pc1 = full_pc_basis(:,1) / norm(full_pc_basis(:,1));
        est_std_pc1 = bootCoeff_sGPA(:,1) / norm(bootCoeff_sGPA(:,1));
        est_con_pc1 = bootCoeff_cGPA(:,1) / norm(bootCoeff_cGPA(:,1));

        % Sign flip for consistency
        if dot(gt_pc1, est_std_pc1) < 0, est_std_pc1 = -est_std_pc1; end
        if dot(gt_pc1, est_con_pc1) < 0, est_con_pc1 = -est_con_pc1; end

        boot_cos_sim_std1(b) = dot(gt_pc1, est_std_pc1);
        boot_cos_sim_con1(b) = dot(gt_pc1, est_con_pc1);
    end

    boot_results.std1.mean = mean(boot_cos_sim_std1);
    boot_results.std1.ci = prctile(boot_cos_sim_std1, [2.5 97.5]);
    boot_results.con1.mean = mean(boot_cos_sim_con1);
    boot_results.con1.ci = prctile(boot_cos_sim_con1, [2.5 97.5]);
end
