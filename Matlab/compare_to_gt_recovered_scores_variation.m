clear all;
clc;

% Pfadangaben DK
addpath(genpath('E:\2025_SSM_ArmaSuisse\Skripte\DK'))
tutorialPath = 'E:\2025_SSM_ArmaSuisse';

% Condition: Feet with or without sole (sole removed in external R script 'ply_obj_delete_sole.R')
sole_condition = 'sole'; % 'sole' - with sole points, 'no_sole' - w/o sole points
if strcmp(sole_condition, 'sole')
    DataPath = 'E:\2025_SSM_ArmaSuisse\daten_matlab_punktkorrespondenz_obj';
else
    DataPath = 'E:\2025_SSM_ArmaSuisse\daten_matlab_punktkorrespondenz_obj_ohne_sohle';
end

chdir(tutorialPath);

temporary_number_objects = 2356; % Gesamtanzahl files: 2356
GPA = "noGPA"; % Variable for file (export) naming

%% Load all faces
% Hier werden jetzt die zuvor im Skript correspondence_meshmonk_batch.m
% erstellten (korrespondieren) Pointclouds geladen
objs = dir(strcat(DataPath,filesep,'*.obj'));

% remove any files prefaced with a '._' they are not real
isreal = ~contains({objs.name},'._');
objs = objs(isreal);

for i = 1:temporary_number_objects
   filePath = fullfile(DataPath, objs(i).name);
   
   % import as a shape 3D
   obj = shape3D;
   importWavefront(obj, objs(i).name, objs(i).folder, []); % Lädt das eigentliche obj file
   if i==1
        % initialise 3D matrix to contain vertices of all faces
        DataMatrix3D = zeros(obj.nVertices,3,temporary_number_objects); % Variante mit obj's
   end
   DataMatrix3D(:,:,i) = obj.Vertices;
end

%% Filter by foot length
foot_lengths = zeros(temporary_number_objects,1);

for i = 1:temporary_number_objects
    foot = DataMatrix3D(:, :, i); % Get current foot [2000 x 3]
    
    % Compute foot length along x-axis (x is anterior-posterior)
    x_coords = foot(:, 1);
    foot_lengths(i) = max(x_coords) - min(x_coords);
end

% Define a threshold range around a desired foot length in mm
desiredLength = 270; % 270 mm corresponds to (roughly) EU size 42
tolerance = 0.5; % Possible values (in mm): 0.5, 2.5, 5, 10, 20, 50

keepSlices = abs(foot_lengths - desiredLength) <= tolerance;
DataMatrix3D = DataMatrix3D(:,:,keepSlices);
temporary_number_objects = size(DataMatrix3D,3);

% Get names of files within the subset
files = dir(strcat(DataPath, '/*.obj')); % Struct containing files names
fileNames = {files.name};
keptFileNames = fileNames(keepSlices);

%% Translation‐only alignment (remove horizontal offsets only)

nShapes = size(DataMatrix3D,3);
for f = 1:nShapes
    pts = DataMatrix3D(:,:,f);
    
    % compute mean in x and z only
    mean_xz = mean(pts(:,[1,3]), 1);
    
    % subtract it so each shape is centered in the horizontal plane
    pts(:,[1,3]) = pts(:,[1,3]) - mean_xz;
    
    DataMatrix3D(:,:,f) = pts;
end

% new mean shape (will still reflect vertical differences)
meanPoints = mean(DataMatrix3D,3);

%% Füße plotten
n_feet = size(DataMatrix3D, 3); % number of feet to visualize

mFace = shape3D;
v=viewer(mFace);
for f = 1:n_feet
    shp = shape3D; % create shape 3D
    shp.Vertices = DataMatrix3D(:,:,f);

    viewer(shp,v);
end

%% Principal Components Analysis

% Represent each foot as a row of a 2D matrix
Y = arrayStructure2Vector(DataMatrix3D);

% Y is already mean-centered data (Skript DK)!
[coeff, score, latent] = pca(Y);

% Percentage of total variance explained by each PC
explained = latent / sum(latent) * 100;

%% Generate synthetic shapes from PCA variances
Y = arrayToRowVector(DataMatrix3D);

meanShapeVec = mean(Y, 1);
nPCsToUse = 5;

nVertices = size(DataMatrix3D,1);

variation_levels = [2 5 10 20 30 40 50 100 200 400];
results = struct();

for variation_idx = 1:length(variation_levels)

    % Build synthetic shapes by sampling PC1 & PC2 weights ----------------
    L = variation_levels(variation_idx);

    useDeterministicWeights = true;

    if useDeterministicWeights
        % Generate 25 deterministic weight combinations (5×5 grid)
        range = linspace(-L, L, 5);
        [W1, W2] = meshgrid(range, range);
        W = [W1(:), W2(:)]; % 25 × 2 matrix
        nSynthetic = size(W, 1); % Number of synthetic feet
    else
        % Random weights in [-L, L] for PC1 and PC2
        W = (2*rand(nSynthetic,2) - 1) * L;
    end 
    
    % Zero out all other PCs
    W_full = zeros(nSynthetic, nPCsToUse);
    W_full(:,1:2) = W;  
    
    SyntheticData3D = zeros(nVertices, 3, nSynthetic);
    ground_truth_pcs = coeff(:,1:nPCsToUse);
    
    for i = 1:nSynthetic
        weights = W_full(i,:); % 1×nPCsToUse
        shapeVec = meanShapeVec + weights * ground_truth_pcs';
        SyntheticData3D(:,:,i) = rowVectorToArray(shapeVec);
    end

    %% Füße plotten
    n_feet = size(SyntheticData3D, 3); % number of feet to visualize

    mFace = shape3D;
    v=viewer(mFace);

    for f = 1:n_feet
        shp = shape3D;
        shp.Vertices = SyntheticData3D(:,:,f);

        viewer(shp,v);
    end

    set(gcf, 'Color', 'w'); % white background
    set(gca, 'Color', 'w'); % white axes background
    axis tight off; % remove axes and fit bounds
    camproj('orthographic'); % optional: parallel view
    set(gca, 'Position', [0 0 1 1]);       % remove all margins

    % Save figure
    output_dir = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\simulated_point_clouds';
    variation_val = variation_levels(variation_idx);
    figName = sprintf('pointcloud_variation_%d', variation_val);

    % Save as clean PNG using exportgraphics
    exportgraphics(gca, fullfile(output_dir, [figName '.png']), ...
        'BackgroundColor', 'white', ...
        'ContentType', 'image', ...
        'Resolution', 300);  % Optional: high DPI
    
    % Close the figure if desired
    close(gcf);
    
    %% Standard GPA
    
    SyntheticData3D_standardGPA = SyntheticData3D;
    
    GPA = "standardGPA"; % for export naming
    
    Template = shape3D;
    Template.Vertices = SyntheticData3D(:,:,1);
    
    % Center each shape, record its centroid‐size, and normalize to unit‐size
    nShapes = size(SyntheticData3D_standardGPA,3);
    centroid_sizes = zeros(nShapes,1);
    
    for f = 1:nShapes
        pts = SyntheticData3D_standardGPA(:,:,f);
        c = mean(pts,1);
        centered = pts - c;    
        centroid_sizes(f) = norm(centered,'fro');
        % normalize to unit‐size for template selection & initial GPA
        SyntheticData3D_standardGPA(:,:,f) = centered / centroid_sizes(f);
    end
    
    % Compute the crude mean (unit‐size)
    crudeMean = mean(SyntheticData3D_standardGPA,3);
    
    % Pick the shape closest (in Frobenius distance) to that mean
    dists = zeros(nShapes,1);
    for f = 1:nShapes
        diff = SyntheticData3D_standardGPA(:,:,f) - crudeMean;
        dists(f) = norm(diff,'fro')^2;
    end
    [~, bestIdx] = min(dists);
    
    % log which file is your starting template
    templateFileName = keptFileNames{bestIdx};
    fprintf('Using %s as the starting template (unit‐normalized).\n', templateFileName);
    
    % Initialize meanPoints from that selected unit‐size shape
    meanPoints = SyntheticData3D_standardGPA(:,:,bestIdx);
    
    % Iterative GPA, aligning & re‐scaling each to C_ref every iteration
    iter  = 10;   
    for i = 1:iter
        for f = 1:nShapes
            % align to current meanPoints (which is at size = C_ref)
            faceVerts = SyntheticData3D_standardGPA(:,:,f);           
            T = computeTransform(faceVerts, meanPoints, true);
            faceVerts = applyTransform(faceVerts, T);
            
            % re‐center (no scaling)
            centered = faceVerts - mean(faceVerts,1);
            SyntheticData3D_standardGPA(:,:,f) = centered;
        end
        
        % update the mean shape and re‐enforce size = C_ref
        meanPoints = mean(SyntheticData3D_standardGPA,3);
        meanPoints = meanPoints / norm(meanPoints,'fro');
    end
    
    %% Füße plotten
    % n_feet = size(SyntheticData3D_standardGPA, 3); % number of feet to visualize
    % 
    % mFace = shape3D;
    % v=viewer(mFace);
    % for f = 1:n_feet
    %     shp = shape3D; % create shape 3D
    %     shp.Vertices = SyntheticData3D_standardGPA(:,:,f);
    % 
    %     viewer(shp,v);
    % end
    
    %% Principal Components Analysis
    
    % Represent each foot as a row of a 2D matrix
    Y = arrayStructure2Vector(SyntheticData3D_standardGPA);
    
    % Y is already mean-centered data (Skript DK)!
    [coeff_sGPA, score_sGPA, latent_sGPA] = pca(Y); % sGPA := standard GPA
    
    % Percentage of total variance explained by each PC
    explained_sGPA = latent_sGPA / sum(latent_sGPA) * 100;
    
    %% Constrained Generalized Procrustes Analysis (DK)
    
    % Start from your synthetic shapes (with true length variation)
    SyntheticData3D_constrainedGPA = SyntheticData3D;
    
    GPA = "constrainedGPA";  % for export naming
    
    nShapes = size(SyntheticData3D_constrainedGPA, 3);
    centroid_sizes = zeros(nShapes,1);
    
    % Center each shape and record its true centroid‐size,
    % but do NOT normalize to unit‐size here.
    for f = 1:nShapes
        pts = SyntheticData3D_constrainedGPA(:,:,f);
        c = mean(pts,1);
        centered = pts - c;
        centroid_sizes(f) = norm(centered,'fro');       % real size in mm
        SyntheticData3D_constrainedGPA(:,:,f) = centered;
    end
    
    crudeMean = mean(SyntheticData3D_constrainedGPA, 3);
    
    % 3) Choose the template closest to that mean
    dists = zeros(nShapes,1);
    for f = 1:nShapes
        diff = SyntheticData3D_constrainedGPA(:,:,f) - crudeMean;
        dists(f) = norm(diff,'fro')^2;
    end
    [~, bestIdx] = min(dists);

    templateFileName = keptFileNames{bestIdx};
    fprintf('Using %s as starting template.\n', templateFileName);
    
    % Initialize meanPoints from that real‐sized template
    meanPoints = SyntheticData3D_constrainedGPA(:,:,bestIdx);
    
    % Find sole indices once (on centered shapes)
    vertical_height_weights = 1;  
    sole_indices = cell(1, nShapes);
    for f = 1:nShapes
        tmp = SyntheticData3D_constrainedGPA(:,:,f);
        sole_indices{f} = find(tmp(:,2) < vertical_height_weights);
    end
    
    % First iteration: zero‐out sole height
    for f = 1:nShapes
        verts = SyntheticData3D_constrainedGPA(:,:,f);
        sole_y = mean(verts(sole_indices{f},2));
        verts(:,2) = verts(:,2) - sole_y;
        SyntheticData3D_constrainedGPA(:,:,f) = verts;
    end
    meanPoints = mean(SyntheticData3D_constrainedGPA,3);
    
    % Main constrained‐GPA loop
    iter  = 10;
    scale = true;  % keep horizontal (anisotropic) scaling ON
    
    for i = 1:iter
        for f = 1:nShapes
            faceVerts = SyntheticData3D_constrainedGPA(:,:,f);
            
            % apply only horizontal scale + rotation about y + horizontal translation
            T = computeConstrainedTransform(faceVerts, meanPoints, scale, sole_indices{f});
            faceVerts = applyTransformConstrained(faceVerts, T);
            
            % keep sole at y=0
            sole_y_after = mean(faceVerts(sole_indices{f},2));
            faceVerts(:,2) = faceVerts(:,2) - sole_y_after;
            
            % store the transformed shape (no further scaling)
            SyntheticData3D_constrainedGPA(:,:,f) = faceVerts;
        end
        
        % update the mean shape (retaining true scale variation)
        meanPoints = mean(SyntheticData3D_constrainedGPA,3);
    end
    
    %% Füße plotten
    % n_feet = size(SyntheticData3D_constrainedGPA, 3); % number of feet to visualize
    % 
    % mFace = shape3D;
    % v=viewer(mFace);
    % for f = 1:n_feet
    %     shp = shape3D; % create shape 3D
    %     shp.Vertices = SyntheticData3D_constrainedGPA(:,:,f);
    % 
    %     viewer(shp,v);
    % end
    
    %% Principal Components Analysis
    
    % Represent each foot as a row of a 2D matrix
    Y = arrayStructure2Vector(SyntheticData3D_constrainedGPA);
    
    % Y is already mean-centered data (Skript DK)!
    [coeff_cGPA, score_cGPA, latent_cGPA] = pca(Y); % cGPA := constrained GPA
    
    % Percentage of total variance explained by each PCverify if PC1 really captures foot length
    explained_cGPA = latent_cGPA / sum(latent_cGPA) * 100;
    
    %% PCs der drei Methoden vergleichsweise plotten
    
    % --- Plot Mean Shape + PC1 Arrows for GT, Std‐GPA & Constrained‐GPA --

    % Compute the mean‐shape cloud and its centroid
    meanShapeVec = mean(arrayToRowVector(SyntheticData3D), 1);
    meanVerts = rowVectorToArray(meanShapeVec); % N×3
    centroid = mean(meanVerts,1); % 1×3
    
    % Extract the three PC1 vectors and reshape to N×3
    Dgt  = reshape(ground_truth_pcs(:,1), [], 3);
    Dstd = reshape(coeff_sGPA(:,1), [], 3);
    Dcon = reshape(coeff_cGPA(:,1), [], 3);
    
    % Collapse each field to one 1×3 direction and normalize
    v_gt  = mean(Dgt,  1); v_gt  = v_gt  / norm(v_gt);
    v_std = mean(Dstd, 1); v_std = v_std / norm(v_std);
    v_con = mean(Dcon, 1); v_con = v_con / norm(v_con);
    
    arrowLength = 50; % mm
    
    % Define contact region (e.g., y < 1 mm)
    isContact = meanVerts(:,2) < 1; % Logical index of contact points
    
    % Plot without title or legend
    hFig = figure('Color', 'w');
    hold on;
    
    % Plot contact region points (e.g., in gray)
    scatter3(meanVerts(isContact,1), meanVerts(isContact,2), meanVerts(isContact,3), ...
             5, [0.6 0.6 0.6], 'filled');
    
    % Plot non-contact points (in black)
    scatter3(meanVerts(~isContact,1), meanVerts(~isContact,2), meanVerts(~isContact,3), ...
             5, 'k', 'filled');
    
    % PC1 arrows
    quiver3(centroid(1), centroid(2), centroid(3), ...
            v_gt(1)*arrowLength, v_gt(2)*arrowLength, v_gt(3)*arrowLength, ...
            0, 'r', 'LineWidth',2, 'MaxHeadSize',0.8);
    quiver3(centroid(1), centroid(2), centroid(3), ...
            v_std(1)*arrowLength, v_std(2)*arrowLength, v_std(3)*arrowLength, ...
            0, 'g', 'LineWidth',2, 'MaxHeadSize',0.8);
    quiver3(centroid(1), centroid(2), centroid(3), ...
            v_con(1)*arrowLength, v_con(2)*arrowLength, v_con(3)*arrowLength, ...
            0, 'b', 'LineWidth',2, 'MaxHeadSize',0.8);
    
    % View and minimalistic display
    view([0 -30 10]);          
    camproj('orthographic');   
    rotate3d on;
    
    % Tighten axes around point cloud — ADD THIS
    xlim([min(meanVerts(:,1)) max(meanVerts(:,1))]);
    ylim([min(meanVerts(:,2)) max(meanVerts(:,2))]);
    zlim([min(meanVerts(:,3)) max(meanVerts(:,3))]);
    
    % Make axes fill the figure area — ADD THIS
    set(gca, 'Position', [0 0 1 1]);
    
    % Hide everything
    axis equal off;
    set(gca, 'Visible', 'off', 'Box', 'off');
    set(hFig, 'MenuBar', 'none', 'ToolBar', 'none');

    % Save figure
    output_dir = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\PC_vectors_simulation';
    variation_val = variation_levels(variation_idx);
    figName = sprintf('meanShape_PC1_variation_%d', variation_val);

    % Save as MATLAB figure
    savefig(hFig, fullfile(output_dir, [figName '.fig']));

    % Save as PNG
    saveas(hFig, fullfile(output_dir, [figName '.png']));

    close(hFig);

    %% --- EVALUATION METRICS -----------------------------------------

    % Precompute means for Procrustes once
    mean_ground_truth  = rowVectorToArray(meanShapeVec);     % original mean
    mean_std = mean(SyntheticData3D_standardGPA,3);          % std‐GPA mean
    mean_con = mean(SyntheticData3D_constrainedGPA,3);       % con‐GPA mean

    % Compute ground-truth explained variance in SHAPE space
    v_weights  = var(W_full);                               % variance of PC weights
    v_pcnorms  = vecnorm(ground_truth_pcs).^2;              % squared norm of each PC in shape space
    v_shape    = v_weights .* v_pcnorms;                    % variance contribution per PC in shape space
    v_shape_total = sum(v_shape);                           % total shape variance from all varied PCs
    var_ground_truth_percent = v_shape / v_shape_total * 100;         % explained by each ground-truth PC

    fprintf('\n=== Recovery Metrics for variation L=%g ===\n', variation_levels(variation_idx));
    for k = 1:nPCsToUse
        % ground‐truth PC
        ground_truth_k      = ground_truth_pcs(:,k);
        ground_truth_k      = ground_truth_k / norm(ground_truth_k);

        % recovered PC
        est_std_k = coeff_sGPA(:,k);
        est_std_k = est_std_k / norm(est_std_k);
        est_con_k = coeff_cGPA(:,k);
        est_con_k = est_con_k / norm(est_con_k);

        % cosine similarity (absolute)
        cosine_similarity_std_k = abs(dot(ground_truth_k, est_std_k));
        cosine_similarity_con_k = abs(dot(ground_truth_k, est_con_k));

        % variance explained in percent
        % Only for PCs actually varied, compare to synthetic shape-space variance)
        if v_weights(k) > 0
            var_ground_truth_k = var_ground_truth_percent(k);  % ground-truth variance
            var_std_k = latent_sGPA(k) / sum(latent_sGPA) * 100;
            var_con_k = latent_cGPA(k) / sum(latent_cGPA) * 100;
            delta_var_std = var_std_k - var_ground_truth_k;
            delta_var_con = var_con_k - var_ground_truth_k;
        else
            var_ground_truth_k = NaN;
            delta_var_std = NaN;
            delta_var_con = NaN;
        end

        % Procrustes distance only for k==1
        if k==1
            [d_std, ~, ~] = procrustes(mean_ground_truth, mean_std, 'Scaling', false, 'Reflection', false);
            [d_con, ~, ~] = procrustes(mean_ground_truth, mean_con, 'Scaling', false, 'Reflection', false);
        end

        fprintf('\nPC%-2d recovery:\n', k);
        fprintf('  cosine sim:    std = %.3f   con = %.3f\n', cosine_similarity_std_k, cosine_similarity_con_k);
        fprintf('  Δvariance (%%): std = %+5.2f   con = %+5.2f\n', delta_var_std, delta_var_con);
        if k==1
            fprintf('  mean‐shape dist: std = %.4f   con = %.4f\n', d_std, d_con);
        end

        % Store
        results(variation_idx).(['cosine_similarity_std_k' num2str(k)]) = cosine_similarity_std_k;
        results(variation_idx).(['cosine_similarity_con_k' num2str(k)]) = cosine_similarity_con_k;
        results(variation_idx).(['delta_var_std' num2str(k)]) = delta_var_std;
        results(variation_idx).(['delta_var_con' num2str(k)]) = delta_var_con;
        if k==1
            results(variation_idx).d_std = d_std;
            results(variation_idx).d_con = d_con;
        end

        % Save as csv
        output_dir = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\PC_vectors_simulation';
        results_table = struct2table(results);
    
        % Use the variation level directly in the filename
        variation_val = variation_levels(variation_idx);
        csv_filename = fullfile(output_dir, sprintf('results_variation_%d.csv', variation_val));
        writetable(results_table, csv_filename);
    end

end % Ende der variation_levels Schleife










