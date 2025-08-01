clear all;
clc;

% Pfadangaben DK
addpath(genpath('E:\2025_SSM_ArmaSuisse\Skripte\DK'))
tutorialPath = 'E:\2025_SSM_ArmaSuisse';
DataPath = 'E:\2025_SSM_ArmaSuisse\daten_matlab_punktkorrespondenz_obj_ohne_sohle';

chdir(tutorialPath);

temporary_number_objects = 2356; % Gesamtanzahl files: 2356

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
tolerance = 100; % Possible values (in mm): 0.5, 2.5, 5, 10, 20, 50

keepSlices = abs(foot_lengths - desiredLength) <= tolerance;
DataMatrix3D = DataMatrix3D(:,:,keepSlices);
temporary_number_objects = size(DataMatrix3D,3);

% Get names of files within the subset
files = dir(strcat(DataPath, '/*.obj')); % Struct containing files names
fileNames = {files.name};
keptFileNames = fileNames(keepSlices);

%% Loop over each desired sample size
sampleSizes = [25, 50, 100, 500, 1000, temporary_number_objects];
% rng(0); % or rng('shuffle') if you want differdraws each run

for s = 1:numel(sampleSizes)
    n = sampleSizes(s);
    
    % % Skip if the requested n is bigger than what you have
    % if n > temporary_number_objects
    %     warning('Requested %d samples, but only %d objects available—skipping.', ...
    %             n, temporary_number_objects);
    %     continue;
    % end
    
    % Randomly pick n distinct indices from 1:Nkept
    idx = randperm(temporary_number_objects, n);
    
    % Extract your subset
    DataMatrixSample = DataMatrix3D(:, :, idx);        % [nVerts×3×n]
    FileNamesSample    = keptFileNames(idx);           % 1×n cell
    
    % — now you can run whatever analysis you need on this sample —
    % e.g. compute mean shape:
    meanShape = mean(DataMatrixSample, 3);             % [nVerts×3]
    
    % or save them, plot them, etc.:
    fprintf('Sample size = %d: mean foot length = %.1f mm\n', ...
            n, mean(meanShape(:,1)) - min(meanShape(:,1)));
    
    SyntheticData3D = DataMatrixSample;

    % Compute the shared mean shape ONCE after generating the synthetic data
    meanShapeVec = mean(arrayToRowVector(SyntheticData3D), 1); 
    sharedMean = rowVectorToArray(meanShapeVec);  % N × 3

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

    % %% Füße plotten
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

    %% Create and store shape deformation plots
    PCnum = 1;                % first principal component
    k     = 3;                % ±3σ
    variation_idx = 2;        % whichever variation you’re visualizing

    exportPath = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\visualize_PCs_sample_size';

    visualizePCMorph(coeff_sGPA, latent_sGPA, PCnum, sharedMean, k, ...
                 'standardGPA', temporary_number_objects, variation_idx, ...
                 exportPath);

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
    scale = true; % keep horizontal (anisotropic) scaling ON
    
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

    %% Principal Components Analysis

    % Represent each foot as a row of a 2D matrix
    Y = arrayStructure2Vector(SyntheticData3D_constrainedGPA);
    
    % Y is already mean-centered data (Skript DK)!
    [coeff_cGPA, score_cGPA, latent_cGPA] = pca(Y); % cGPA := constrained GPA
    
    % Percentage of total variance explained by each PCverify if PC1 really captures foot length
    explained_cGPA = latent_cGPA / sum(latent_cGPA) * 100;

    %% Create and store shape deformation plots
    PCnum = 1;                % first principal component
    k     = 3;                % ±3σ
    variation_idx = 2;        % whichever variation you’re visualizing

    exportPath = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\visualize_PCs_sample_size';

    visualizePCMorph(coeff_cGPA, latent_sGPA, PCnum, sharedMean, k, ...
                 'constrainedGPA', temporary_number_objects, variation_idx, ...
                 exportPath);

    %% PCs der drei Methoden vergleichsweise plotten

    % --- Plot Mean Shape + PC1 Arrows for GT, Std‐GPA & Constrained‐GPA --
    
    % Compute the mean‐shape cloud and its centroid
    meanShapeVec = mean(arrayToRowVector(SyntheticData3D), 1);
    meanVerts = rowVectorToArray(meanShapeVec); % N×3
    centroid = mean(meanVerts,1); % 1×3
    
    % Extract the three PC1 vectors and reshape to N×3
    Dstd = reshape(coeff_sGPA(:,1), [], 3);
    Dcon = reshape(coeff_cGPA(:,1), [], 3);
    
    % Collapse each field to one 1×3 direction and normalize
    v_std = mean(Dstd, 1); v_std = v_std / norm(v_std);
    v_con = mean(Dcon, 1); v_con = v_con / norm(v_con);
    
    arrowLength = 50;
    
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
    output_dir = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\PC_vectors_simulation_sample_size';
    figName = sprintf('meanShape_PC1_variation_%d', s);
   
    % Save as PNG
    saveas(hFig, fullfile(output_dir, [figName '.png']));
end