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
tolerance = 0.1; % Possible values (in mm): 0.5, 2.5, 5, 10, 20, 50

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

mFace = shape3D; % shape3d = meshmonk object
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
nVertices = size(DataMatrix3D,1);

% Define clean foot-length mode (subsequently used as PC2) ----------------

meanShapeArray = rowVectorToArray(meanShapeVec); % Convert mean shape to array

% Define symmetric foot-length mode (to be used as PC2)
x_coords = meanShapeArray(:,1);
x_min = min(x_coords);
x_max = max(x_coords);

% Normalize x-coordinates to range [-1, 1] for uniform proportional stretching
x_centered = 2 * (x_coords - x_min) / (x_max - x_min) - 1;

footLengthVec = zeros(nVertices, 3);
footLengthVec(:,1) = x_centered;

% Flatten and normalize
footLengthVec_flat = arrayToRowVector(footLengthVec);
footLengthVec_flat = footLengthVec_flat / norm(footLengthVec_flat);

% Calibrate PC2 so that 1 unit = 1 mm foot length change
% Temporary ground truth matrix (coeff(:,1) is still PC1)
ground_truth_pcs_temp = [coeff(:,1), footLengthVec_flat'];

% Generate shapes with ±1 unit of PC2
shape_plus1  = meanShapeVec + [0, 1] * ground_truth_pcs_temp';
shape_minus1 = meanShapeVec + [0, -1] * ground_truth_pcs_temp';

% Convert to 3D shape arrays
foot_plus1  = rowVectorToArray(shape_plus1);
foot_minus1 = rowVectorToArray(shape_minus1);

% Compute foot length as x-range
len_plus1  = max(foot_plus1(:,1)) - min(foot_plus1(:,1));
len_minus1 = max(foot_minus1(:,1)) - min(foot_minus1(:,1));

% Foot length difference caused by unit PC2
delta_length_per_unit = (len_plus1 - len_minus1) / 2;
fprintf('[PC2 calibration] Symmetric PC2 causes %.2f mm change per unit weight.\n', delta_length_per_unit);

% Rescale to make 1 unit weight = 1 mm change in foot length
scaling_factor = 1 / delta_length_per_unit;
footLengthVec_flat_scaled = footLengthVec_flat * scaling_factor;

% Final ground truth PCs (PC1 from real PCA, PC2 = calibrated foot length)
ground_truth_pcs = [coeff(:,1), footLengthVec_flat_scaled'];

% Use the next three PCs from the original pca(Y) as low-variance noise modes
noise_pcs = coeff(:,3:5); % size 3n×3

% Build a 3n×5 basis: [PC1, PC2, PC3, PC4, PC5]
full_pc_basis = [ground_truth_pcs, noise_pcs];

foot_length_levels = [0, 2, 5, 10, 20, 30, 40, 50]; % in mm

nPCsToUse = size(ground_truth_pcs, 2); % 5;

results = struct();

for variation_idx = 1:length(foot_length_levels)

    % Build synthetic shapes by sampling PC1 & PC2 weights ----------------
    
    % Define the variation amplitude for PC1 (arch height) and PC2 (foot length)
    L_shape = 40; % fixed shape variation along PC1
    L_length = foot_length_levels(variation_idx); % variable: foot length variation along PC2

    gridSize = ceil(sqrt(1000)); % 5 (25), 10 (100), ceil(sqrt(1000)) (1024), 
    range_shape  = linspace(-L_shape,  L_shape,  gridSize); % arch variation
    range_length = linspace(-L_length, L_length, gridSize); % foot length variation
    [W1, W2] = meshgrid(range_shape, range_length);     % W2 now controls foot length
    W = [W1(:), W2(:)]; % 25 × 2 matrix
    
    nSynthetic = size(W, 1);

    syntheticNames = arrayfun(@(i) sprintf('synthetic_%03d',i), 1:nSynthetic, ...
                          'UniformOutput',false);
    
    noise_sd = 15; % mm
    
    rng(2, 'twister');
    noise_W = noise_sd * randn(nSynthetic, 3);  % 25×3

    W_full = [ W, noise_W ];

    % Variante mit noise
    SyntheticData3D = zeros(nVertices, 3, nSynthetic);

    for i = 1:nSynthetic
        all_weights = W_full(i,:);

        shapeVec = meanShapeVec + all_weights * full_pc_basis';
        
        SyntheticData3D(:,:,i) = rowVectorToArray(shapeVec);
    end

    %% What does a foot length variation level of 2 correspond to in terms of EU shoe sizes?
    % % Generate foot with variation level 0 (mean shape)
    % foot0 = meanShapeVec;
    % 
    % % Generate foot with variation level 2 (lengthened)
    % foot2 = meanShapeVec + foot_length_levels(variation_idx) * footLengthVec_flat_scaled;
    % 
    % % Reshape and extract x-length
    % foot0_array = rowVectorToArray(foot0);
    % foot2_array = rowVectorToArray(foot2);
    % 
    % len0 = max(foot0_array(:,1)) - min(foot0_array(:,1));
    % len2 = max(foot2_array(:,1)) - min(foot2_array(:,1));
    % 
    % delta_mm = len2 - len0;

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
    set(gca, 'Position', [0 0 1 1]); % remove all margins

    % % Save figure
    % output_dir = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\simulated_point_clouds';
    % variation_val = foot_length_levels(variation_idx);
    % figName = sprintf('pointcloud_variation_%d', variation_val);
    % 
    % % Save as clean PNG using exportgraphics
    % exportgraphics(gca, fullfile(output_dir, [figName '.png']), ...
    %     'BackgroundColor', 'white', ...
    %     'ContentType', 'image', ...
    %     'Resolution', 300);
    % 
    % close(gcf);

    % %% Get the range of x-coordinates as a proxy for foot length from a 3D point cloud
    % 
    % % x := length axis of the foot
    % foot_lengths = zeros(n_feet,1);
    % 
    % for i = 1:n_feet
    %     x_coords = DataMatrix3D(:,1,i); % x-axis slice for foot i
    %     foot_lengths(i) = max(x_coords) - min(x_coords);
    % end
    % 
    % range_length = [min(foot_lengths), max(foot_lengths)];
    % disp(['Foot length range: ', num2str(range_length(1)), '–', num2str(range_length(2)), ' mm']);


    %% Visualise principal component by transforming average face

    % Compute the shared mean shape ONCE after generating the synthetic data
    meanShapeVec = mean(arrayToRowVector(SyntheticData3D), 1); 
    sharedMean = rowVectorToArray(meanShapeVec);  % N × 3

    GPA = 'groundtruth';
    
    PCnum = 1; % which PC to plot
    
    Transform = arrayVector2Structure(coeff(:,PCnum)'); % Reshapes that 1×(3N) row into an N×3 structure (or array of structs) so you now have a 3D displacement vector at each vertex.
    sigma = sqrt(latent(PCnum)); % latent holds the variances (eigenvalues) of each principal component
    k = 3; % Scaling factor
    
    AverageShape = shape3D;
    AverageShape.Vertices = sharedMean;
    
    % Apply PC to visualize -3σ and +3σ deformations
    PCminusMorph = clone(AverageShape); % clone makes a deep copy of the mean shape
    PCminusMorph.Vertices = sharedMean + Transform*sigma*-k;
    
    PCmeanMorph = clone(AverageShape); % Mean shape, no deformation
    PCmeanMorph.Vertices = sharedMean;
    
    PCplusMorph = clone(AverageShape);
    PCplusMorph.Vertices = sharedMean + Transform*sigma*k;
    
    % Set the mesh color to black
    PCminusMorph.SingleColor = [0 0 0];
    PCmeanMorph.SingleColor = [0 0 0];
    PCplusMorph.SingleColor = [0 0 0];
    
    % Create three plots (-3*sigma, mean, 3*sigma). Plot each in a separate
    % viewer.
    
    % Viewer 1 (mean - 3*sigma) -----------------------------------------------
    v1 = viewer(PCminusMorph);
    set(gcf, 'Color', [1 1 1]); % Set white background
    set(gca, 'Color', [1 1 1]); % Set white background
    % title('PCminusMorph');
    v1.SceneLightLinked = true;
    v1.SceneLightVisible = true;
    
    % Export as .fig
    cd 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\visualize_PCs';
    
    se0 = "n";
    se1 = num2str(temporary_number_objects);
    se2 = "_";
    se3 = num2str(variation_idx);
    se4 = "_PC";
    se5 = num2str(PCnum);
    se6 = "_minus_morph_";
    se7 = GPA;
    se8 = ".fig";
    sges= strcat(se0,se1,se2,se3,se4,se5,se6,se7,se8);
    
    % savefig(gcf, sges);
    
    % Export as .png
    
    % Variante 1 (mit top, bottom, und medial view)
    % Define the specific desired views
    angles = {
        [0, 1, 0];  % Top view (bird perspective, transversal plane)
        [0, 0, 1];  % Medial side view (sagittal plane, medial view)
        [0, -1, 0];  % Bottom view (observer below the foot)
    };
    
    % Corresponding camera-up vectors
    camups = {
        [0, 0, 0];  % top view: heel bottom, toes top
        [0, 1, 0];  % medial side view: vertical axis up
        [0, 0, 0];  % bottom view: heel top, toes bottom
    };
    
    view_names = {'top', 'medial', 'bottom'};
    
    for i = 1:length(angles)
        figHandle = v1.Figure;
        figHandle.Color = [1 1 1]; % White background
        ax = findall(figHandle, 'Type', 'Axes');
        set(ax, 'Color', [1 1 1]);
    
        % Set view angles explicitly
        view(angles{i});
        camup(camups{i});
    
        % Construct a clear filename for export
        filename = sprintf('%s%s_%s_PC%d_minus_%s_%s.png', se0, se1, se3, PCnum, GPA, view_names{i});
    
        % Export as PNG
        exportgraphics(figHandle, filename, 'Resolution', 300, 'BackgroundColor', 'current');
    end
    
    
    % Viewer 2 (mean) ---------------------------------------------------------
    v2 = viewer(PCmeanMorph);
    
    % Set white background
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'Color', [1 1 1]);
    
    % title('PCmeanMorph');
    v2.SceneLightLinked = true;
    v2.SceneLightVisible = true;
    
    se8 = ".fig";
    sges= strcat(se0,se1,se2,se3,se4,se5,se6,se7,se8);
    
    % Export as .fig
    savefig(gcf, sges);
    
    % Export as .png
    for i = 1:length(angles)
        figHandle = v2.Figure;
        figHandle.Color = [1 1 1]; % White background
        ax = findall(figHandle, 'Type', 'Axes');
        set(ax, 'Color', [1 1 1]);
    
        % Set view angles explicitly
        view(angles{i});
        camup(camups{i});
    
        % Construct a clear filename for export
        filename = sprintf('%s%s_%s_PC%d_mean_%s_%s.png', se0, se1, se3, PCnum, GPA, view_names{i});
    
        % Export as PNG
        exportgraphics(figHandle, filename, 'Resolution', 300, 'BackgroundColor', 'current');
    end
    
    % Viewer 3 (mean + 3*sigma) -----------------------------------------------
    v3 = viewer(PCplusMorph);
    
    % Set white background
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'Color', [1 1 1]);
    
    % title('PCplusMorph');
    v3.SceneLightLinked = true;
    v3.SceneLightVisible = true;
    
    % Export as .fig
    % se8 = ".fig";
    % sges= strcat(se0,se1,se2,se3,se4,se5,se6,se7,se8);
    % savefig(gcf, sges);
    
    % Export as .png
    for i = 1:length(angles)
        figHandle = v3.Figure;
        figHandle.Color = [1 1 1]; % White background
        ax = findall(figHandle, 'Type', 'Axes');
        set(ax, 'Color', [1 1 1]);
    
        % Set view angles explicitly
        view(angles{i});
        camup(camups{i});
    
        % Construct a clear filename for export
        filename = sprintf('%s%s_%s_PC%d_plus_%s_%s.png', se0, se1, se3, PCnum, GPA, view_names{i});
    
        % Export as PNG
        exportgraphics(figHandle, filename, 'Resolution', 300, 'BackgroundColor', 'current');
    end
    
    % %% Standard GPA
    % 
    % SyntheticData3D_standardGPA = SyntheticData3D;
    % 
    % GPA = "standardGPA"; % for export naming
    % 
    % Template = shape3D;
    % Template.Vertices = SyntheticData3D(:,:,1);
    % 
    % % Center each shape, record its centroid‐size, and normalize to unit‐size
    % nShapes = size(SyntheticData3D_standardGPA,3);
    % centroid_sizes = zeros(nShapes,1);
    % 
    % for f = 1:nShapes
    %     pts = SyntheticData3D_standardGPA(:,:,f);
    %     c = mean(pts,1);
    %     centered = pts - c;    
    %     centroid_sizes(f) = norm(centered,'fro');
    %     % normalize to unit‐size for template selection & initial GPA
    %     SyntheticData3D_standardGPA(:,:,f) = centered / centroid_sizes(f);
    % end
    % 
    % % Compute the crude mean (unit‐size)
    % crudeMean = mean(SyntheticData3D_standardGPA,3);
    % 
    % % Pick the shape closest (in Frobenius distance) to that mean
    % dists = zeros(nShapes,1);
    % for f = 1:nShapes
    %     diff = SyntheticData3D_standardGPA(:,:,f) - crudeMean;
    %     dists(f) = norm(diff,'fro')^2;
    % end
    % [~, bestIdx] = min(dists);
    % 
    % % log which file is your starting template
    % templateFileName = syntheticNames{bestIdx};
    % % templateFileName = keptFileNames{bestIdx};
    % fprintf('Using %s as the starting template (unit‐normalized).\n', templateFileName);
    % 
    % % Initialize meanPoints from that selected unit‐size shape
    % meanPoints = SyntheticData3D_standardGPA(:,:,bestIdx);
    % 
    % % Iterative GPA, aligning & re‐scaling each to C_ref every iteration
    % iter  = 10;   
    % for i = 1:iter
    %     for f = 1:nShapes
    %         % align to current meanPoints (which is at size = C_ref)
    %         faceVerts = SyntheticData3D_standardGPA(:,:,f);           
    %         T = computeTransform(faceVerts, meanPoints, true);
    %         faceVerts = applyTransform(faceVerts, T);
    % 
    %         % re‐center (no scaling)
    %         centered = faceVerts - mean(faceVerts,1);
    %         SyntheticData3D_standardGPA(:,:,f) = centered;
    %     end
    % 
    %     % update the mean shape and re‐enforce size = C_ref
    %     meanPoints = mean(SyntheticData3D_standardGPA,3);
    %     meanPoints = meanPoints / norm(meanPoints,'fro');
    % end

    %% Standard GPA (neue Version 01.08.2025)

    SyntheticData3D_standardGPA = SyntheticData3D;

    GPA = "standardGPA"; % for export naming

    % % Temporäre Zwischenvariablen --> später löschen
    meanPoints = mean(SyntheticData3D_standardGPA, 3); % mean(DataMatrix3D,3);
    % nShapes = size(DataMatrix3D,3);
    nShapes = size(SyntheticData3D_standardGPA,3);
    % % meanShapeVec = mean(arrayToRowVector(DataMatrix3D), 1); 
    % % sharedMean = rowVectorToArray(meanShapeVec);  % N × 3

    % Parameters
    maxIter = 10;       % same number of iterations as before
    tol     = 1e-6;     % convergence threshold
    
    AlignedData = zeros(size(SyntheticData3D_standardGPA));
    
    for iterGPA = 1:maxIter
        % Align each shape to the current meanPoints via Procrustes
        for f = 1:nShapes
            pts = SyntheticData3D_standardGPA(:,:,f);
            % [~, Z, ~] returns the aligned points Z
            [~, Z, ~] = procrustes(meanPoints, pts, ...
                                   'scaling',    true, ...
                                   'reflection', false);
            AlignedData(:,:,f) = Z;
        end
    
        % Compute updated mean shape
        Mnew = mean(AlignedData, 3);
        Mnew = Mnew - mean(Mnew, 1);        % re‐center
        Mnew = Mnew / norm(Mnew, 'fro');    % unit‐normalize
    
        % Check convergence
        if norm(Mnew - meanPoints, 'fro') < tol
            break;
        end
    
        % Update for next iteration
        meanPoints = Mnew;
        SyntheticData3D_standardGPA = AlignedData;
    end

    % sharedMean = meanPoints;
    meanShapeVec = mean(arrayToRowVector(SyntheticData3D_standardGPA), 1); 
    sharedMean = rowVectorToArray(meanShapeVec);  % N × 3
    
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
    % 
    % % Save the figure
    % output_dir = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK';
    % fig_name = 'standardGPA_foot_overlay';
    % 
    % % Get current figure handle
    % hFig = gcf;
    % 
    % % Save as high-res .png
    % exportgraphics(hFig, fullfile(output_dir, [fig_name '.png']), 'Resolution', 300, 'BackgroundColor', 'white');

    % %% Test
    % X = reshape(SyntheticData3D_standardGPA, [], nShapes)';
    % [coeff, ~, ~] = pca(X);
    % 
    % PC1 = reshape(coeff(:,1), [], 3); % N × 3
    % 
    % % Compute average PC1 direction vector
    % v_pc1 = mean(PC1, 1);
    % v_pc1 = v_pc1 / norm(v_pc1); % normalize
    % disp('First PC average direction:');
    % disp(v_pc1);
    % 
    % % Clean foot-length mode (already normalized in your code)
    % footLengthVec_clean = arrayToRowVector(footLengthVec); 
    % footLengthVec_clean = footLengthVec_clean / norm(footLengthVec_clean);
    % footLengthVec_clean = reshape(footLengthVec_clean, [], 3);
    % 
    % % Compute average direction of the clean foot-length mode
    % v_clean = mean(footLengthVec_clean, 1);
    % v_clean = v_clean / norm(v_clean);
    % 
    % cos_sim = dot(v_pc1, v_clean);
    % fprintf('Cosine similarity between PC1 and clean foot-length mode: %.4f\n', cos_sim);
    % 
    % angle_deg = acosd(cos_sim);
    % fprintf('Angle between PC1 and clean foot-length mode: %.2f degrees\n', angle_deg);
    % 
    % % Visualize both vectors in 3D
    % meanShape = mean(SyntheticData3D_standardGPA, 3);  % Mean foot shape
    % centroid = mean(meanShape,1);                      % Origin for arrows
    % arrowLength = 40;                                  % Adjust length if needed
    % 
    % figure;
    % scatter3(meanShape(:,1), meanShape(:,2), meanShape(:,3), 10, 'k'); hold on;
    % quiver3(centroid(1), centroid(2), centroid(3), ...
    %         arrowLength*v_pc1(1), arrowLength*v_pc1(2), arrowLength*v_pc1(3), ...
    %         'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    % quiver3(centroid(1), centroid(2), centroid(3), ...
    %         arrowLength*v_clean(1), arrowLength*v_clean(2), arrowLength*v_clean(3), ...
    %         'b', 'LineWidth', 2, 'MaxHeadSize', 2);
    % 
    % legend('Mean Shape', 'PC1 Direction', 'Clean Foot Length Mode');
    % title(sprintf('PC1 vs Clean Foot Length Mode (Angle = %.2f°)', angle_deg));
    % axis equal; grid on; view(3);

    % %% === Restore original scale of standard GPA shapes ===
    % SyntheticData3D_standardGPA_rescaled = SyntheticData3D_standardGPA;
    % 
    % for f = 1:nShapes
    %     SyntheticData3D_standardGPA_rescaled(:,:,f) = SyntheticData3D_standardGPA(:,:,f) * centroid_sizes(f);
    % end

    %% Principal Components Analysis
    
    % Represent each foot as a row of a 2D matrix
    Y = arrayStructure2Vector(SyntheticData3D_standardGPA);
    
    % Y is already mean-centered data (Skript DK)!
    [coeff_sGPA, score_sGPA, latent_sGPA] = pca(Y); % sGPA := standard GPA
    
    % Percentage of total variance explained by each PC
    explained_sGPA = latent_sGPA / sum(latent_sGPA) * 100;

    %% Create and store shape deformation plots
    % PCnum = 1;                % first principal component
    % k     = 3;                % ±3σ
    % variation_idx = s;        % whichever variation you’re visualizing
    % 
    % exportPath = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\visualize_PCs_sample_size';
    % 
    % visualizePCMorph(coeff_sGPA, latent_sGPA, PCnum, sharedMean, k, ...
    %              'standardGPA', temporary_number_objects, variation_idx, ...
    %              exportPath);

    %% Visualise principal component by transforming average face

    PCnum = 1;
    k = 3; % Scaling factor

    Transform = arrayVector2Structure(coeff_sGPA(:,PCnum)');
    sigma = sqrt(latent_sGPA(PCnum));
    
    AverageShape = shape3D;
    AverageShape.Vertices = sharedMean;
    
    % Apply PC to visualize -3σ and +3σ deformations
    PCminusMorph = clone(AverageShape);
    PCminusMorph.Vertices = sharedMean + Transform*sigma*-k;
    
    PCmeanMorph = clone(AverageShape);
    PCmeanMorph.Vertices = sharedMean;
    
    PCplusMorph = clone(AverageShape);
    PCplusMorph.Vertices = sharedMean + Transform*sigma*k;
    
    % Set the mesh color to black
    PCminusMorph.SingleColor = [0 0 0];
    PCmeanMorph.SingleColor = [0 0 0];
    PCplusMorph.SingleColor = [0 0 0];
    
    % Create three plots (-3*sigma, mean, 3*sigma). Plot each in a separate
    % viewer.
    
    % Viewer 1 (mean - 3*sigma) -----------------------------------------------
    v1 = viewer(PCminusMorph);
    set(gcf, 'Color', [1 1 1]); % Set white background
    set(gca, 'Color', [1 1 1]); % Set white background
    % title('PCminusMorph');
    v1.SceneLightLinked = true;
    v1.SceneLightVisible = true;
    
    % Export as .fig
    cd 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\visualize_PCs';
    
    se0 = "n";
    se1 = num2str(temporary_number_objects);
    se2 = "_";
    se3 = num2str(variation_idx);
    se4 = "_PC";
    se5 = num2str(PCnum);
    se6 = "_minus_morph_";
    se7 = GPA;
    se8 = ".fig";
    sges= strcat(se0,se1,se2,se3,se4,se5,se6,se7,se8);
    
    % savefig(gcf, sges);
    
    % Export as .png
    
    % Variante 1 (mit top, bottom, und medial view)
    % Define the specific desired views
    angles = {
        [0, 1, 0];  % Top view (bird perspective, transversal plane)
        [0, 0, 1];  % Medial side view (sagittal plane, medial view)
        [0, -1, 0];  % Bottom view (observer below the foot)
    };
    
    % Corresponding camera-up vectors
    camups = {
        [0, 0, 0];  % top view: heel bottom, toes top
        [0, 1, 0];  % medial side view: vertical axis up
        [0, 0, 0];  % bottom view: heel top, toes bottom
    };
    
    view_names = {'top', 'medial', 'bottom'};
    
    for i = 1:length(angles)
        figHandle = v1.Figure;
        figHandle.Color = [1 1 1]; % White background
        ax = findall(figHandle, 'Type', 'Axes');
        set(ax, 'Color', [1 1 1]);
    
        % Set view angles explicitly
        view(angles{i});
        camup(camups{i});
    
        % Construct a clear filename for export
        filename = sprintf('%s%s_%s_PC%d_minus_%s_%s.png', se0, se1, se3, PCnum, GPA, view_names{i});
    
        % Export as PNG
        exportgraphics(figHandle, filename, 'Resolution', 300, 'BackgroundColor', 'current');
    end
    
    
    % Viewer 2 (mean) ---------------------------------------------------------
    v2 = viewer(PCmeanMorph);
    
    % Set white background
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'Color', [1 1 1]);
    
    % title('PCmeanMorph');
    v2.SceneLightLinked = true;
    v2.SceneLightVisible = true;
    
    se8 = ".fig";
    sges= strcat(se0,se1,se2,se3,se4,se5,se6,se7,se8);
    
    % Export as .fig
    savefig(gcf, sges);
    
    % Export as .png
    for i = 1:length(angles)
        figHandle = v2.Figure;
        figHandle.Color = [1 1 1]; % White background
        ax = findall(figHandle, 'Type', 'Axes');
        set(ax, 'Color', [1 1 1]);
    
        % Set view angles explicitly
        view(angles{i});
        camup(camups{i});
    
        % Construct a clear filename for export
        filename = sprintf('%s%s_%s_PC%d_mean_%s_%s.png', se0, se1, se3, PCnum, GPA, view_names{i});
    
        % Export as PNG
        exportgraphics(figHandle, filename, 'Resolution', 300, 'BackgroundColor', 'current');
    end
    
    % Viewer 3 (mean + 3*sigma) -----------------------------------------------
    v3 = viewer(PCplusMorph);
    
    % Set white background
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'Color', [1 1 1]);
    
    % title('PCplusMorph');
    v3.SceneLightLinked = true;
    v3.SceneLightVisible = true;
    
    % Export as .fig
    % se8 = ".fig";
    % sges= strcat(se0,se1,se2,se3,se4,se5,se6,se7,se8);
    % savefig(gcf, sges);
    
    % Export as .png
    for i = 1:length(angles)
        figHandle = v3.Figure;
        figHandle.Color = [1 1 1]; % White background
        ax = findall(figHandle, 'Type', 'Axes');
        set(ax, 'Color', [1 1 1]);
    
        % Set view angles explicitly
        view(angles{i});
        camup(camups{i});
    
        % Construct a clear filename for export
        filename = sprintf('%s%s_%s_PC%d_plus_%s_%s.png', se0, se1, se3, PCnum, GPA, view_names{i});
    
        % Export as PNG
        exportgraphics(figHandle, filename, 'Resolution', 300, 'BackgroundColor', 'current');
    end
    
    %% Constrained Generalized Procrustes Analysis (DK)
    
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

    templateFileName = syntheticNames{bestIdx};
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
    
    % sharedMean = meanPoints;
    meanShapeVec = mean(arrayToRowVector(SyntheticData3D_constrainedGPA), 1); 
    sharedMean = rowVectorToArray(meanShapeVec);  % N × 3

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
    % 
    % % Save the figure
    % output_dir = 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK';
    % fig_name = 'constrainedGPA_foot_overlay';
    % 
    % % Get current figure handle
    % hFig = gcf;
    % 
    % % Save as high-res .png
    % exportgraphics(hFig, fullfile(output_dir, [fig_name '.png']), 'Resolution', 300, 'BackgroundColor', 'white');

    %% Principal Components Analysis
    
    % Represent each foot as a row of a 2D matrix
    Y = arrayStructure2Vector(SyntheticData3D_constrainedGPA);
    
    % Y is already mean-centered data (Skript DK)!
    [coeff_cGPA, score_cGPA, latent_cGPA] = pca(Y); % cGPA := constrained GPA
    
    % Percentage of total variance explained by each PCverify if PC1 really captures foot length
    explained_cGPA = latent_cGPA / sum(latent_cGPA) * 100;

    %% === Visualize one overlay: rescaled standard GPA vs. constrained GPA ===
    % f = 1; % Foot index to visualize
    % 
    % pts_con = SyntheticData3D_constrainedGPA(:,:,f);
    % pts_std = SyntheticData3D_standardGPA_rescaled(:,:,f);
    % 
    % figure;
    % hold on; axis equal;
    % title(sprintf('Foot %d overlay: Constrained GPA (blue) vs. Standard GPA (rescaled, red)', f));
    % 
    % plot3(pts_con(:,1), pts_con(:,2), pts_con(:,3), 'b.');
    % plot3(pts_std(:,1), pts_std(:,2), pts_std(:,3), 'r.');
    % 
    % xlabel('X'); ylabel('Y'); zlabel('Z');
    % grid on; view(3);
    % legend('Constrained GPA', 'Standard GPA (rescaled)');

    % rmse = sqrt(mean(sum((pts_con - pts_std).^2, 2)));
    % fprintf('RMSE between constrained and rescaled standard GPA shapes (foot %d): %.4f mm\n', f, rmse);
    
    %% Visualise principal component by transforming average face
    
    Transform = arrayVector2Structure(coeff_cGPA(:,PCnum)');
    sigma = sqrt(latent_cGPA(PCnum));
    
    AverageShape = shape3D;
    AverageShape.Vertices = meanPoints;
    
    % add transform scaled according to  z-score of +/- 3 onto average vertices
    % to visualise PC
    PCminusMorph = clone(AverageShape); % clone makes a deep copy of the mean shape
    PCminusMorph.Vertices = sharedMean + Transform*sigma*-k;
    
    PCmeanMorph = clone(AverageShape);
    PCmeanMorph.Vertices = sharedMean;

    PCplusMorph = clone(AverageShape);
    PCplusMorph.Vertices = sharedMean + Transform*sigma*k;
    
    % Set the mesh color to black
    PCminusMorph.SingleColor = [0 0 0];
    PCmeanMorph.SingleColor = [0 0 0];
    PCplusMorph.SingleColor = [0 0 0];
    
    % Create three plots (-3*sigma, mean, 3*sigma). Plot each in a separate
    % viewer.
    
    % Viewer 1 (mean - 3*sigma) -----------------------------------------------
    v1 = viewer(PCminusMorph);
    set(gcf, 'Color', [1 1 1]); % Set white background
    set(gca, 'Color', [1 1 1]); % Set white background
    % title('PCminusMorph');
    v1.SceneLightLinked = true;
    v1.SceneLightVisible = true;
    
    % Export as .fig
    cd 'E:\2025_SSM_ArmaSuisse\Skripte\results_DK\visualize_PCs';
    
    se0 = "n";
    se1 = num2str(temporary_number_objects);
    se2 = "_";
    se3 = num2str(variation_idx);
    se4 = "_PC";
    se5 = num2str(PCnum);
    se6 = "_minus_morph_";
    se7 = GPA;
    se8 = ".fig";
    sges= strcat(se0,se1,se2,se3,se4,se5,se6,se7,se8);
    
    % savefig(gcf, sges);
    
    % Export as .png
    
    % Variante 1 (mit top, bottom, und medial view)
    % Define the specific desired views
    angles = {
        [0, 1, 0];  % Top view (bird perspective, transversal plane)
        [0, 0, 1];  % Medial side view (sagittal plane, medial view)
        [0, -1, 0];  % Bottom view (observer below the foot)
    };
    
    % Corresponding camera-up vectors
    camups = {
        [0, 0, 0];  % top view: heel bottom, toes top
        [0, 1, 0];  % medial side view: vertical axis up
        [0, 0, 0];  % bottom view: heel top, toes bottom
    };
    
    view_names = {'top', 'medial', 'bottom'};
    
    for i = 1:length(angles)
        figHandle = v1.Figure;
        figHandle.Color = [1 1 1]; % White background
        ax = findall(figHandle, 'Type', 'Axes');
        set(ax, 'Color', [1 1 1]);
    
        % Set view angles explicitly
        view(angles{i});
        camup(camups{i});
    
        % Construct a clear filename for export
        filename = sprintf('%s%s_%s_PC%d_minus_%s_%s.png', se0, se1, se3, PCnum, GPA, view_names{i});
    
        % Export as PNG
        exportgraphics(figHandle, filename, 'Resolution', 300, 'BackgroundColor', 'current');
    end
    
    
    % Viewer 2 (mean) ---------------------------------------------------------
    v2 = viewer(PCmeanMorph);
    
    % Set white background
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'Color', [1 1 1]);
    
    % title('PCmeanMorph');
    v2.SceneLightLinked = true;
    v2.SceneLightVisible = true;
    
    se8 = ".fig";
    sges= strcat(se0,se1,se2,se3,se4,se5,se6,se7,se8);
    
    % Export as .fig
    savefig(gcf, sges);
    
    % Export as .png
    for i = 1:length(angles)
        figHandle = v2.Figure;
        figHandle.Color = [1 1 1]; % White background
        ax = findall(figHandle, 'Type', 'Axes');
        set(ax, 'Color', [1 1 1]);
    
        % Set view angles explicitly
        view(angles{i});
        camup(camups{i});
    
        % Construct a clear filename for export
        filename = sprintf('%s%s_%s_PC%d_mean_%s_%s.png', se0, se1, se3, PCnum, GPA, view_names{i});
    
        % Export as PNG
        exportgraphics(figHandle, filename, 'Resolution', 300, 'BackgroundColor', 'current');
    end
    
    % Viewer 3 (mean + 3*sigma) -----------------------------------------------
    v3 = viewer(PCplusMorph);
    
    % Set white background
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'Color', [1 1 1]);
    
    % title('PCplusMorph');
    v3.SceneLightLinked = true;
    v3.SceneLightVisible = true;
    
    % Export as .fig
    % se8 = ".fig";
    % sges= strcat(se0,se1,se2,se3,se4,se5,se6,se7,se8);
    % savefig(gcf, sges);
    
    % Export as .png
    for i = 1:length(angles)
        figHandle = v3.Figure;
        figHandle.Color = [1 1 1]; % White background
        ax = findall(figHandle, 'Type', 'Axes');
        set(ax, 'Color', [1 1 1]);
    
        % Set view angles explicitly
        view(angles{i});
        camup(camups{i});
    
        % Construct a clear filename for export
        filename = sprintf('%s%s_%s_PC%d_plus_%s_%s.png', se0, se1, se3, PCnum, GPA, view_names{i});
    
        % Export as PNG
        exportgraphics(figHandle, filename, 'Resolution', 300, 'BackgroundColor', 'current');
    end

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

    % Force v_con into the same hemisphere as v_gt
    if dot(v_con, v_gt) < 0
        v_con = -v_con;
    end
    
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
    variation_val = foot_length_levels(variation_idx);
    figName = sprintf('meanShape_PC1_variation_%d', variation_val);

    % Save as MATLAB figure
    savefig(hFig, fullfile(output_dir, [figName '.fig']));

    % Save as PNG
    saveas(hFig, fullfile(output_dir, [figName '.png']));

    close(hFig);

    %% ------- EVALUATION METRICS -----------------------------------------

    % Precompute means for Procrustes once
    mean_ground_truth = rowVectorToArray(meanShapeVec);     % original mean
    mean_std = mean(SyntheticData3D_standardGPA,3);          % std‐GPA mean
    mean_con = mean(SyntheticData3D_constrainedGPA,3);       % con‐GPA mean

    % Compute ground-truth explained variance in SHAPE space
    v_weights  = var(W_full);                  % 1×5
    v_pcnorms  = vecnorm(full_pc_basis).^2;    % 1×5  ← match v_weights
    v_shape    = v_weights .* v_pcnorms;       % 1×5
    v_shape_total           = sum(v_shape);
    var_ground_truth_percent = v_shape / v_shape_total * 100;

    fprintf('\n=== Recovery Metrics for variation L=%g ===\n', foot_length_levels(variation_idx));
    for k = 1:nPCsToUse
        % ground‐truth PC
        ground_truth_k      = ground_truth_pcs(:,k);
        ground_truth_k      = ground_truth_k / norm(ground_truth_k);

        % recovered PC
        est_std_k = coeff_sGPA(:,k);
        est_std_k = est_std_k / norm(est_std_k);
        est_con_k = coeff_cGPA(:,k);
        est_con_k = est_con_k / norm(est_con_k);

        % Flip sign if necessary to align with ground-truth
        if dot(ground_truth_k, est_std_k) < 0
            est_std_k = -est_std_k;
        end
        cosine_similarity_std_k = dot(ground_truth_k, est_std_k);

        if dot(ground_truth_k, est_con_k) < 0
            est_con_k = -est_con_k;
        end
        cosine_similarity_con_k = dot(ground_truth_k, est_con_k);

        % cosine similarity (absolute)
        % cosine_similarity_std_k = abs(dot(ground_truth_k, est_std_k));
        % cosine_similarity_con_k = abs(dot(ground_truth_k, est_con_k));

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
        variation_val = foot_length_levels(variation_idx);
        csv_filename = fullfile(output_dir, sprintf('results_variation_%d.csv', variation_val));
        writetable(results_table, csv_filename);
    end

    close all;

end % Ende der foot_length_levels Schleife










