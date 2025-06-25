% Add MeshMonk's toolbox to the working path and setup current folder
clear all;
% close all;

addpath(genpath('C:\Users\algona_workstation\Documents\GitHub\meshmonk\meshmonk\')) % Set to location of meshmonk

studypath = 'C:\Users\algona_workstation\Documents\SSM_Arma\Daten_ArmaSuisse\processed_R\';   % Set to location of demo material
cd(studypath);

% List of OBJ files (foot scans)
fileList = {'fuss1.obj', 'fuss2.obj', 'fuss5.obj', 'fuss7.obj', ...
    'fuss9.obj', 'fuss11.obj', 'fuss13.obj', 'fuss17.obj', 'fuss19.obj', ...
    'fuss21.obj', 'fuss23.obj', 'fuss25.obj', 'fuss27.obj', 'fuss35.obj', ...
    'fuss37.obj', 'fuss39.obj', 'fuss41.obj', 'fuss45.obj', 'fuss49.obj', ...
    'fuss51.obj', 'fuss53.obj', 'fuss57.obj', 'fuss59.obj', 'fuss61.obj', ...
    'fuss63.obj', 'fuss65.obj', 'fuss67.obj', 'fuss73.obj', 'fuss75.obj', ...
    'fuss87.obj', 'fuss91.obj', 'fuss97.obj', 'fuss99.obj', 'fuss101.obj'};
path = pwd;

% Load Template Mesh (reference)
template = shape3D; % Erzeugt ein leeres shape3D Objekt % Als Template/Target wird aktuell fuss6 verwendet
importWavefront(template, 'template_foot.obj', path, []); % Lädt das eigentliche obj file

% Loop through each target foot scan
results = {};
for i = 1:length(fileList)
    fprintf('Processing target scan: %s\n', fileList{i});
    
    % Load Target Mesh
    target = shape3D;
    importWavefront(target, fileList{i}, path, []);

    % Initialize Shape Mapper
    obj = ShapeMapper;
    obj.FloatingShape = clone(template); % Template (Floating)
    obj.TargetShape = clone(target);     % Target shape (Scan)

    %% Rigid Registration
    obj.Display = false;
    obj.NumIterations = 80;
    obj.InlierKappa = 3;
    obj.TransformationType = 'rigid';
    obj.UseScaling = true;
    obj.CorrespondencesNumNeighbours = 3;
    obj.CorrespondencesFlagThreshold = 0.9;
    obj.CorrespondencesSymmetric = true;
    obj.InlierUseOrientation = true;
    obj.FlagFloatingBoundary = true;
    obj.FlagTargetBoundary = true;
    obj.FlagTargetBadlySizedTriangles = true;
    obj.TriangleSizeZscore = 6;
    obj.UpSampleTarget = false;

    % Perform Rigid Registration
    tic; map(obj); time = toc;
    fprintf('Rigid Registration took %f seconds.\n', time);

    %% Non-Rigid Registration
    obj.Display = false;
    obj.NumIterations = 200;
    obj.TransformNumNeighbors = 80;
    obj.CorrespondencesNumNeighbours = 3;
    obj.CorrespondencesSymmetric = true;
    obj.CorrespondencesFlagThreshold = 0.9;
    obj.CorrespondencesUseOrientation = true;
    obj.InlierKappa = 12;
    obj.FlagFloatingBoundary = true;
    obj.FlagTargetBoundary = true;
    obj.FlagTargetBadlySizedTriangles = true;
    obj.TriangleSizeZscore = 6;
    obj.UpSampleTarget = false;
    obj.UseScaling = 1;
    obj.TransformationType = 'nonrigid';
    obj.TransformSigma = 3;
    obj.TransformNumViscousIterationsStart = 200;
    obj.TransformNumViscousIterationsEnd = 1;
    obj.TransformNumElasticIterationsStart = 200;
    obj.TransformNumElasticIterationsEnd = 1;

    % Perform Non-Rigid Registration
    tic; map(obj); time = toc;
    fprintf('Non-Rigid Registration took %f seconds.\n', time);

    % Save the registered mesh (aligned to template)
    saveFile = sprintf('aligned_%s.mat', fileList{i});
    save(saveFile, 'obj');

    % % Visualize results (optional)
    % vref = viewer3D;
    % obj.TargetShape.ViewMode = 'points';
    % obj.FloatingShape.ViewMode = 'points';
    % obj.TargetShape.VertexSize = 10;
    % obj.FloatingShape.VertexSize = 12;
    % viewer(obj.FloatingShape, vref);
    % viewer(obj.TargetShape, vref);

    tmp = pointCloud(obj.FloatingShape.Vertices, color = [1 1 0]);

    results(i) = {tmp};
end

%% Plausbilitaetscheck
template_pt = pointCloud(template.Vertices, color = [1 1 0]);

% Indizes plotten
idx = randi([1,500],1,15); % Vektor mit Zufallsvariablen

color_plot = ["r","g","b","y","m","c","k","w","#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];

figure(1);
pcshow(template_pt);
hold on;
xlabel('X');
ylabel('Y');
zlabel('Z');
for j = 1:15
    pcshow(template_pt.Location(idx(1,j), :), color_plot(j), 'MarkerSize', 100);
end

for i=1:length(results)
    pcshow(results{1,i});
    for j=1:15
        pcshow(results{1,i}.Location(idx(1,j),:),color_plot(j),'Markersize',100)
    end
end