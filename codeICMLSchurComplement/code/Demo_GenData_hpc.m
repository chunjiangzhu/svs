%---------------------------
% Author : 
% Year   : 
%---------------------------
% CLUSTERING OPTIONS
%---------------------------

%-------------------------------------------------------------------------
% Generate the Graphs and save in T.txt and W.mat. 
%-------------------------------------------------------------------------


% Clear before start
clear;
seed = 15;

% Choose similarity graph type
% - 1 for full similarity graph
% - 2 for epsilon neighborhood
% - 3 for k-Nearest neighbors
% - 4 for mutual k-Nearest neighbors
% SimGraphType = 3;

% Choose spectral clustering algorithm
% - 1 for unnormalized Clustering
% - 2 for normalized   Clustering (Shi and Malik)
% - 3 for normalized   Clustering (Jordan and Weiss)
ClusteringType = 3;

% Choose a test mode
% - 1 for four Gaussians
% - 2 for two half moons
% - 3 for sculpture
% - 4 for two circles
% - 5 for camera
test_mode = 3;

% Choose whether to draw the similarity graph and if so,
% whether to color it (Only use with few data points!)
bplotSimGraph   = 0;
plotSimGraphMap = 1;

%---------------------------
% CREATE DATA
%---------------------------
time_data = tic;

if test_mode == 1
    
    % Four Gaussians
    T = [GenData_Gaussian(500, [0 0], 0.1*eye(2)) ...
         GenData_Gaussian(500, [2 0], 0.1*eye(2)) ...
         GenData_Gaussian(500, [1 1], 0.1*eye(2)) ...
         GenData_Gaussian(500, [1 -1], 0.1*eye(2))];
    
    SimGraphType = 4;
    num_clust = 4;
    Paramk    = 16;
    ParamEps  = 0.15;
    Sigma = 10;
    
elseif test_mode == 2
    
    % Two half moons
    T = [GenData_Ellipse(600, 1, 0.01, [0 0], 0, pi) ...
         GenData_Ellipse(600, 1, 0.01, [1 0.5], pi, 2*pi)];
    
    SimGraphType   = 2;
    num_clust = 2;
    Paramk    = 30;
    ParamEps  = 0.1; 
    Sigma = 0.5;
    
elseif test_mode == 4
    
    % 2 circles
    T = [GenData_Ellipse(1000, 1, 0, [0 0], 0, 2*pi) ...
         GenData_Ellipse(1000, 1.05, 0, [0 0], 0, 2*pi)];
    
    SimGraphType   = 4;
    num_clust = 2;
    Paramk    = 20;
    ParamEps  = 0.1; 
    Sigma = 20;
    
elseif test_mode == 3
    % Sculpture
        % Sculpture
    SimGraphType   = 4;
    
    I = imread('sculpture.jpg');
    a = size(I);
    width = size(I,1)
    height = size(I,2)
    num_clust = 3;
    Paramk    = 37; 
    ParamEps  = 5.5; 
    Sigma = 5;
    
    x_gap = 4;
    y_gap = 4;
    x_range = ceil(width/x_gap);
    y_range = ceil(height/y_gap);
    I_small(1:x_range,1:y_range,:) = I(1:x_gap:width,1:y_gap:height,:);
%     imshow(I_small);
    T = zeros(5,x_range*y_range);
    for i=1:x_range
        for j=1:y_range
            col = (i-1)*y_range+j;
            T(1,col) = i;
            T(2,col) = j;
            T(3,col) = I_small(i,j,1);
            T(4,col) = I_small(i,j,2);
            T(5,col) = I_small(i,j,3);
        end
    end
elseif test_mode == 5
    % tiger
    SimGraphType   = 2;
    I = imread('camera.jpg');
    a = size(I);
    width = size(I,1)
    height = size(I,2)
    num_clust = 3;
    Paramk    = 14; 
    ParamEps  = 5; 
    Sigma = 3.8;
    
    x_gap = 8;
    y_gap = 8;
    x_range = ceil(width/x_gap);
    y_range = ceil(height/y_gap);
    I_small(1:x_range,1:y_range,:) = I(1:x_gap:width,1:y_gap:height,:);
%     imshow(I_small);
    T = zeros(5,x_range*y_range);
    for i=1:x_range
        for j=1:y_range
            col = (i-1)*y_range+j;
            T(1,col) = i;
            T(2,col) = j;
            T(3,col) = I_small(i,j,1);
            T(4,col) = I_small(i,j,2);
            T(5,col) = I_small(i,j,3);
        end
    end
end

time_data = toc(time_data);
fprintf('Generating Data: %.2fs\n', time_data);

%---------------------------
% EXECUTE ALGORITHM AND PLOT
%---------------------------
 
% Parameter corresponding to SimGraphType (k or Eps)
switch SimGraphType
    case 1% the value for case 1 is arbitrary
        Param = 0;
    case 2
        Param = ParamEps;
    case {3, 4}
        Param = Paramk;
end

% calculate adjacency matrix
time_graph = tic;

W = SimGraph(T, SimGraphType, Param, Sigma);
% load W_circle
% sum(W1~=W)

m = nnz(W)/2
n = size(W,1)

% [x,y] = find(W);
% G = digraph(x,y);
% plot(G,'Layout','layered')

time_graph = toc(time_graph);
fprintf('Distance Matrix and Similarity Graph: %.2fs\n', ...
    time_graph);

n = size(W,1);


if test_mode == 1
    W_gaussian = W;
    save W_gaussian;
    dlmwrite('T_gaussian.txt',T);
elseif test_mode == 2
    W_halfmoon = W;
    save W_halfmoon;
    dlmwrite('T_halfmoon.txt',T);
elseif test_mode == 3
    save W_sculpture;
    dlmwrite('T_sculpture.txt',T);
elseif test_mode == 4
    W_circle = W;
    save W_circle;
    dlmwrite('T_circle.txt',T);
elseif test_mode == 5
%     W_circle = W;
    save W_camera;
    dlmwrite('T_camera.txt',T);
end
fprintf('Writed T\n');
fprintf('Writed sparse W\n');

% check for connected components (needs Bioinformatics Toolbox)
try
    conncomps = graphconncomp(W, 'Directed', false);
    if conncomps > 2
        warning(['Found %d connected components, maybe try ' ...
                 'changing the parameters.\n'], conncomps);
    else
        fprintf('Number of connected components: %d\n', ...
            conncomps);
    end
catch exception
    disp('Could not check for number of connected components.');
end

%compute clusters using spectral clustering
time_algorithm = tic;

A = SpectralClustering(W, num_clust, num_clust, ClusteringType);
time_algorithm = toc(time_algorithm);
fprintf('Clustering Algorithm: %.2fs\n', time_algorithm);

% plot result
scrsz = get(0,'ScreenSize');
hFig  = figure('Position', [scrsz(3)/10 scrsz(4)/3 ...
                            scrsz(3)/1.25 scrsz(4)/2]);
set(gca, 'Color', 'none');

if bplotSimGraph == 1
    subplot(1, 2, 1);
end
hold on;
col = [A zeros(size(A, 1), 1)];
% sFigureTitle = sprintf(...
%     'Clustered Data (Algorithm Time: %fs)', time_algorithm);
% title(sFigureTitle)
cols = 'brgcmykw';

for ii = 1:num_clust
    curcol = [cols(mod(ii, size(cols, 2))) '*'];
    scatter(T(1, A(:, ii) == 1), T(2, A(:, ii) == 1), 7, curcol);
end
hold off;

% plot similarity graph
if bplotSimGraph == 1
    hSG = subplot(1, 2, 2);
    title('Similarity Graph');
    plotSimGraph(T, W, plotSimGraphMap);
end

print('Output.pdf', '-dpdf', '-bestfit');
% Commands to export with export_fig
% export_fig 'output.pdf' -q101 -transparent
% export_fig(hSG, 'output.pdf', '-q101', '-transparent')
