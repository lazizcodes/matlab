% ### SOLUTION ###

[tri, x, y] = generateTriangularMesh(60);
trimesh(tri, x, y); axis('equal');

% Set up parameters
numNodes = numel(x);
numTriangles = size(tri, 1);
dt = 0.0005; % time step
finalTime = 1; % final simulation time
programLifeTime = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

% Initialize solution
u = zeros(numNodes, 1);

% Set initial condition inside the domain
for i = 1:numNodes
    u(i) = sin(pi * (x(i)^2 + y(i)^2)) * ((x(i) - 1)^2 + y(i)^2 - 9);
end

% Enforce boundary condition (u = 0 on the boundary)
u = enforceBoundaryCondition(u, tri, x, y);

% Velocity field
v = [1, 2];

% Create the initial figure
figure;
h = trisurf(tri, x, y, u, 'EdgeColor', 'k', 'FaceColor', 'interp');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Solution (u)');
title(['Numerical Solution at t = 0']);
drawnow;

% Time-stepping loop
for t = dt:dt:finalTime
    % Laplacian term
    laplacian_u = computeLaplacian(tri, x, y, u);
    
    % Loop over triangles for FVM
    for i = 1:numTriangles
        nodes = tri(i, :);
        
        % Calculate face normals
        normal = computeFaceNormals(x(nodes), y(nodes));
        
        % Compute flux term
        flux = (v * normal') / norm(normal) * u(nodes);

        
        % Update solution using FVM
        u(nodes) = u(nodes) + (dt / polyarea(x(nodes), y(nodes))) * (0.1 * laplacian_u(i) - flux);
    end

    % Enforce boundary condition (u = 0 on the boundary)
    u = enforceBoundaryCondition(u, tri, x, y);

    if ismember(t, finalTime * programLifeTime)
        % Update the figure
        set(h, 'Vertices', [x(:), y(:), u(:)]);
        title(['Numerical Solution at t = ' num2str(t)]);
        drawnow;
    end
end

% ### HELPER FUNCTIONS ###

function u = enforceBoundaryCondition(u, tri, x, y)
    numNodes = numel(x);
    numTriangles = size(tri, 1);

    % Loop over triangles to identify boundary nodes
    for i = 1:numTriangles
        nodes = tri(i, :);
        
        % Identify nodes on the boundary
        boundaryNodes = nodes(x(nodes).^2 + y(nodes).^2 == 1 | ...
                               (x(nodes)-1).^2 + y(nodes).^2 == 9);
        
        % Set solution to zero on boundary nodes
        u(boundaryNodes) = 0;
    end
end

function laplacian_u = computeLaplacian(tri, x, y, u)
    % Compute Laplacian of the solution at each triangle using linear shape functions
    
    numTriangles = size(tri, 1);
    laplacian_u = zeros(numTriangles, 1);

    for i = 1:numTriangles
        nodes = tri(i, :);
        area = polyarea(x(nodes), y(nodes));
        
        % Linear shape functions for the Laplacian
        N1 = 1 / (2 * area) * (y(nodes(2)) - y(nodes(3)));
        N2 = 1 / (2 * area) * (y(nodes(3)) - y(nodes(1)));
        N3 = 1 / (2 * area) * (y(nodes(1)) - y(nodes(2)));
        
        % Calculate Laplacian using linear shape functions
        laplacian_u(i) = N1 * u(nodes(1)) + N2 * u(nodes(2)) + N3 * u(nodes(3));
    end
end

function faceNormals = computeFaceNormals(x, y)
    % Compute outward-facing normal vectors for each face of a triangle
    
    numNodes = length(x);
    faceNormals = zeros(numNodes, 2);

    for i = 1:numNodes
        j = mod(i, numNodes) + 1; % Next node in the sequence
        dx = x(j) - x(i);
        dy = y(j) - y(i);
        
        % Compute outward-facing normal vector
        faceNormals(i, :) = [-dy, dx] / norm([dx, dy]);
    end
end


function centroids = getTriangularMeshCentroids(tri, x, y)
centroids = zeros(size(tri, 1), 2);

% Calculate centroids
for i = 1:size(tri, 1)
    % Extract vertices of the current triangle
    x_coors = x(tri(i, :));
    y_coors = y(tri(i, :));

    % Calculate centroid as the average of the vertices
    centroids(i, :) = [mean(x_coors), mean(y_coors)];
end
end

function [tri, x, y] = generateTriangularMesh(N)
delta = 6 / N;
counter = 1;

% Generate initial grid
for i = 1:N
    for j = 1:N
        x0(counter) = -2 + i * delta;
        y0(counter) = -3 + j * delta;
        counter = counter + 1;
    end
end

counter = 1;

% Extract in-domain nodes
for i = 1:N*N
    if norm([x0(i)-1 y0(i)]) < 3-delta/2 && norm([x0(i) y0(i)]) > 1+delta/2
        x(counter) = x0(i);
        y(counter) = y0(i);
        counter = counter + 1;
    end
end

counter = counter - 1;

% Populate boundary nodes
bnodes_out = 3 * N;
for i = 1:bnodes_out
    h = 1 + 3 * cos(2 * pi / bnodes_out * i);
    v = 3 * sin(2 * pi / bnodes_out * i);
    x(i+counter) = h; y(i+counter) = v;
end

bnodes_in = N;
for i = 1:bnodes_in
    h = cos(2 * pi / bnodes_in * i);
    v = sin(2 * pi / bnodes_in * i);
    x(i+bnodes_out+counter) = h; y(i+bnodes_out+counter) = v;
end

% Add a central node
x(bnodes_in+bnodes_out+counter+1) = 0;
y(bnodes_in+bnodes_out+counter+1) = 0;

% Perform Delaunay triangulation
tri = delaunay(x, y);

% Remove elements connected with dummy central node
[numRows, ~] = size(tri);
elemsToRemove = [];

for k = 1:numRows
    if ismember(bnodes_in + bnodes_out + counter + 1, tri(k, :))
        elemsToRemove = [elemsToRemove, k];
    end
end

tri = tri(setdiff(1:numRows,elemsToRemove),:);

x = x(1:end-1);
y = y(1:end-1);
end