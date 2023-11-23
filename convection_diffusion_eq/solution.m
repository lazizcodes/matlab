[tri, x, y] = generateTriangularMesh(30);
centroids = getTriangularMeshCentroids(tri, x, y);
trimesh(tri, x, y); axis('equal');

% Add centroids as points
hold on;
scatter(centroids(:, 1), centroids(:, 2), 5, 'r', 'filled');
hold off;

% Set up parameters
mu = 0.1;
numNodes = numel(x);
numTriangles = size(tri, 1);
dt = 0.01; % time step
finalTime = 1; % final simulation time

% Initialize solution
u = zeros(numTriangles, 1);

% Set initial condition inside the domain
for i = 1:numTriangles
    centroid = centroids(i, :);
    u(i) = sin(pi * (centroid(1)^2 + centroid(2)^2)) * ((centroid(1) - 1)^2 + centroid(2)^2 - 9);
end

% Velocity field
v = [1, 2];

% Create the initial figure
figure;
% h = trisurf(tri, centroids(:, 1), centroids(:, 2), u, 'EdgeColor', 'k', 'FaceColor', 'interp');
scatterObj = scatter(centroids(:, 1), centroids(:, 2), 70, u, 'filled');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Solution (u)');
title(['Numerical Solution at t = 0']);
colorbar;
drawnow;

t = dt;
while t < finalTime
    disp(t)
    % Loop over triangles for FVM
    for i = 1:numTriangles
        nodes = tri(i, :);

        % Calculate the centroid of the current triangle
        centroid = centroids(i, :);

        % Calculate the area of the current triangle
        triangleArea = polyarea(x(nodes), y(nodes));

        % Initialize the summation term
        summationTerm = 0;
        flag = 0;
        % Iterate over edges of the current triangle
        for j = 1:3
            edge = nodes([j, mod(j, 3) + 1]); % Current edge as [start, end]
            edgeLength = norm([x(edge(1)) - x(edge(2)), y(edge(1)) - y(edge(2))]);

            % Find the index of the adjacent triangle sharing the current edge
            % adjacentTriangle = find(sum(ismember(tri(:, 1:3), edge), 2) == 2 & ~ismember(1:numTriangles, i));
            adj = findAdjTri(tri, edge, i);

            if adj < 0
                flag = 1;
                break;
            end


            nrm = centroids(adj, :) - centroids(i, :);
            vij = dot(v, nrm / norm(nrm));

            % Accumulate the summation term
            summationTerm = summationTerm + edgeLength * (max(vij, 0) * u(i) + min(vij, 0) * u(adj) - mu * (u(adj) - u(i)) / norm(nrm));
        end

        if flag < 1
            % Update the solution using the given update rule
            u(i) = u(i) - (dt / triangleArea) * summationTerm;
        end
    end

    scatter(centroids(:, 1), centroids(:, 2), 70, u, 'filled');
    title(['Numerical Solution at t = ' num2str(t)]);
    colorbar
    drawnow;
    t = t + dt;
end

% ### HELPER FUNCTIONS ###

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

function matchingRow = findAdjTri(tri, edge, j)
% Initialize matching row as an empty array
matchingRow = -1;

target1 = edge(1);
target2 = edge(2);

% Iterate over each row
for i = 1:size(tri, 1)
    % Check if both targets are present in the current row
    if i ~= j && ismember(target1, tri(i, :)) && ismember(target2, tri(i, :))
        % Set matchingRow to the current row index
        matchingRow = i;
        % Break the loop since we found a match
        break;
    end
end
end