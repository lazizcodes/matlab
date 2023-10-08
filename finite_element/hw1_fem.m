% Define problem parameters
Lx = 1;   % Domain length in x-direction
Ly = 1;   % Domain length in y-direction
Nx = 100;  % Number of elements in x-direction
Ny = 100;  % Number of elements in y-direction

% Generate mesh
[x, y] = meshgrid(linspace(0, Lx, Nx+1), linspace(0, Ly, Ny+1));

% Define boundary conditions
% Dirichlet boundary conditions
u_x0_y = zeros(Ny+1, 1);   % u(x, 0) = 0
u_0_y = zeros(Ny+1, 1);    % u(0, y) = 0

% Neumann boundary conditions
% (Note: Neumann BCs are applied to the stiffness matrix directly)

% Define the source term
f = @(x, y) 400 * (x.^4 - y.^4) .* sin(20 * x .* y);

% Initialize the stiffness matrix and load vector
K = zeros((Nx+1)*(Ny+1), (Nx+1)*(Ny+1));
F = zeros((Nx+1)*(Ny+1), 1);

% Loop over elements and assemble the global stiffness matrix and load vector
for i = 1:Nx
    for j = 1:Ny
        % Define local element coordinates
        x1 = x(i, j);
        x2 = x(i+1, j);
        x3 = x(i, j+1);
        x4 = x(i+1, j+1);
        y1 = y(i, j);
        y2 = y(i+1, j);
        y3 = y(i, j+1);
        y4 = y(i+1, j+1);

        % Compute element stiffness matrix and load vector
        [Ke, Fe] = computeElementMatrices(x1, x2, x3, x4, y1, y2, y3, y4, f);

        % Assemble into global stiffness matrix and load vector
        % (using element connectivity)
        indices = elementConnectivity(i, j, Nx, Ny);
        K(indices, indices) = K(indices, indices) + Ke;
        F(indices) = F(indices) + Fe;
    end
end

% Apply Neumann boundary conditions directly to the load vector
x_values = x(end, :);
y_values = y(:, end);
F = applyNeumannBCs(F, x_values, y_values, Nx, Ny);

% Apply Dirichlet boundary conditions
K(1:Ny+1, :) = 0;
K(1:Ny+1, 1:Ny+1) = eye(Ny+1);  % u(0, y) = 0
F(1:Ny+1) = u_0_y;

% Solve the linear system of equations
u = K\F;

% Reshape the solution to a 2D grid
u_grid = reshape(u, Ny+1, Nx+1);

% Visualize the solution
contourf(x, y, u_grid, 20);
colorbar;
xlabel('x');
ylabel('y');
title('Numerical Solution of PDE');
