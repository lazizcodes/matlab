% Parameters
Lx = 1; % interval length
Ly = 1; % interval length
Nx = 64; % Number of grid points in the x-direction
Ny = 64; % Number of grid points in the y-direction
dx = Lx / (Nx - 1); % Grid spacing in the x-direction
dy = Ly / (Ny - 1); % Grid spacing in the y-direction

% Create grid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% Exact solution
u_exact = @(x, y) (x.^2 - y.^2) .* sin(20 .* x .* y);
u_exact_values = u_exact(X, Y);

% RHS
f = @(x, y) 400 * (x.^4 - y.^4) .* sin(20 * x .* y);

% Neumann Boundary Conditions
bc1 = @(x) 20 * x .* (x.^2 - 1) .* cos(20 * x) - 2 * sin(20 * x);
bc2 = @(y) 20 * y .* (1 - y.^2) .* cos(20 * y) + 2 * sin(20 * y);

% Initialize the solution matrix
u = zeros(Ny, Nx);

u(:, 1) = 0; % u(0, y) = 0
u(1, :) = 0; % u(x, 0) = 0

eplsilon = 1e-8;
maxIterations = 10000;
k = 0;

while k < maxIterations
    uPrev = u;
    
    % Update interior grid point value
    for i = 2:Ny-1
        for j = 2:Nx-1
            u(i, j) = (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) + dx * dy * f(x(j), y(i))) / 4;
        end
    end

    % Neumann boundary conditions
    % u_x at x = 1 (right boundary)
    for i = 2:Ny-1
        u(i, Nx) = u(i, Nx-1) + dx * bc2(y(i));
    end
    
    % u_y at y = 1 (top boundary)
    for j = 2:Nx-1
        u(Ny, j) = u(Ny-1, j) + dy * bc1(x(j));
    end 

    if max(abs(u - u_exact_values)) < eplsilon
        break;
    end
    
    k = k + 1;
end

fprintf('Number of iterations: %d\n', k);

%% Test


figure;
soln_diff = abs(u_exact_values - u);
surf(X, Y, soln_diff);
colorbar;
xlabel('x');
ylabel('y');
zlabel('|u - u_{exact}|')
title('Abs difference');
view(3);


