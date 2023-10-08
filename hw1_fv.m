function [error, x, y] = hw1_fv(Nx, Ny)
    Lx = 1;   % Length of the domain in the x-direction
    Ly = 1;   % Length of the domain in the y-direction
    dx = Lx / Nx; % Grid spacing in the x-direction
    dy = Ly / Ny; % Grid spacing in the y-direction
    
    % Create grid
    x = linspace(0, Lx, Nx+1); % Including boundary points
    y = linspace(0, Ly, Ny+1); % Including boundary points
    
    % Initialize solution matrix
    u = zeros(Nx+1, Ny+1);
    u_old = zeros(Nx+1, Ny+1);
    
    % Define the right-hand side function
    f = @(x, y) 400 * (x.^4 - y.^4) .* sin(20 * x .* y);
    
    % Apply boundary conditions
    u(:, 1) = 0; % u(x,0) = 0
    u(1, :) = 0; % u(0,y) = 0
    
    k = 1;
    maxIterations = 1000;
    tolerance = 1e-3;
    
    while k <= maxIterations
        % Discretize the Laplacian operator using central differences
        for i = 2:Nx
            for j = 2:Ny
                % Calculate second derivatives using central differences
                u_xx = (u(i+1, j) - 2*u(i, j) + u(i-1, j)) / dx^2;
                u_yy = (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / dy^2;
    
                % Calculate source term
                source = integral2(f, x(i-1), x(i), y(j-1), y(j));
    
                % Update the solution
                u(i, j) = u(i, j) + (u_xx + u_yy - source) * dx^2 * dy^2;
            end
        end
    
        % Apply Neumann boundary conditions
        for i = 2:Nx
            source_x = 20 * x(i) * (x(i)^2 - 1) * cos(20 * x(i)) - 2 * sin(20 * x(i));
            u(i, Ny+1) = u(i, Ny) + source_x * dx;
        end
    
        for j = 2:Ny
            source_y = 20 * y(j) * (1 - y(j)^2) * cos(20 * y(j)) + 2 * sin(20 * y(j));
            u(Nx+1, j) = u(Nx, j) + source_y * dy;
        end
    
        if max(abs(u - u_old)) < tolerance
            break
        end
    
        u_old = u;
        k = k + 1;
    end
    
    
    
    % Define the exact solution function
    u_exact = @(x, y) (x.^2 - y.^2) .* sin(20 .* x .* y);
    
    % Compute the exact solution on the grid
    u_exact_values = u_exact(x, y);
    
    % Calculate the difference between the numerical and exact solutions
    error = abs(u - u_exact_values);
    
    % Plot the difference/error
    [X, Y] = meshgrid(x, y);
    surf(X, Y, error');
    xlabel('x');
    ylabel('y');
    zlabel('Error');
    title('Error Between Numerical and Exact Solutions');
end