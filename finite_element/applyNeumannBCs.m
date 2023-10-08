function F = applyNeumannBCs(F, x_values, y_values, Nx, Ny)
% Apply Neumann boundary conditions to the load vector (F)

% Neumann boundary condition on the top boundary (y = 1)
for i = 1:Nx
    x_mid = (x_values(i) + x_values(i+1)) / 2;
    F(end - Nx + i) = F(end - Nx + i) + 20 * x_mid * (x_mid^2 - 1) * cos(20 * x_mid) - 2 * sin(20 * x_mid);
end

% Neumann boundary condition on the right boundary (x = 1)
for j = 1:Ny
    y_mid = (y_values(j) + y_values(j+1)) / 2;
    F(j * (Nx + 1) + Nx + 1) = F(j * (Nx + 1) + Nx + 1) + 20 * y_mid * (1 - y_mid^2) * cos(20 * y_mid) + 2 * sin(20 * y_mid);
end
end