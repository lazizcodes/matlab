function [Ke, Fe] = computeElementMatrices(x1, x2, x3, x4, y1, y2, y3, y4, f)
% Compute the element stiffness matrix (Ke) and load vector (Fe)

% Compute the area of the element
area = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

% Define shape function derivatives
N1_x = (y2 - y4) / (2 * area);
N2_x = (y3 - y1) / (2 * area);
N3_x = (y4 - y2) / (2 * area);
N4_x = (y1 - y3) / (2 * area);

N1_y = (x4 - x2) / (2 * area);
N2_y = (x1 - x3) / (2 * area);
N3_y = (x2 - x4) / (2 * area);
N4_y = (x3 - x1) / (2 * area);

% Define the element stiffness matrix (4x4 matrix)
Ke = (1/(4 * area)) * [
    N1_x * N1_x + N1_y * N1_y, N1_x * N2_x + N1_y * N2_y, N1_x * N3_x + N1_y * N3_y, N1_x * N4_x + N1_y * N4_y;
    N2_x * N1_x + N2_y * N1_y, N2_x * N2_x + N2_y * N2_y, N2_x * N3_x + N2_y * N3_y, N2_x * N4_x + N2_y * N4_y;
    N3_x * N1_x + N3_y * N1_y, N3_x * N2_x + N3_y * N2_y, N3_x * N3_x + N3_y * N3_y, N3_x * N4_x + N3_y * N4_y;
    N4_x * N1_x + N4_y * N1_y, N4_x * N2_x + N4_y * N2_y, N4_x * N3_x + N4_y * N3_y, N4_x * N4_x + N4_y * N4_y
    ];

% Define the element load vector (4x1 vector)
x_mid = (x1 + x2 + x3 + x4) / 4;
y_mid = (y1 + y2 + y3 + y4) / 4;
Fe = area * f(x_mid, y_mid) * [1; 1; 1; 1];
end