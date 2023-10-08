function indices = elementConnectivity(i, j, Nx, Ny)
    % Compute the global indices of the nodes of element (i, j)
    
    % Assuming that nodes are numbered from left to right, bottom to top
    % (i = 1, j = 1) corresponds to the bottom-left node
    % (i = Nx, j = Ny) corresponds to the top-right node
    
    % Calculate global indices
    nodes_per_row = Nx + 1;
    bottom_left = (j - 1) * nodes_per_row + i;
    bottom_right = bottom_left + 1;
    top_left = bottom_left + nodes_per_row;
    top_right = top_left + 1;
    
    % Store the global indices in an array
    indices = [bottom_left, bottom_right, top_left, top_right];
    end