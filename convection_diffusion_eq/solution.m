[tri, x, y] = generateTriangularMesh(60);
trimesh(tri, x, y); axis('equal');

function [tri, x, y] = generateTriangularMesh(N)
    delta = 6 / N;
    counter = 1;

    % Generate initial grid
    for i = 1:60
        for j = 1:60
            x0(counter) = -2 + i * delta;
            y0(counter) = -3 + j * delta;
            counter = counter + 1;
        end
    end

    counter = 1;

    % Extract in-domain nodes
    for i = 1:3600
        if norm([x0(i)-1 y0(i)]) < 3-delta/2 && norm([x0(i) y0(i)]) > 1+delta/2
            x(counter) = x0(i);
            y(counter) = y0(i);
            counter = counter + 1;
        end
    end

    counter = counter - 1;

    % Populate boundary nodes
    bnodes_out = 90;
    for i = 1:bnodes_out
        h = 1 + 3 * cos(2 * pi / bnodes_out * i);
        v = 3 * sin(2 * pi / bnodes_out * i);
        x(i+counter) = h; y(i+counter) = v;
    end
    
    bnodes_in = 60;
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