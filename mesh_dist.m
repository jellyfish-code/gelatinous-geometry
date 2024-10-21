%{
======================================================================
    Function that (?)
======================================================================

INPUT:
        jelly (?): 
        row_start (?):
        row_end (?):

OUTPUT:
        dist (?):
%}

function dist = mesh_dist(jelly, row_start, row_end)
    mesh_size = size(jelly);
    dist = zeros(mesh_size(1), mesh_size(2), 3); %mm

    for i = 1:mesh_size(1)
        for j = row_start(i):row_end(i)
            if j < row_end(i)
                x1 = jelly(i, j+1, 1) - jelly(i, j, 1);
                y1 = jelly(i, j+1, 2) - jelly(i, j, 2);
                dist(i, j, 1) = (x1^2 + y1^2)^(1/2); %horiz
            end
            if i < mesh_size(1) && j+1 >= row_start(i+1) && j+1 <= row_end(i+1)
                x2 = jelly(i+1, j+1, 1) - jelly(i, j, 1);
                y2 = jelly(i+1, j+1, 2) - jelly(i, j, 2);
                dist(i,j,2) = (x2^2 + y2^2)^(1/2); %up right
            end
            if i < mesh_size(1) && j >= row_start(i+1) && j <= row_end(i+1)
                x3 = jelly(i+1, j, 1) - jelly(i, j, 1);
                y3 = jelly(i+1, j, 2) - jelly(i, j, 2);
                dist(i,j,3) = (x3^2 + y3^2)^(1/2); %up left
            end
        end
    end