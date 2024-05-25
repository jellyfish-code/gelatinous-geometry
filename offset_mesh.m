function [jelly, row_start, row_end] = offset_mesh(offset)
    top_start = [1, 2, 3, 4, 5, 6];
    bottom_start = offset + [1, 1, 1, 1, 1];
    row_start = cat(2, bottom_start, top_start);
    width = [5, 6, 7, 8, 9, 10 + offset, 9, 8, 7, 6, 5];
    row_end = row_start + width;
    mesh_x = zeros(length(row_start), max(row_end));
    mesh_y = zeros(length(row_start), max(row_end));
    edge_row = [3*sqrt(7)/2, 5*sqrt(3)/2, sqrt(83)/2, sqrt(83)/2, 5*sqrt(3)/2, 3*sqrt(7)/2];
    edge2_row = [sqrt(39)/2, 2*sqrt(3), sqrt(13), 2*sqrt(3), sqrt(39)/2];
    edge_column = [sqrt(21)/2, 3, sqrt(57)/2, 3*sqrt(2), 4.5, sqrt(21), 4.5, 3*sqrt(2), sqrt(57)/2, 3, sqrt(21)/2];
    edge2_column = [sqrt(13)/2, 2.5, sqrt(10), 3.5, sqrt(13), 3.5, sqrt(10), 2.5, sqrt(13)/2];
    
    %center of mesh
    for i = 1:length(row_start)
        s = row_start(i);
        t = width(i);
        for j = s:s+t
            mesh_x(i, j) = 2.5 - 0.5*i + j;
        end
        mesh_y(i, s:s+t) = i*sqrt(3)/2;
    end
    
        %first row
    for i = 1
        s = row_start(i);
        t = width(i);
        for j = s:s+t
            mesh_y(i, j) = 6*sqrt(3)/2 - edge_row(j-s+1);
        end
    end
    
    %second row
    for i = 2
        s = row_start(i);
        t = width(i);
        for j = s+1:(s+t)-1
            mesh_y(i, j) = 6*sqrt(3)/2 - edge2_row(j-s);
        end
    end
    
    %last row 
    for i = length(row_start)
        s = row_start(i);
        t = width(i);
        for j = s:s+t
            mesh_y(i, j) = 6*sqrt(3)/2 + edge_row(j-s+1);
        end
    end
    
    %second to last row
    for i = length(row_start)-1
        s = row_start(i);
        t = width(i);
        for j = s+1:s+t-1
            mesh_y(i, j) = 6*sqrt(3)/2 + edge2_row(j-s);
        end
    end
    
    %last column
    for i = 1:length(row_start)
        s = row_start(i) + width(i);
        if i <= length(bottom_start)+1
            mesh_x(i, s) = 5.5 + edge_column(i) + offset;
        else
            mesh_x(i, s) = 5.5 + edge_column(i);
        end
        if i == length(bottom_start)+1
            mesh_x(i, s-offset) = 5.5 + edge_column(i);
        end
    end
    
    %second to last column
    for i = 2:length(row_start)-1
        s = row_start(i) + width(i);
        if i <= length(bottom_start)+1
            mesh_x(i, s-1) = 5.5 + edge2_column(i-1) + offset;
        else
            mesh_x(i, s-1) = 5.5 + edge2_column(i-1);
        end
        if i == length(bottom_start)+1
            mesh_x(i, s-offset-1) = 5.5 + edge2_column(i-1);
        end
    end
    
    %first column
    for i = 1:length(row_start)
        s = row_start(i);
        if i <= length(bottom_start)
            mesh_x(i, s) = 5.5 - edge_column(i) + offset;
        else
            mesh_x(i, s) = 5.5 - edge_column(i);
        end
        if i == length(bottom_start)+1
            mesh_x(i, s+offset) = 5.5 - edge_column(i) + offset;
        end
    end
    
    %second to last 
    for i = 2:length(row_start)-1
        s = row_start(i);
        if i <= length(bottom_start)
            mesh_x(i, s+1) = 5.5 - edge2_column(i-1) + offset;
        else
            mesh_x(i, s+1) = 5.5 - edge2_column(i-1);
        end
        if i == length(bottom_start)+1
            mesh_x(i, s+1+offset) = 5.5 - edge2_column(i-1) + offset;
        end
    end
        
    jelly = cat(3, mesh_x, mesh_y);
    
    
        
    
