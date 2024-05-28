function [jelly_initial, row_start, row_end] = butterfly_mesh()
    wing_width = 11;
    body_height = 2;
    cut_width = 5;

    body_start = [7, 6, 6, 7, 9];
    body_end = [9, 11, 11, 12, 11];

    row_start = [1, 2, 3, 4, 5, 6, 6, 7, 7, 7, 7, 7, 7];
    row_end = [11, 11, 11, 11, 11, 11, 11, 12, 13, 14, 15, 16, 17];
    h = length(row_start);

    mesh_x = zeros(length(row_start), max(row_end));
    mesh_y = zeros(length(row_start), max(row_end));
    body_x = mesh_x;
    body_y = mesh_y;

    edge_row = [3*sqrt(7)/2, 5*sqrt(3)/2, sqrt(83)/2, sqrt(83)/2, 5*sqrt(3)/2, 3*sqrt(7)/2];
    edge2_row = [sqrt(39)/2, 2*sqrt(3), sqrt(13), 2*sqrt(3), sqrt(39)/2];
    edge_column = [sqrt(21)/2, 3, sqrt(57)/2, 3*sqrt(2), 4.5, sqrt(21), 4.5, 3*sqrt(2), sqrt(57)/2, 3, sqrt(21)/2];
    edge2_column = [sqrt(13)/2, 2.5, sqrt(10), 3.5, sqrt(13), 3.5, sqrt(10), 2.5, sqrt(13)/2];

    %% Define the centers of all of the pieces
    wing1 = [6, 1]; %lower
    wing2 = [6, 1+wing_width];
    body = [6 + body_height, 1+wing_width/2];

    %% Fill in the basic mesh of the wings
    %wing1
    c_x = mean([row_start(1), row_end(1)]);
    for i = 1:6
        for j = row_start(i):row_end(i)
            mesh_x(i,j) = wing1(1) + (j - c_x) - 0.5*(i-1);
        end
        mesh_y(i,row_start(i):row_end(i)) = wing1(2)+(i-1)*sqrt(3)/2;
    end

    %wing2
    c_x = mean([row_start(h), row_end(h)]);
    for i = 1:6
        for j = row_start(h-i+1):row_end(h-i+1)
            mesh_x(h-i+1,j) = wing2(1) + (j - c_x) + 0.5*(i-1);
        end
        mesh_y(h-i+1,row_start(h-i+1):row_end(h-i+1)) = wing2(2)-(i-1)*sqrt(3)/2;
    end

    %% Edges and muscles of wings
    %wing1
    for j = row_start(6):row_end(6)
        mesh_y(6, j) = wing1(2) + edge_row(j-row_start(6)+1);
    end

    for j = row_start(5)+1:row_end(5)-1
        mesh_y(5, j) = wing1(2) + edge2_row(j-row_start(5));
    end

    for i = 1:6
        mesh_x(i,row_end(i)) = edge_column(i+5) + wing1(1);
        mesh_x(i,row_start(i)) = wing1(1) - edge_column(i+5);
    end
    for i = 1:5
        mesh_x(i,row_end(i)-1) = edge2_column(i+4) + wing1(1);
        mesh_x(i,row_start(i)+1) = wing1(1) - edge2_column(i+4);
    end

    %wing2
    for j = row_start(h-5):row_end(h-5)
        mesh_y(h-5, j) = wing2(2) - edge_row(j-row_start(h-5)+1);
    end

    for j = row_start(h-4)+1:row_end(h-4)-1
        mesh_y(h-4, j) = wing2(2) - edge2_row(j-row_start(h-4));
    end

    for i = h-5:h
        mesh_x(i,row_end(i)) = edge_column(i-(h-6)) + wing2(1);
        mesh_x(i,row_start(i)) = wing2(1) - edge_column(i-(h-6));
    end
    for i = h-4:h
        mesh_x(i,row_end(i)-1) = edge2_column(i-(h-5)) + wing2(1);
        mesh_x(i,row_start(i)+1) = wing2(1) - edge2_column(i-(h-5));
    end

    c_x = 11;
    c_y = 7;
    %body center mesh
    for i=5:5+cut_width-1
        for j = body_start(i-4):body_end(i-4)
            %for i = c_y, j=c_x, body_x = body(1)
            body_x(i,j) = body(1) + (j - c_x) - 0.5*(i-c_y);
        end
        body_y(i,body_start(i-4):body_end(i-4)) = body(2)+(i-c_y)*sqrt(3)/2;
    end

    %body edges
    for i = 6:8
        body_x(i,body_start(i-4)) = body(1) - edge_column(i-1);
        body_x(i,body_start(i-4)+1) = body(1) - edge2_column(i-2);
    end
    body_x(5,body_start(5-4)) = body(1) - edge2_column(5-2);
    body_x(9,body_start(9-4)) = body(1) - edge2_column(9-2);


    %% integrate the pieces
    for i=5:5+cut_width-1
        for j = body_start(i-4):body_end(i-4)
            if mesh_x(i,j) ~= 0 || mesh_y(i,j) ~= 0
                mesh_x(i,j) = (mesh_x(i,j) + body_x(i,j))/2;
                mesh_y(i,j) = (mesh_y(i,j) + body_y(i,j))/2;
            else
                mesh_x(i,j) = body_x(i,j);
                mesh_y(i,j) = body_y(i,j);
            end
        end
    end


    jelly_initial = cat(3, mesh_x, mesh_y);

%     for n = 1:length(row_start)      %beads
%         figure1 = scatter(jelly_initial(n,row_start(n):row_end(n), 1), jelly_initial(n,row_start(n):row_end(n), 2), 20, 'filled');
%         hold on
%     end
%     hold off
