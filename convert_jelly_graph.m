function jelly = convert_jelly_graph(jelly_initial, row_start, row_end)
    %% Convert array coordinates to something readable by a graph
    a = size(jelly_initial);
    
    x_coord = []; y_coord = [];
    for i = 1:length(row_start)
        cur_x = jelly_initial(i, row_start(i):row_end(i), 1);
        cur_y = jelly_initial(i, row_start(i):row_end(i), 2);
        x_coord = cat(2, x_coord, cur_x);
        y_coord = cat(2, y_coord, cur_y);
    end

    %% Convert that matrix to a graph
    num = 0;

    %Count the number of nodes
    for i = 1:a(1)
        num = num + (row_end(i) - row_start(i) + 1);
    end

    %Initialize a matrix to define edges
    mat_A = zeros(num, num);
    n_names = zeros(1,2); %names of nodes (row and column)
    %Define the edges (connections) in the matrix
    num = 0;
    for i = 1:a(1)
        for j = row_start(i):row_end(i)
            n_names = cat(1, n_names, [i,j]);
            if j < row_end(i)
                nodeid = num+(j-row_start(i)+1);
                mat_A(nodeid, nodeid+1) = 1;
                mat_A(nodeid+1, nodeid) = 1;
            end
            if i < a(1) && j+1 >= row_start(i+1) && j+1 <= row_end(i+1)
                %find ID of node up and right
                nodeid = num+(row_end(i)-row_start(i)+1)+(j+1 - row_start(i+1)+1);
                mat_A(num+(j-row_start(i)+1), nodeid) = 1;
                mat_A(nodeid ,num+(j-row_start(i))+1) = 1;    
            end
            if i < a(1) && j >= row_start(i+1) && j <= row_end(i+1)
                %find ID of node up and left
                nodeid = num+(row_end(i)-row_start(i)+1)+(j - row_start(i+1)+1);
                mat_A(num+(j-row_start(i)+1), nodeid) = 1;
                mat_A(nodeid ,num+(j-row_start(i))+1) = 1;    
            end
        end
        num = num + (row_end(i) - row_start(i) + 1);
    end

    %Actually convert to graph
    jelly = graph(mat_A);

    %Add variables to jelly graph
    jelly.Nodes.Node_Names = n_names(2:numnodes(jelly)+1, :);  
    jelly.Nodes.x_coord = x_coord.';
    jelly.Nodes.y_coord = y_coord.';

    jelly.Edges.d_current = zeros(numedges(jelly), 1);
    %%Find distance of each edge
    for i = 1:numedges(jelly)
        node1 = jelly.Edges.EndNodes(i,1);
        node2 = jelly.Edges.EndNodes(i,2);
        dx = jelly.Nodes.x_coord(node1) - jelly.Nodes.x_coord(node2);
        dy = jelly.Nodes.y_coord(node1) - jelly.Nodes.y_coord(node2);
        dist_current = (dx^2 + dy^2)^(1/2);
        jelly.Edges.d_current(i) = dist_current;
    end