%{
======================================================================
    Function that finds edges for the butterfly mesh.
======================================================================

INPUT:
        row_start (array):  For each row of the jellyfish, position(or coordinate?) of where row starts.  
        row_end (array):    For each row of the jellyfish, position(?) of where row ends.

OUTPUT:
        edges (array):      Describes the edge nodes in clockwise order. 
%}

function edges = find_edges_butterfly(row_start, row_end)
    
    edges = cat(2, row_start, flip(row_end,2), row_start(1));
    x = zeros(1, length(edges));
    edges = cat(1, edges, x);
    insert = zeros(2,1);

    i = 1;
    a=0;
    while i <= length(row_start)
        edges(2,i+a) = i;
        if abs(edges(1,i+1+a) - edges(1,i+a)) <= 1
            i = i+1;
        else
            insert(1,1) = edges(1,i+a) + (edges(1,i+1+a) - edges(1,i+a))/abs(edges(1,i+1+a) - edges(1,i+a));
            insert(2,1) = i;
            edges = cat(2, edges(:,1:i+a), insert, edges(:,i+a+1:length(edges)));
            a = a+1;
        end
    end
    countfrom = length(row_start);
    while i <= 2*length(row_start)
        edges(2,i+a) = countfrom;
        if abs(edges(1,i+1+a) - edges(1,i+a)) <= 1
            i = i+1;
            countfrom = countfrom - 1;
        else
            insert(1,1) = edges(1,i+a) + (edges(1,i+1+a) - edges(1,i+a))/abs(edges(1,i+1+a) - edges(1,i+a));
            insert(2,1) = countfrom;
            edges = cat(2, edges(:,1:i+a), insert, edges(:,i+a+1:length(edges)));
            a = a+1;
        end
    end
    
    edges(2,length(edges)) = 1;
    %made a dumb mistake, the row number is the first row and column number
    %is the second row
    edge = edges;
    edge(1,:) = edges(2,:);
    edge(2,:) = edges(1,:);
    edges = edge;