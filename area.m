%{
======================================================================
    Function that calculates jellyfish area.
======================================================================

INPUT:
        jelly (graph): Graph representation of jellyfish.

OUTPUT:
        jelly_area (float): Area of jellyfish in millimeter squared.
        edge_idx (array):   Array of edge notes.    
%}


function [jelly_area, edge_idx] = area(jelly)
    [jelly1, idx] = sortrows(jelly.Nodes, {'edges'});
    for i = 1:numnodes(jelly)
        node1 = idx(i);
        if jelly.Nodes.edges(node1) > 0
            s = i;
            break
        end
    end
    edge_idx = cat(1,idx(s:numnodes(jelly)), idx(s));
    
    jelly_area = 0;

    % Implement shoelace method to calculate jellyfish area
    for i = 1:length(edge_idx)-1
        jelly_area = jelly_area + (jelly.Nodes.x_coord(edge_idx(i))*jelly.Nodes.y_coord(edge_idx(i+1))) - (jelly.Nodes.y_coord(edge_idx(i))*jelly.Nodes.x_coord(edge_idx(i+1)));
    end
    
    jelly_area = abs(jelly_area);
end
    
            
        