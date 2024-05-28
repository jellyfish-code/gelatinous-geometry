function F_muscle = contraction3(contraction_strength, jelly)

    %Reset F_muscle
    F_muscle = zeros(numnodes(jelly), 2);

    for i = 1:numedges(jelly)
        if jelly.Edges.muscle(i) ~= 0
            strength = jelly.Edges.muscle(i);
            node1 = jelly.Edges.EndNodes(i,1);
            node2 = jelly.Edges.EndNodes(i,2);
            
            if jelly.Nodes.outmus(node1)>jelly.Nodes.outmus(node2) || jelly.Nodes.inmus(node1)>jelly.Nodes.inmus(node2)
                %already clockwise
                %dx, dy are directed toward node 1
                dx = (jelly.Nodes.x_coord(node1) - jelly.Nodes.x_coord(node2))/jelly.Edges.d_current(i);
                dy = (jelly.Nodes.y_coord(node1) - jelly.Nodes.y_coord(node2))/jelly.Edges.d_current(i);
            else
                dx = (jelly.Nodes.x_coord(node2) - jelly.Nodes.x_coord(node1))/jelly.Edges.d_current(i);
                dy = (jelly.Nodes.y_coord(node2) - jelly.Nodes.y_coord(node1))/jelly.Edges.d_current(i);
            end
            
            F_muscle(node1,1) = dy*contraction_strength*strength + F_muscle(node1, 1);
            F_muscle(node1,2) = -1*dx*contraction_strength*strength + F_muscle(node1, 2);

            F_muscle(node2,1) = dy*contraction_strength*strength + F_muscle(node2 ,1);
            F_muscle(node2,2) = -1*dx*contraction_strength*strength + F_muscle(node2, 2);
        end
    end
            
            
   
    