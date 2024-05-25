function jelly = SLM_elastic(jelly, A0, A1)

    %% Update the strain on each of the springs. 
    jelly.Edges.strain0 = (jelly.Edges.d_current - jelly.Edges.d_rel0)./jelly.Edges.d_rel0;
    jelly.Edges.strain1 = (jelly.Edges.d_current - jelly.Edges.d_rel1)./jelly.Edges.d_rel1;

    %% Find the elastic force on each node 
    %F_elastic(N) = scale(m^2)*strain(unitless)*stiffness(Pa)
    jelly.Nodes.F_elastic = zeros(numnodes(jelly), 2);
    for i = 1:numedges(jelly)
        %For each edge, find the nodes that the edge acts on, and the
        %direction of the force dx, dy
        [node1, node2] = findedge(jelly, i);
        dx = jelly.Nodes.x_coord(node1) - jelly.Nodes.x_coord(node2); 
        dy = jelly.Nodes.y_coord(node1) - jelly.Nodes.y_coord(node2);

        %positive dx, dy means toward node 1. When strain is positive,
        %that means edges is under tension. Direction of force should
        %therefore be toward the middle of the edge. For node 1 that is
        %negative, for node 2 that is positive.
        jelly.Nodes.F_elastic(node1, 1) = jelly.Nodes.F_elastic(node1, 1) - (jelly.Edges.strain0(i)*A0 + jelly.Edges.strain1(i)*A1)*dx/jelly.Edges.d_current(i);
        jelly.Nodes.F_elastic(node1, 2) = jelly.Nodes.F_elastic(node1, 2) - (jelly.Edges.strain0(i)*A0 + jelly.Edges.strain1(i)*A1)*dy/jelly.Edges.d_current(i);

        jelly.Nodes.F_elastic(node2, 1) = jelly.Nodes.F_elastic(node2, 1) + (jelly.Edges.strain0(i)*A0 + jelly.Edges.strain1(i)*A1)*dx/jelly.Edges.d_current(i);
        jelly.Nodes.F_elastic(node2, 2) = jelly.Nodes.F_elastic(node2, 2) + (jelly.Edges.strain0(i)*A0 + jelly.Edges.strain1(i)*A1)*dy/jelly.Edges.d_current(i);
    end