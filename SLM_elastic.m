%{
======================================================================
    Function that computes the strain between edges and elastic stress acting on each node.
======================================================================

INPUT:
        jelly (graph):         Graph containing information on Nodes and Edges of jellyfish.
        elast0 (float):        Elasticity of spring (in Pascals).
        elast1 (float):        Elasticity of spring (in Pascals).

OUTPUT:
        jelly (graph):         Updated graph containing information on Nodes and Edges of jellyfish.
%}

function jelly = SLM_elastic(jelly, elast0, elast1)

    %% Update the strain on each of the springs. 
    jelly.Edges.strain0 = (jelly.Edges.d_current - jelly.Edges.d_rel0)./jelly.Edges.d_rel0;
    jelly.Edges.strain1 = (jelly.Edges.d_current - jelly.Edges.d_rel1)./jelly.Edges.d_rel1;

    %% Find the elastic stress on each node 
    %stress_elastic(N) = scale(m^2)*strain(unitless)*stiffness(Pa)
    jelly.Nodes.stress_elastic = zeros(numnodes(jelly), 2);
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
        jelly.Nodes.stress_elastic(node1, 1) = jelly.Nodes.stress_elastic(node1, 1) - (jelly.Edges.strain0(i)*elast0 + jelly.Edges.strain1(i)*elast1)*dx/jelly.Edges.d_current(i);
        jelly.Nodes.stress_elastic(node1, 2) = jelly.Nodes.stress_elastic(node1, 2) - (jelly.Edges.strain0(i)*elast0 + jelly.Edges.strain1(i)*elast1)*dy/jelly.Edges.d_current(i);

        jelly.Nodes.stress_elastic(node2, 1) = jelly.Nodes.stress_elastic(node2, 1) + (jelly.Edges.strain0(i)*elast0 + jelly.Edges.strain1(i)*elast1)*dx/jelly.Edges.d_current(i);
        jelly.Nodes.stress_elastic(node2, 2) = jelly.Nodes.stress_elastic(node2, 2) + (jelly.Edges.strain0(i)*elast0 + jelly.Edges.strain1(i)*elast1)*dy/jelly.Edges.d_current(i);
    end