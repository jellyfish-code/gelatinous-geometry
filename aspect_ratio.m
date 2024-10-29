%{
======================================================================
    Function to calculate aspect ratio of jellyfish.
======================================================================
    INPUT:
        jelly (graph):          Graph containing information on Nodes and Edges of jellyfish.

    OUTPUT:
        a_r (double):            Ratio of major axis to minor axis.                      
        major_axis (double):     Longest distance found between opposite boundary nodes.
        minor_axis (double):     Distance between boundary nodes perpendicular to major axis.
%}
function [a_r, major_axis, minor_axis] = aspect_ratio(jelly)

[trash, edge_idx] = area(jelly);

x_s = jelly.Nodes.x_coord(edge_idx);
y_s = jelly.Nodes.y_coord(edge_idx);

a = floor(length(edge_idx)/2);

axes = [];
ang = [];

%% We want to find the two nodes that are the furthest away from each other who are on the edge
%%Making the assumption that there is some symmetry is probably good

for i = 1:a
    dx = x_s(i) - x_s(i+a);
    dy = y_s(i) - y_s(i+a);
    
    axes = cat(1, axes, (dx^2 + dy^2)^(1/2));
    ang = cat(1, ang, atan(dy/dx));
end

%%Find the major axis, meaning the longest
major_axis = max(axes);
major_ang = ang(axes == major_axis);

%%Find the minor axis, meaning perpendicular to major axis
find_minor = abs(abs(ang - major_ang) - pi/2);
minor_axis = axes(find_minor == min(find_minor));
a_r = major_axis/minor_axis;