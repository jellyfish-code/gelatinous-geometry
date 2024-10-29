%{
======================================================================
    Function that computes pressure in jellyfish.
======================================================================

INPUT:
        jelly (graph):          Graph containing information on Nodes and Edges of jellyfish.
        area_relaxed (float):   Relaxed area used to calculate internal pressure (in mm^2). Calculated by
                                multiplying initial graft area with percentage factor (greater than 1) to account for release
                                of internal elastic stress after dissection.
        bulk_modulus (float):   Bulk modulus of jellyfish (in Pascals). 

OUTPUT:
        pressure (float):   Pressure inside jellyfish (in Pascals).
%}

function pressure = find_f_pressure(jelly, area_relaxed, bulk_modulus)

    [area_current, edge_idx] = area(jelly);
    
    d_area = (area_current - area_relaxed)/area_relaxed;    % Dimensionless
    
    pressure = zeros(numnodes(jelly), 2); % Initialise pressure variable

    %d_pressure = bulk_mod*dV/V
    %Pressure(Pa) = bulk_mod(Pa) * -d_area(unitless)
    %Assumes thickness doesn't change
    d_pressure = bulk_modulus*-1*d_area;
    
    %Pa = N/m^2
    
    %Force(N) = pressure(Pa)*area(m^2)
    %So find the area and the direction of the force
    %direction of force is perpendicular to the cross_area, outward if
    %d_area is negative
    
    %edges must go clockwise around the jelly
    for i = 1:length(edge_idx)-1
        
        %define the side that pressure is acting on              
%         edge_cur = findedge(jelly, edge_idx(i), edge_idx(i+1));
%         l = jelly.Edges.d_current(edge_cur);
        
        dx = jelly.Nodes.x_coord(edge_idx(i+1)) - jelly.Nodes.x_coord(edge_idx(i));
        dy = jelly.Nodes.y_coord(edge_idx(i+1)) - jelly.Nodes.y_coord(edge_idx(i));
        
        pressure(edge_idx(i), 1) = -1*dy*d_pressure + pressure(edge_idx(i), 1);
        pressure(edge_idx(i), 2) = dx*d_pressure + pressure(edge_idx(i), 2);
        
        pressure(edge_idx(i+1),1) = -1*dy*d_pressure + pressure(edge_idx(i+1), 1);
        pressure(edge_idx(i+1),2) = dx*d_pressure + pressure(edge_idx(i+1), 2);
    end
    