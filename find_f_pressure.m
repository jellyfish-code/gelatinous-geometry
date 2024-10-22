%{
======================================================================
    Function that (?)
======================================================================

INPUT:
        muscle_outer (?):

OUTPUT:
        edges (?):
%}

function F_pressure = find_f_pressure(jelly, area_relaxed, bulk_modulus)

    [area_current, edge_idx] = area(jelly);
    
    d_area = (area_current - area_relaxed)/area_relaxed;    % Dimensionless
    
    F_pressure = zeros(numnodes(jelly), 2); % Initialise pressure variable

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
        
        F_pressure(edge_idx(i), 1) = -1*dy*d_pressure + F_pressure(edge_idx(i), 1);
        F_pressure(edge_idx(i), 2) = dx*d_pressure + F_pressure(edge_idx(i), 2);
        
        F_pressure(edge_idx(i+1),1) = -1*dy*d_pressure + F_pressure(edge_idx(i+1), 1);
        F_pressure(edge_idx(i+1),2) = dx*d_pressure + F_pressure(edge_idx(i+1), 2);
    end
    