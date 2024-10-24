%{
======================================================================
    Function that (?)
======================================================================

INPUT:
        jelly (graph):          Graph containing information on Nodes and Edges of jellyfish.
        area_0 (float):     
        j_area (float):
        bulk_modulus (float):   Bulk modulus of jellyfish. 
        edges (float):

OUTPUT:
        pressure (float):       Pressure inside jellyfish (in Pascals).
%}

function pressure = find_f_pressure_initial(jelly, area_0, j_area, bulk_modulus, edges)
    
    d_area = (j_area - area_0)/area_0; % dimensionless
    mesh_size = size(jelly);
    pressure = zeros(mesh_size(1), mesh_size(2), 2);
    
    d_pressure = bulk_modulus*-1*d_area; % In Pascals.
    
    %Pa = N/m^2
    
    %Force(N) = pressure(Pa)*area(m^2)
    %So find the area and the direction of the force
    %direction of force is perpendicular to the cross_area, outward if
    %d_area is negative
    
    %edges must go clockwise around the jelly
    for i = 1:length(edges)-1
        %define the side that pressure is acting on
        dx = jelly(edges(1,i+1), edges(2,i+1),1) - jelly(edges(1,i),edges(2,i),1);
        dy = jelly(edges(1,i+1), edges(2,i+1),2) - jelly(edges(1,i),edges(2,i),2);
        l = (dx^2 + dy^2)^(1/2);
        pressure(edges(1,i), edges(2,i),1) = -1*dy/l*d_pressure + pressure(edges(1,i), edges(2,i),1);
        pressure(edges(1,i), edges(2,i),2) = dx/l*d_pressure + pressure(edges(1,i), edges(2,i),2);
        
        pressure(edges(1,i+1), edges(2,i+1),1) = -1*dy/l*d_pressure + pressure(edges(1,i+1), edges(2,i+1),1);
        pressure(edges(1,i+1), edges(2,i+1),2) = dx/l*d_pressure + pressure(edges(1,i+1), edges(2,i+1),2);
    end
    