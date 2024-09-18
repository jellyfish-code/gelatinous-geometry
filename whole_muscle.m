%{
======================================================================
    Function that initialises a set of nodes for the muscle layers.
======================================================================

INPUT:
        None.

OUTPUT:
        muscle_outer (matrix):             A matrix storing index of jelly array corresponding to outer muscles.
        muscle_inner (matrix):             A matrix storing index of jelly array corresponding to inner muscles.
%}


function [muscle_outer, muscle_inner] = whole_muscle()
    muscle_top2 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 11, 11, 11, 11];
    muscle_top1 = [6, 7, 8, 9, 10, 11, 11, 11, 11, 11, 11, 10, 9, 8, 7, 6];
    muscle_bottom2 = [1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
    muscle_bottom1 = [6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6];
    muscle1 = cat(1, muscle_top1, muscle_top2);
    muscle2 = cat(1, flip(muscle_bottom1), flip(muscle_bottom2));
    
    muscle2_top2 = [2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 10];
    muscle2_top1 = [6, 7, 8, 9, 10, 10, 10, 10, 10, 9, 8, 7, 6];
    muscle2_bottom2 = [2, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    muscle2_bottom1 = [6, 5, 4, 3, 2, 2, 2, 2, 2, 3, 4, 5, 6];
    muscle3 = cat(1, muscle2_top1, muscle2_top2);
    muscle4 = cat(1, flip(muscle2_bottom1), flip(muscle2_bottom2));
    
    muscle_outer = cat(3, muscle1, muscle2);
    muscle_inner = cat(3, muscle3, muscle4);
    