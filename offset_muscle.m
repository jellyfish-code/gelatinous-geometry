%{
======================================================================
    Function that (?)
======================================================================

INPUT:
        offset (?):

OUTPUT:
        muscle_outer (?): 
        muscle_inner (?):
%}

function [muscle_outer, muscle_inner] = offset_muscle(offset)
    muscle_top2 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 11, 11, 11, 11];
    muscle_top1 = [6, 7, 8, 9, 10, 11, 11, 11, 11, 11, 11, 10, 9, 8, 7, 6];
    muscle_bottom2 = offset + [1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
    muscle_bottom1 = [6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6];
    muscle1 = cat(1, muscle_top1, muscle_top2);
    muscle2 = cat(1, flip(muscle_bottom1), flip(muscle_bottom2));
    
    muscle2_top2 = [2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 10];
    muscle2_top1 = [6, 7, 8, 9, 10, 10, 10, 10, 10, 9, 8, 7, 6];
    muscle2_bottom2 = offset + [2, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    muscle2_bottom1 = [6, 5, 4, 3, 2, 2, 2, 2, 2, 3, 4, 5, 6];
    muscle3 = cat(1, muscle2_top1, muscle2_top2);
    muscle4 = cat(1, flip(muscle2_bottom1), flip(muscle2_bottom2));
    
    muscle_outer = cat(3, muscle1, muscle2);
    muscle_inner = cat(3, muscle3, muscle4);
    