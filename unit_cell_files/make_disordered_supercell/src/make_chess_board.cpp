// File: make_chess_board.cpp
// Author:Tom Ostler
// Created: 22 Nov 2014
// Last-modified: 22 Nov 2014 12:43:55

//The purpose of this section of code is to create a unit cell
//file for use with the main program. The specific type of unit
//cell that this code creates is a disordered system of two
//species (the initial idea was to create Gd/Fe). Furthermore,
//the chess board idea comes from creating areas with different
//concentrations
//    .     |    .     |    .     |    .     |    .     |
//    .     |    .     |    .     |    .     |    .     |
//    .     |    .     |    .     |    .     |    .     |
//--------------------------------------------------------------
//          |          |          |          |          | ....
//  20% Gd  |  30% Gd  |  20% Gd  |  30% Gd  |  20 %Gd  | ....
//          |          |          |          |          | ....
//--------------------------------------------------------------
//          |          |          |          |          | ....
//  30% Gd  |  20% Gd  |  30% Gd  |  20% Gd  |  30 %Gd  | ....
//          |          |          |          |          | ....
//--------------------------------------------------------------
//          |          |          |          |          | ....
//  20% Gd  |  30% Gd  |  20% Gd  |  30% Gd  |  20 %Gd  | ....
//          |          |          |          |          | ....
//--------------------------------------------------------------
//    .     |    .     |    .     |    .     |    .     |
//    .     |    .     |    .     |    .     |    .     |
//    .     |    .     |    .     |    .     |    .     |
#include <cmath>
#include <iostream>
#include <libconfig.h++>
#include "../../../inc/error.h"
#include "../../../inc/array.h"
#include "../../../inc/array2d.h"
#include "../../../inc/array3d.h"
#include "../../../inc/array4d.h"
int main(int argc,char *argv[])
{

    return(EXIT_SUCCESS);
}
