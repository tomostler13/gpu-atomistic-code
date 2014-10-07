//File matrix_conv.cpp
// Author: Tom Ostler
// Created: 07 Jan 2014
// Last-modified: 07 Oct 2014 10:35:18
// The routines within this file convert a 2D matrix to
// a number of formats depending on the routine used. The
// return structures depend on the storage format

// This calculate the diagonal offsets for the DIA format
// The first arguement is the array to store the offsets
// The second is the J matrix (4 dimensions because it
// represents all elements of the interactio matrix).
void dia_offsets(Array<unsigned int>& diaoffset,Array4D<double>& JMat,unsigned int N)
{

