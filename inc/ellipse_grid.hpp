#include <string>
#ifndef _ELLIPSE_GRID_H_
#define _ELLIPSE_GRID_H_
double *ellipse_grid ( int n, double r[2], double c[2], int ng );
int ellipse_grid_count ( int n, double r[2], double c[2] );
int i4_ceiling ( double x );
void r82vec_print_part ( int n, double a[], int max_print, std::string title );
void r8mat_write ( std::string output_filename, int m, int n, double table[] );
void timestamp ( );
#endif/*_ELLIPSE_GRID_H_*/
