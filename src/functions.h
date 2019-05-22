#include "atom.h"
#include "emptyBin.h"
#include "voids.h"
#include "settings.h"
#include<iostream>
#include<sstream>
#include<fstream>

using namespace std;

//const int max_neigh =8 ;

#ifndef FIND_NO_LINES_H
#define FIND_NO_LINES_H

void find_no_lines(int &nlines1, int &natoms1, int &nframes1, string filename);

#endif

#ifndef FIND_NO_CONNT_LINES_H
#define FIND_NO_CONNT_LINES_H

int find_no_connt_lines(istream &myfile);

#endif

#ifndef READ_INPUT_XMOLOUT_H
#define READ_INPUT_XMOLOUT_H

int read_input_xmolout(istream &myfile,int nframes1, int natoms1,int periodic, atom **at_list1, string *cell_parameter_lines1,
                        double **lx1, double **ly1, double **lz1, double *alpha1, double *beta1, double *gama1);                        

#endif

#ifndef FIND_CELL_PARAMETERS_H
#define FIND_CELL_PARAMETERS_H

void find_cell_parameters(string s1, int frame_no, double **lx1, double **ly1, double **lz1, double *alpha1, double *beta1, double *gama1);

#endif

#ifndef GET_SETTINGS_H
#define GET_SETTINGS_H

void get_settings(settings& settings_list1);

#endif

#ifndef FIND_DISTANCE_H
#define FIND_DISTANCE_H

double find_distance(int nframes1,int target_atoms1,int object_atoms1,
                   int *t_atoms1, int *o_atoms1,atom **at_list1,
                   double **lx1,double **ly1,double **lz1,double *alpha1,
                   double *beta1,double *gama1, int frame1);

#endif

#ifndef STRING_TO_INTEGER_H
#define STRING_TO_INTEGER_H

int string_to_integer(char *s1);

#endif

#ifndef STRING_TO_DOUBLE_H
#define STRING_TO_DOUBLE_H

double string_to_double(char *s1);

#endif

#ifndef FLOAT_TO_STRING_H
#define FLOAT_TO_STRING_H

string float_to_string(int a);

#endif

#ifndef WRITE_LAMMPS_INPUT_H
#define WRITE_LAMMPS_INPUT_H

void write_lammps_input(int natoms, int no_atom_types, int lines_per_section, atom **at_list,double **lx, double **ly, double **lz,
                    double *alpha, double *beta, double *gama, string *atom_types, string *ip_atom_types, double *mass);

#endif

#ifndef LAMMPS_OUTPUT_H
#define LAMMPS_OUTPUT_H

void lammps_output(string *ip_atom_types);

#endif

#ifndef MAKE_CONNECTION_TABLE_H
#define MAKE_CONNECTION_TABLE_H

void make_connection_table(ofstream &file1, int nframes, int natoms, atom **at_list,
                           double **lx, double **ly, double **lz, double *alpha, double *beta, double *gama,
                           int ****grid_list1, int ***grid_atom_count1, int nx1, int ny1, int nz1, settings settings_list1);
                          
#endif

#ifndef IDENTIFY_HEADER_LINES_H
#define IDENTIFY_HEADER_LINES_H

int identify_header_lines(istream &myfile);

#endif

#ifndef INPUT_CONNECTION_TABLE_H
#define INPUT_CONNECTION_TABLE_H

void input_connection_table(istream &myfile, int frmae_no, int natoms, int hlines, atom **at_list);

#endif

#ifndef READ_TABLE_H
#define READ_TABLE_H

void read_table(istream& myfile, int frame_no, int natoms1, int hlines, double cut_off,atom **at_list1);

#endif

#ifndef READ_REAX_TABLE_H
#define READ_REAX_TABLE_H

void read_reax_table(istream& myfile, int frame_no, int natoms1, double cut_off, atom **at_list1);

#endif


#ifndef BIN_ATOMS_H
#define BIN_ATOMS_H

void bin_atoms(int nframes1, int natoms1, double *grid_dr1, int grid_max_atoms1, atom **at_list1, double **lx1, double **ly1, double **lz1,int ****grid_list, int ***grid_atom_count, int nx1, int ny1, int nz1, settings settings_list1);

#endif

#ifndef RING_ANALYSIS_H
#define RING_ANALYSIS_H
void ring_analysis(int act_frame1, int natoms1, atom **at_list1, int nx1, int ny1, int nz1, double **lx1, double **ly1, double **lz1,settings settings_list1);
#endif

#ifndef IDENTIFY_VOIDS_H
#define IDENTIFY_VOIDS_H
void identify_voids(int act_frame1, int natoms1,double *grid_dr1,int no_empty_bins, emptyBin *bin_list1, int ***grid_atom_count, int ****grid_list1, atom **at_list1, double **lx1, double **ly1, double **lz1, int nx1, int ny1, int nz1);
#endif

#ifndef EQUATION_OF_PLANE_H
#define EQUATION_OF_PLANE_H
void equation_of_plane(double (&points)[3][3], double (&params)[4]);
#endif
