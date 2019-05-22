#include<string>
#include "atom.h"

atom::atom(double xx, double yy, double zz, string at_name){
                    x =xx;
                    y =yy;
                    z =zz;
                    name = at_name;
                    mass = 0.0;
                    }
                    
/*void atom::set_coordinates(double xx, double yy, double zz, string at_name, int n1){
                     x = xx;
                     y = yy;
                     z = zz;
                     name = at_name;
                     mol_no = n1;
                     }*/
 
void atom::set_coordinates(double xx, double yy, double zz, string name1){
                     x = xx;
                     y = yy;
                     z = zz;
                     name = name1;
                     }
 
void atom::set_atom_characteristics(double xx, double yy, double zz, string at_name, int n1){
                     x = xx;
                     y = yy;
                     z = zz;
                     n1 = mol_no;
                     name = at_name;
                     }
                     
void atom::set_atom_state(int s1){
                     atom_state = s1;
                     }

void atom::set_mol_no(int i1){
                     mol_no = i1;
                     }
					 
void atom::set_ring_defect_flag(int i1){
                      ring_defect_flag = i1;
                     }

void atom::set_atom_bins(int xbin1, int ybin1, int zbin1){
                     xbin = xbin1;
                     ybin = ybin1;
                     zbin = zbin1;
                     }

void atom::set_num_neighbors(int num){
                     num_neighbors = num;
                     }

void atom::set_neigh_index(int count, int ind){
                     neigh_indexes[count] = ind;
                    }
void atom::set_total_bo(double value){
                     total_bo = value;
}
