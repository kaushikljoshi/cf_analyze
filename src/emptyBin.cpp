#include<string>
#include "emptyBin.h"

emptyBin::emptyBin(int xx, int yy, int zz){
                    xbin =xx;
                    ybin =yy;
                    zbin =zz;
                    }
                    
void emptyBin::set_bin_indexes(int xx, int yy, int zz){
                     xbin = xx;
                     ybin = yy;
                     zbin = zz;
                     }
 
void emptyBin::set_bin_coordinates(double xx, double yy, double zz){
                     xc = xx;
                     yc = yy;
                     zc = zz;
                     }
                     
void emptyBin::set_bin_state(int s1){
                     bin_state = s1;
                     }

void emptyBin::set_void_no(int i1){
                     void_no = i1;
                     }

void emptyBin::set_bin_vol(double temp_vol){
                     vol = temp_vol;
                     }
					 
void emptyBin::set_num_neighbors(int num){
                     num_neighbors = num;
                     }

void emptyBin::set_neigh_index(int count, int ind){
                     neigh_indexes[count] = ind;
                    }
