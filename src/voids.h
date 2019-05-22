#ifndef VOIDS_H
#define VOIDS_H
#include <string>
#include <vector>
#include "functions.h"
using namespace std;
const int max_void_types = 1000;

class voids {
      protected:
              int no_bins;
              int void_number;
              double void_vol;
      public:
             vector < int > vbin_list;
             vector < int > edge_list;
             vector < int > edge_atom_list;
             vector < double > face_xc;
             vector < double > face_yc;
             vector < double > face_zc;
             double max_len,index1,index2;
             voids(){
                        no_bins = 0;
                        void_number = 0;
                        void_vol = 0.0;
                        }
             void set_void_number( int n1){ void_number = n1;}
             
             void increase_bin_nos(){
                  no_bins = no_bins + 1;
                  }

             int return_no_bins(){ return no_bins;}
             
             void set_void_vol(double m1){void_vol = m1;}
 
             void set_void_no_bins(int m1){no_bins = m1;}
             
             double return_void_vol(){return void_vol;}
             
             int return_void_number(){return void_number;}
      
      };

#endif
