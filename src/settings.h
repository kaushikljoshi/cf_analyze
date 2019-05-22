#ifndef settings_H
#define settings_H
#include <string>

using std::string;
class settings{
      public:
             int table_flag;
             int ring_flag;
             int void_flag;
             int periodic_flag;
             int empty_bin_cut_off;
             double grid_size;
             double bo_cut_off;
             double distance_cut_off;
             string struct_name;

         settings(){
             table_flag = ring_flag = void_flag = periodic_flag = 0;
             empty_bin_cut_off = 1;
             grid_size = 4.0;
             bo_cut_off = 0.3;
             distance_cut_off = 1.7;
             struct_name = "structure.xyz";
         }
             
         int return_ring_flag(){ return ring_flag;}
         int return_void_flag(){ return void_flag;}
         int return_periodic_flag(){ return periodic_flag;}
         int return_empty_bin_cut_off(){ return empty_bin_cut_off;}
         double return_grid_size(){ return grid_size;}
         double return_distance_cut_off(){ return distance_cut_off;}
         
};
      
#endif
