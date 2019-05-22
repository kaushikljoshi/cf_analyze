#ifndef emptyBin_H
#define emptyBin_H
#include <string>

using std::string;
const int bin_max_neigh = 6;
class emptyBin{
      protected:
              int xbin,ybin,zbin;
              double xc,yc,zc;
              int bin_state; //0 means yes, 1 means No
              int void_no;
              double vol;
   	      //int xbin,ybin,zbin; //required if neighbors are calculated based on distances
      public:
              int neigh_indexes[bin_max_neigh];
              int non_neigh_indexes[bin_max_neigh];
              int num_neighbors;
              int num_non_neighbors;
             emptyBin(){
                    xbin = 0;
                    ybin = 0;
                    zbin = 0;
                    xc = yc = zc = 0.0;
                    void_no = 0;
                    num_neighbors = 0;
                    vol = 0.0;
             }
         emptyBin(int xx, int yy, int zz);
         int return_xbin(){ return xbin;}
         int return_ybin(){return ybin;}
         int return_zbin() {return zbin;}
         double return_xc(){ return xc;}
         double return_yc(){ return yc;}
         double return_zc(){ return zc;}
         int return_bin_state() { return bin_state;}
         int return_void_no(){ return void_no;}
         double return_bin_vol(){ return vol;}
         
         void set_bin_indexes(int xx, int yy, int zz);
         void set_bin_coordinates(double xx, double yy, double zz);
         void set_bin_state(int n1);
         void set_void_no(int i1);
         void set_num_neighbors(int num);
         void set_neigh_index(int at, int index);
         void set_bin_vol(double temp_vol);
};
      
#endif
