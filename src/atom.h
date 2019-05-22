#ifndef ATOM_H
#define ATOM_H
#include <string>

using std::string;
const int max_neigh = 16;
class atom{
      protected:
              double x,y,z;
              double mass;
              int mol_no;
              int flag1;
              string name;
              int atom_state; //0 means yes, 1 means No
   	      int xbin,ybin,zbin; //required if neighbors are calculated based on distances
              double total_bo;
              int ring_defect_flag; //0 means part of only 6-member rings, < 0 not initialized, > 0 means member of defective ring
      public:
              int neigh_indexes[max_neigh];
              int num_neighbors;
             atom(){
                    x = 0.0;
                    y = 0.0;
                    z = 0.0;
                    mass = 0.0;
                    mol_no = 0;
                    name = "No";
  		    xbin = 0;
	            ybin = 0;
		    zbin = 0;
                    num_neighbors = 0;
                    ring_defect_flag = -1;
                    }
         atom(double xx, double yy, double zz, string at_name);
         double return_xcord(){ return x;}
         double return_ycord(){return y;}
         double return_zcord() {return z;}
         string return_atomname() { return name;}
         int return_atom_state() { return atom_state;}
         int return_mol_no(){ return mol_no;}
	 int return_xbin(){ return xbin;}
	 int return_ybin(){ return ybin;}
	 int return_zbin(){ return zbin;}
	 int return_ring_defect_flag(){ return ring_defect_flag;}
         double return_total_bo(){return total_bo;}
         
         void set_coordinates(double xx, double yy, double zz, string at_name, int n1);
         void set_coordinates(double xx, double yy, double zz, string name1);
         void set_atom_characteristics(double xx, double yy, double zz, string at_name, int n1);
         void set_atom_state(int n1);
         void set_mol_no(int i1);
	 void set_atom_bins (int xbin1, int ybin1, int zbin1);
         void set_num_neighbors(int num);
         void set_neigh_index(int at, int index);
         void set_total_bo(double value);
         void set_ring_defect_flag(int num);
};
      
#endif
