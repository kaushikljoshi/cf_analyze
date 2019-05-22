#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include "atom.h"
#include "emptyBin.h"
#include "functions.h"
#include "settings.h"

const double fill_cut_off = 0.1; //cell needs to have atleast 10% atoms of average atomic density

void find_voids(int act_frame1, int natoms, double *grid_dr1, int ***grid_atom_count1,int ****grid_list1, atom **at_list1, double **lx1, double **ly1, double **lz1,int nx1, int ny1, int nz1, int empty_bins);

void bin_atoms(int nframes1, int natoms1, double *grid_dr1, int grid_max_atoms1, atom **at_list1,double **lx1,double **ly1,double **lz1,int ****grid_list1, int ***grid_atom_count1,int nx1, int ny1, int nz1, settings settings_list1){
	int i,j,k,l;
	double x,y,z;
        double xtemp,ytemp,ztemp;
        double xlow,ylow,zlow,xhigh,yhigh,zhigh;
        int small_vol_flag;
	int xb,yb,zb;
        int ind1,ind2;
        int empty_bins;
        double avg_atoms_per_bin;
        ofstream ofile1,ofile2;
        
       
        cout << "starting binning atoms" << endl;      
        //ofile1.open("bin_data.txt");
        //ofile2.open("bin.xyz");
	
	for (i=0;i<1;i++){
		for (j=0;j<natoms1;j++){
	                //Following if-else is necessary for situations when division coordinate is perfect intergral multiple of grid_dr	
			if (int((at_list1[i][j].return_xcord()-lx1[i][0])/grid_dr1[0]) == nx1){
                            xb = nx1 -1;
                        }else{
                            xb = int((at_list1[i][j].return_xcord()-lx1[i][0])/grid_dr1[0]);
                        }
			if (int((at_list1[i][j].return_ycord()-ly1[i][0])/grid_dr1[1]) == ny1){
                            yb = ny1 -1;
                        }else{
			    yb = int((at_list1[i][j].return_ycord()-ly1[i][0])/grid_dr1[1]);
                        }
                        if (int((at_list1[i][j].return_zcord()-lz1[i][0])/grid_dr1[2]) == nz1){
                            zb = nz1 -1 ;
                        }else{
			    zb = int((at_list1[i][j].return_zcord()-lz1[i][0])/grid_dr1[2]);
                        }
			
                        //cout << "done atom " << j << " " << xb << " " << yb << " " << zb << endl;
			//cout << j+1 << " " << at_list1[i][j].return_xcord()-lx1[i][0] << endl;
                        //cout << "done atom " << j << " " << at_list1[i][j].return_xcord() << " " << at_list1[i][j].return_ycord() << " " << at_list1[i][j].return_zcord() << endl;
			if (xb > nx1-1 || yb > ny1-1 || zb > nz1-1){
				cout << "calculated bin size of atom " << j << " exceeds grid size" << endl;
                                cout << xb << " " << yb << " " << zb << endl;
                                cout << at_list1[i][j].return_xcord() << " " << at_list1[i][j].return_ycord() << " " << at_list1[i][j].return_zcord() << endl;
				exit(EXIT_FAILURE);
			}else{
				at_list1[i][j].set_atom_bins(xb,yb,zb);
                                ind1 = grid_atom_count1[xb][yb][zb];
				grid_list1[xb][yb][zb][ind1] = j;
				grid_atom_count1[xb][yb][zb] = grid_atom_count1[xb][yb][zb] + 1;
				if (grid_atom_count1[xb][yb][zb] > grid_max_atoms1){
				cout << "Number of atoms in bin " << xb << " " << yb << " " << zb << " exceeded grid_max_atoms" << endl;
					exit(EXIT_FAILURE);
				}
                        //cout << "done atom " << j << " " << xb << " " << yb << " " << zb << endl;
			}
		}
                //write bins summary
                avg_atoms_per_bin = natoms1/(nx1*ny1*nz1);
                empty_bins = 0;
                //ofile2 << nx1*ny1*nz1 << endl;
                //ofile2 << lx1[i][0] << " " << lx1[i][1] << " " << ly1[i][0] << " " << ly1[i][1] << " " << lz1[i][0] << " " << lz1[i][1] << " 90 90 90" << endl;
                for (j=0;j<nx1;j++){
                     xtemp = (j+0.5)*grid_dr1[0] + lx1[0][0];
                    for (k=0;k<ny1;k++){
                         ytemp = (k+0.5)*grid_dr1[1] + ly1[0][0];
                         for (l=0;l<nz1;l++){
                              ztemp = (l+0.5)*grid_dr1[2]+lz1[0][0];
                              small_vol_flag = 0;
                              xlow = j*grid_dr1[0] + lx1[0][0];
                              ylow = k*grid_dr1[1] + ly1[0][0];
                              zlow = l*grid_dr1[2] + lz1[0][0];
                              xhigh = xlow + grid_dr1[0];
                              yhigh = ylow + grid_dr1[1];
                              zhigh = zlow + grid_dr1[2];
                              if (xhigh > lx1[0][1]) xhigh = lx1[0][1];
                              if (yhigh > ly1[0][1]) yhigh = ly1[0][1];
                              if (zhigh > lz1[0][1]) zhigh = lz1[0][1];
                              //if ((xhigh -xlow) < grid_dr1[0]/4.0 || (yhigh - ylow)<grid_dr1[1]/4.0 || (zhigh - zlow )<grid_dr1[2]/4.0) small_vol_flag = 1;
                             //if (grid_atom_count1[j][k][l] < avg_atoms_per_bin*fill_cut_off && small_vol_flag == 0) {
                             if (grid_atom_count1[j][k][l] <= settings_list1.empty_bin_cut_off) {
                             //if (grid_atom_count1[j][k][l] <= 1) {
                                 empty_bins = empty_bins + 1;
                                 //ofile2 << "1 " << xtemp << " " << ytemp << " " << ztemp << endl; 
                             }else{
                                 //ofile2 << "2 " << xtemp << " " << ytemp << " " << ztemp << endl; 
                             }
                             //ofile1 << j+1 << " " << k+1 << " " << l+1 << " " << grid_atom_count1[j][k][l] << endl;
                             //for (int m=0;m<grid_atom_count1[j][k][l];m++){
                             //    ofile1 << grid_list1[j][k][l][m]+1 << " " ;
                             //}
                             //ofile1 << endl;
                         }
                    }
                }
	}
        //ofile1.close();
        //ofile2.close();
        cout << "Number of empty bins is " << empty_bins << endl;
        cout << "done binning atoms" << endl;

       if (settings_list1.void_flag == 1) find_voids(nframes1, natoms1, grid_dr1, grid_atom_count1, grid_list1, at_list1, lx1, ly1, lz1, nx1, ny1, nz1, empty_bins);
}
                           
void find_voids(int act_frame1, int natoms1, double *grid_dr1, int ***grid_atom_count1,int ****grid_list1, atom **at_list1,double **lx1, double **ly1, double **lz1,int nx1, int ny1, int nz1, int empty_bin_count){
   int i,j,k;
   int p,q,r;
   int temp_count,flag,temp_neigh;
   int ibin,jbin,kbin;
   int ibin1,jbin1,kbin1;
   int ilow,ihigh,jlow,jhigh,klow,khigh;
   int small_vol_flag;
   double xlen,ylen,zlen;
   double xlow,xhigh,ylow,yhigh,zlow,zhigh;
   double xc,yc,zc;
   emptyBin ebins[empty_bin_count];
   double cell_vol = (lx1[0][2]*ly1[0][2]*lz1[0][2])/(nx1*ny1*nz1);

   cout << "printing grid_dr1 from bin_atoms" << endl;
   cout << grid_dr1[0] << " " << grid_dr1[1] << " " << grid_dr1[2] << endl;
   temp_count = 0;
   for (i=0;i<nx1;i++){
       for (j=0;j<ny1;j++){
           for (k=0;k<nz1;k++){
                if (grid_atom_count1[i][j][k] <= 1){
                    xlow = i*grid_dr1[0] + lx1[0][0];
                    ylow = j*grid_dr1[1] + ly1[0][0];
                    zlow = k*grid_dr1[2] + lz1[0][0];
                    xhigh = xlow + grid_dr1[0];
                    yhigh = ylow + grid_dr1[1];
                    zhigh = zlow + grid_dr1[2];
                    if (xhigh > lx1[0][1]) xhigh = lx1[0][1];
                    if (yhigh > ly1[0][1]) yhigh = ly1[0][1];
                    if (zhigh > lz1[0][1]) zhigh = lz1[0][1];
                    small_vol_flag = 0;
                    //if ((xhigh -xlow) < grid_dr1[0]/4.0 || (yhigh - ylow)<grid_dr1[1]/4.0 || (zhigh - zlow )<grid_dr1[2]/4.0) small_vol_flag = 1;
                    if (small_vol_flag == 0){
                        xc = (xlow+xhigh)/2.0;
                        yc = (ylow+yhigh)/2.0;
                        zc = (zlow+zhigh)/2.0;

                        cell_vol = (xhigh-xlow)*(yhigh-ylow)*(zhigh-zlow);
                        /*if (cell_vol != grid_dr1[0]*grid_dr1[1]*grid_dr1[2]){
                            cout << cell_vol << " " << grid_dr1[0]*grid_dr1[1]*grid_dr1[2] << endl;
                            cout << "different cell volume " << endl;
                            cout << xlow << " " << xhigh << " " << ylow << " " << yhigh << " " << zlow << " " << zhigh << endl;
                            cout << xlow - xhigh << " " << ylow - yhigh << " " << zlow - zhigh << endl;
                            exit(1);
                        }*/
                    
                        ebins[temp_count].set_bin_indexes(i,j,k);
                        ebins[temp_count].set_bin_vol(cell_vol);
                        ebins[temp_count].num_neighbors = 0;
                        ebins[temp_count].set_bin_coordinates(xc,yc,zc);
                        temp_count = temp_count + 1;
                    //find empty neighboring bins
                    }
                }
           }
       }
   }
   cout << "Temp count is " << temp_count << endl;
   //Now assign only those bins as neighbors which are empty
   for (i=0;i<empty_bin_count;i++){
        ibin = ebins[i].return_xbin();
        jbin = ebins[i].return_ybin();
        kbin = ebins[i].return_zbin();

        //cout << i+1 << " " << ibin << " " << jbin << " " << kbin << endl;
 
        ilow  = ibin - 1;
        ihigh = ibin + 1;

        jlow  = jbin - 1;
        jhigh = jbin + 1;

        klow  = kbin - 1;
        khigh = kbin + 1;

        if (ilow < 0) ilow = nx1-1;
        if (ihigh > nx1-1) ihigh = 0;

        if (jlow < 0) jlow = ny1-1;
        if (jhigh > ny1-1) jhigh = 0;

        if (klow < 0) klow = nz1-1;
        if (khigh > nz1-1) khigh = 0;
       
        for (j=i+1;j<empty_bin_count;j++){
             ibin1 = ebins[j].return_xbin();
             jbin1 = ebins[j].return_ybin();
             kbin1 = ebins[j].return_zbin();
             flag =0;
            
             if (ibin1 == ilow || ibin1 == ihigh) {
                 if (jbin1 == jbin && kbin1 == kbin) flag = flag + 1;
             }else if (jbin1 == jlow || jbin1 == jhigh) {
                 if (ibin1 == ibin && kbin1 == kbin) flag = flag + 1;
             }else if (kbin1 == klow || kbin1 == khigh) {
                 if (ibin1 == ibin && jbin1 == jbin) flag = flag + 1;
             }else{
              flag = 0;
             }

             if (flag == 1){
                 temp_neigh = ebins[i].num_neighbors;
                 //cout << temp_neigh << " " << ibin1 << " " << jbin1 << " " << kbin1 << endl;
                 //cout << temp_neigh << endl;
                 if (temp_neigh > 5) {
                     cout << "more than six neighbors" << endl;
		     exit(EXIT_FAILURE);
                 }
                 ebins[i].neigh_indexes[temp_neigh] = j+1;
                 ebins[i].num_neighbors = ebins[i].num_neighbors + 1;

                 temp_neigh = ebins[j].num_neighbors;
                 ebins[j].neigh_indexes[temp_neigh] = i+1;
                 ebins[j].num_neighbors = ebins[j].num_neighbors + 1;
                 
             }
             
        }
   }

   cout << "Done finding empty bins " << endl;
   identify_voids(act_frame1, natoms1,grid_dr1,empty_bin_count, ebins,grid_atom_count1,grid_list1,at_list1,lx1,ly1,lz1,nx1,ny1,nz1);
}
