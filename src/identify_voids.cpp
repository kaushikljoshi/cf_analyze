#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include "atom.h"
#include "voids.h"
#include "functions.h"

using namespace std;
const int max_neigh_size = 5000;
void identify_void_edges(int no_voids,voids *void_list1, emptyBin *bin_list1,int ***grid_atom_count1);
void find_void_surface_atoms(int no_voids, int natoms, double *grid_dr1,voids *void_list1, emptyBin *bin_list1, int ***grid_atom_count1, int ****grid_list1, atom **at_list1, double **lx1, double **ly1, double **lz1, int nx1, int ny1, int nz1);
bool sortcol( const vector<double>& v1,const vector<double>& v2);
void surface_atoms(int dirn,double *grid_dr1, double xc, double yc, double zc,atom **at_list1, voids *void_list1,double **lx1, double **ly1, double **lz1);
void find_max_pore_size(int no_voids,voids *void_list1,double *grid_dr1, double **lx1, double **ly1, double **lz1, int ***grid_atom_count1);
void check_periodic(double& dr, double **len);
void set_array(double array[], double val1, double val2,double len);
void get_distribution(const vector<double>& dist_list1,const vector<double>& dist_list2,double dr,string filename);
bool check_pbc_crossing(const vector< double > temp_vector,double cell_low,double cell_high);

void identify_voids(int act_frame, int natoms1, double *grid_dr1,int no_empty_bins, emptyBin *bin_list1,int ***grid_atom_count1,int ****grid_list1, atom **at_list1, double **lx1, double **ly1, double **lz1, int nx1, int ny1, int nz1){
     cout << "starting void identification " << endl;
     int i,j,k;
     int l,m,n;
     int count1,count2,tmp_no,no_cycles;
     int no1,no2,no3,no4;
     int no_types_voids[max_void_types];         
     double void_vols[max_void_types];
     int pop_voids[max_void_types];
     int no_voids,isolated_bins;
     //int temp_neighbours[natoms1];  //20 is assumed as maximum no of neighbours
     int temp_neighbors,temp_bins;
     double temp_vol;
     int first_neighbours[max_neigh_size];
     int second_neighbours[max_neigh_size];
     string sub,des;
     string flag;
     string s1,s2;
     string lines;
     string temp_lines[50];
     int flags[no_empty_bins]; //0 means yes, 1 means NO
     cout << "starting molecule identification " << endl;
     voids *void_list;
     vector < vector < int > > void_bin_list;
     vector < int > temp_void_bins;
     ifstream myfile;
     ofstream myfile1,myfile2;
     istringstream iss;
     
     //initialize the void status of every bin to No
     for(j=0;j<no_empty_bins;j++){
        bin_list1[j].set_bin_state(1);
     }
     
     //cout << "finished setting atom state" << endl;
     //myfile.open("temp_connection_table.txt");
     
     no_voids = 0;
     isolated_bins = 0;
     
     //identify isolated voids, meaning bins with 0 empty neighbors. Isolated empty bins will not be considered as voids 
     for(j=0;j<no_empty_bins;j++){
         if(bin_list1[j].num_neighbors == 0){
            bin_list1[j].set_bin_state(0);
            //no_voids = no_voids + 1;
            bin_list1[j].set_void_no(-1);
            isolated_bins = isolated_bins + 1;
         }
     }
     //start looking for bins of every void
     cout << "Number of unimolecular pieces are " << isolated_bins << endl;
     for(j=0;j<no_empty_bins;j++){
         temp_neighbors = bin_list1[j].num_neighbors;
         if(bin_list1[j].return_bin_state() == 1 && temp_neighbors != 0){
            no_voids = no_voids + 1;
            //identify first immediate neighbours
            no1= temp_neighbors;
            //cout << "No of neighbours are " << no1 << endl;
            //count1 = 0;
            for(k=0;k<no1;k++){
                first_neighbours[k] = bin_list1[j].neigh_indexes[k]-1;
                //cout << "First neighbors " << first_neighbours[k] << endl;
            }
                                                                                         
            flag = "Yes";
            no_cycles = 0;
            while(flag == "Yes"){
                  for(k=0;k<no_empty_bins;k++){
                      flags[k] = 1;
                  }
                  count2 = 0;
                  count1 = 0;
                  //cout << "for cycle " << no_cycles + 1 << " count1 no is " << count1 << endl;
                  //cout << no_cycles+1 << " started " << endl;
                  for(k=0;k<no1;k++){
                      no2= first_neighbours[k];
                      if(bin_list1[no2].return_bin_state() == 1){
                         bin_list1[no2].set_void_no(no_voids);
                         bin_list1[no2].set_bin_state(0);
                      }else{
                         count2 = count2 + 1;
                      }
                                                                                                            
                      no3 = bin_list1[no2].num_neighbors; //temp_neighbours[no2];
                      //cout << "k with no2 " << k << " with no2 equal to " << no2 << endl;
                                                                                                            
                      for(l=0;l<no3;l++){
                          no4 = bin_list1[no2].neigh_indexes[l]-1; //neighbours_index[no2][l] - 1;
                          if(bin_list1[no4].return_bin_state() == 1 && flags[no4] == 1){
                             second_neighbours[count1] = bin_list1[no2].neigh_indexes[l]-1; //neighbours_index[no2][l] - 1;
                             flags[no4] = 0;
                             count1 = count1 + 1;
                          }
                       }
                       //cout << "k is " << k << "  no2 is " << no2 << "  no4 is " << no4 << endl;
                       //cout << at_list1[i][no2].return_atomname() << " " << no2 << " done" << endl;
                       //cout << count1 << endl;
                  }
                  //cout << no_cycles + 1 << " reached until here " << endl;
                  //cout << count1 << endl;
                  if (count1 > max_neigh_size ){
                      cout << "Increase size of second_neighbour array" << endl;
                      //system("Pause");
                      exit(1);
                  }
                  if(count2 == no1){
                     flag = "No";
                     //break;
                  }else{//transfer second neighbour into the first neighbour
                     flag = "Yes";
                     no1 = count1;
                     for(k=0;k<no1;k++){
                         first_neighbours[k] = second_neighbours[k];
                     }
                  }
                  no_cycles = no_cycles + 1;
                  //cout << "Cycle number is " << no_cycles << " with count1 number "<< count1 << endl;
            }
                                                    
            //cout << "atom " << j << " is done " << endl;
         }//end of if 
     }//end of natoms1 for loop
                             
     cout << "No voids in frame " << act_frame << " is " << no_voids << endl;
     
     //myfile.close();
     
     //allocate array for voids
     void_list = new voids [no_voids];

     //identify size of each void
     for (i=0;i<no_voids;i++){
          temp_bins = 0;
          temp_vol  = 0.0;
          temp_void_bins.clear();
          for (j=0;j<no_empty_bins;j++){
               if (bin_list1[j].return_void_no() == i+1){
                   temp_bins = temp_bins + 1;
                   temp_vol = temp_vol + bin_list1[j].return_bin_vol();
                   void_list[i].vbin_list.push_back(j);
                   //void_list[i].add_bin(j);
                   temp_void_bins.push_back(j);
               }
          }
          void_list[i].set_void_no_bins(temp_bins);
          void_list[i].set_void_vol(temp_vol);
          void_bin_list.push_back(temp_void_bins);
     }
     
     //void population analysis
     for(i=0;i<max_void_types;i++){
         pop_voids[i]=0;
     }

     count1=0;
     for(i=0;i<no_voids;i++){
         s2= "N";
         if(i==0){
            void_vols[count1] = void_list[i].return_void_vol();
            no_types_voids[count1]= void_list[i].return_no_bins();
            pop_voids[count1] = pop_voids[count1] + 1;
            count1 = count1 + 1;
            s2 = "Y";
          }else{
           temp_bins = void_list[i].return_no_bins();
           //temp_vol = void_list[i].return_void_vol();
           for(j=0;j<count1;j++){
               //if(temp_vol == void_vols[j]){
               if(temp_bins == no_types_voids[j]){
                  pop_voids[j] = pop_voids[j] + 1;
                  s2 = "Y";
                  break;
               }
                                                                   }
           }
          if(s2 == "N"){
             void_vols[count1] = void_list[i].return_void_vol();
             no_types_voids[count1] = void_list[i].return_no_bins();
             pop_voids[count1] = pop_voids[count1] + 1;
             count1 = count1 + 1;
          }
     }
    
    if (act_frame == 0){
        myfile1.open("voids.out");
    }else{
        myfile1.open("voids.out",ios::app);
    }
    for(i=0;i<count1;i++){
        //cout << pop_mols[i] << "  " << no_types_mols[i] << endl;
        myfile1 << act_frame << "   " << pop_voids[i] << "  x  " <<  void_vols[i] << endl;
        //cout  << act_frame << "   " << pop_voids[i] << "  x  " <<  void_vols[i] << endl;
    }
    myfile1 << "Total number of voids:           " << no_voids << endl;
    myfile1 << "Total number of emptyBins:         " << no_empty_bins << endl;
    double t_vol = 0.0;
    for(i=0;i<count1;i++){
        t_vol = void_vols[i]*pop_voids[i] + t_vol;
    }
    myfile1 << "Total void vol:   " << t_vol << endl;
    myfile1.close();

    //Following is for debugging purposes only
    /*myfile2.open("combined_empty_bins.xyz");
    myfile2 << no_empty_bins-isolated_bins << endl;
    myfile2 << endl;
    for (i=0;i<count1;i++){
        std::string filename = std::to_string(no_types_voids[i]) +".xyz";
        myfile1.open(filename);
        myfile1 << pop_voids[i]*no_types_voids[i] << endl;
        myfile1 << endl;
        for (j=0;j<no_voids;j++){
             if (no_types_voids[i] == void_bin_list[j].size()){
                 for (k=0;k<void_bin_list[j].size();k++){
                      temp_bins = void_bin_list[j][k];
                      myfile1 << "C " << bin_list1[temp_bins].return_xc() << " " << bin_list1[temp_bins].return_yc() << " " << bin_list1[temp_bins].return_zc() << endl;
                      myfile2 << "C " << bin_list1[temp_bins].return_xc() << " " << bin_list1[temp_bins].return_yc() << " " << bin_list1[temp_bins].return_zc() << endl;
                 }
             }
        }
        myfile1.close();
    }
    myfile2.close();*/
   
    identify_void_edges(no_voids,void_list,bin_list1,grid_atom_count1);
    find_void_surface_atoms(no_voids,natoms1,grid_dr1,void_list,bin_list1,grid_atom_count1,grid_list1,at_list1,lx1,ly1,lz1,nx1,ny1,nz1);
    
     delete[] void_list;
}

//########################################################################################################
//Following function identities edges/edge bins of each void
void identify_void_edges(int no_voids,voids *void_list1,emptyBin *bin_list1,int ***grid_atom_count1){
     int i,j,k;
     int ibin,jbin,kbin;
     int ebin_no,no1,no2;
     int flag;
     int count,global_edge_count;
     ofstream myfile,myfile1,myfile2;

     global_edge_count = 0;
     for (i=0;i<no_voids;i++){
          for (j=0;j<void_list1[i].vbin_list.size();j++){
               ebin_no = void_list1[i].vbin_list[j];
               no1 = bin_list1[ebin_no].num_neighbors;
               count = no1;
               /*for (k=0;k<no1;k++){
                    no2 = bin_list1[ebin_no].neigh_indexes[k]-1;
                    ibin = bin_list1[no2].return_xbin();
                    jbin = bin_list1[no2].return_ybin();
                    kbin = bin_list1[no2].return_zbin();
                    if(grid_atom_count1[ibin][jbin][kbin] == 0) count = count + 1;
               }*/
               if (count != 6){
                   void_list1[i].edge_list.push_back(ebin_no);
                   global_edge_count = global_edge_count + 1;
               }
          }
     }
    /*myfile.open("combined_void_edges.xyz");
    myfile << global_edge_count << endl;
    myfile << endl;
    for (i=0;i<no_voids;i++){
         //std::string filename = std::to_string(void_list1[i].return_void_vol()) +"_" + std::to_string(i+1) + "_edge.xyz";
         std::string filename = std::to_string(i+1) +"_"+ std::to_string(void_list1[i].edge_list.size())+"_edge_bins.xyz";
         std::string filename1 = std::to_string(i+1) + "_"+ std::to_string(void_list1[i].vbin_list.size())+ "_bins.xyz";
         myfile1.open(filename);
         myfile2.open(filename1);
         no1 = void_list1[i].edge_list.size();
         myfile1 << no1 << endl;
         myfile1 << endl;
         myfile2 << void_list1[i].vbin_list.size() << endl;
         myfile2 << endl;
         for (j=0;j<void_list1[i].edge_list.size();j++){
              ebin_no = void_list1[i].edge_list[j];
              myfile << "N " << bin_list1[ebin_no].return_xc() << " " << bin_list1[ebin_no].return_yc() << " " << bin_list1[ebin_no].return_zc() << " " << i+1 << endl;
              myfile1 << "N " << bin_list1[ebin_no].return_xc() << " " << bin_list1[ebin_no].return_yc() << " " << bin_list1[ebin_no].return_zc() << endl;
         }
         for (j=0;j<void_list1[i].vbin_list.size();j++){
              ebin_no = void_list1[i].vbin_list[j];
              myfile2 << "C " << bin_list1[ebin_no].return_xc() << " " << bin_list1[ebin_no].return_yc() << " " << bin_list1[ebin_no].return_zc() << endl;
         }
         myfile1.close();
         myfile2.close();
    }
    myfile.close();*/
}

//#############################################################################################################
void find_void_surface_atoms(int no_voids, int natoms, double *grid_dr1,voids *void_list1, emptyBin *bin_list1, int ***grid_atom_count1, int ****grid_list1, atom **at_list1, double **lx1, double **ly1, double **lz1, int nx1, int ny1, int nz1){
     int i,k,j,l,m,n;
     int ibin1,jbin1,kbin1;
     int ilo,ihi,jlo,jhi,klo,khi;
     int ibin,jbin,kbin;
     double xc,yc,zc;
     double xtemp,ytemp,ztemp;
     double x1,y1,z1,x2,y2,z2;
     double dx,dy,dz,dr,temp_dr;
     vector< vector <double> > at_dist_list;
     vector< double > temp_list;
     vector< double > temp_xlist;
     vector< double > temp_ylist;
     vector< double > temp_zlist;
     vector< double > dist_list1;
     vector< double > dist_list2;
     vector< int > temp_neigh_list;
     int edge_bin_count,edge_bin_no,ebin_no;
     int no1,no2,no3,at_id,temp1;
     int at_tag[natoms];
     int at_num_cutoff = 5;
     int total_edge_atoms = 0;
     int total_edge_faces = 0;
     int face_flag;
     string filename;
     ofstream myfile,myfile1;

     for (i=0;i<natoms;i++){
          at_tag[i] = 0;
     }

     for (i=0;i<no_voids;i++){
          temp_xlist.clear();
          temp_ylist.clear();
          temp_zlist.clear();
          edge_bin_count = void_list1[i].edge_list.size();
          for (j=0;j<edge_bin_count;j++){
               edge_bin_no = void_list1[i].edge_list[j]; //this number correponds to index in empty bin list
               xc  = bin_list1[edge_bin_no].return_xc();
               yc  = bin_list1[edge_bin_no].return_yc();
               zc  = bin_list1[edge_bin_no].return_zc();
               ibin1 = bin_list1[edge_bin_no].return_xbin();
               jbin1 = bin_list1[edge_bin_no].return_ybin();
               kbin1 = bin_list1[edge_bin_no].return_zbin();
               
               ilo = ibin1-1;
               ihi = ibin1+1;
               jlo = jbin1-1;
               jhi = jbin1+1;
               klo = kbin1-1;
               khi = kbin1+1;

               if (ilo < 0)     ilo = nx1-1;
               if (ihi > nx1-1) ihi = 0;

               if (jlo < 0)     jlo = ny1-1;
               if (jhi > ny1-1) jhi = 0;

               if (klo < 0)     klo = nz1-1;
               if (khi > nz1-1) khi = 0;

               temp_neigh_list.clear();
                 
               temp_neigh_list.push_back(ilo);
               temp_neigh_list.push_back(jbin1);
               temp_neigh_list.push_back(kbin1);
                
               temp_neigh_list.push_back(ihi);
               temp_neigh_list.push_back(jbin1);
               temp_neigh_list.push_back(kbin1);

               temp_neigh_list.push_back(ibin1);
               temp_neigh_list.push_back(jlo);
               temp_neigh_list.push_back(kbin1);

               temp_neigh_list.push_back(ibin1);
               temp_neigh_list.push_back(jhi);
               temp_neigh_list.push_back(kbin1);

               temp_neigh_list.push_back(ibin1);
               temp_neigh_list.push_back(jbin1);
               temp_neigh_list.push_back(klo);

               temp_neigh_list.push_back(ibin1);
               temp_neigh_list.push_back(jbin1);
               temp_neigh_list.push_back(khi);

               //no1 = bin_list1[edge_bin_no].num_neighbors;
               no1 = 6;
               for (k=0;k<no1;k++){
                    //no2 = bin_list1[edge_bin_no].neigh_indexes[k]-1;
                    //ibin = bin_list1[no2].return_xbin();
                    //jbin = bin_list1[no2].return_ybin();
                    //kbin = bin_list1[no2].return_zbin();
                    ibin = temp_neigh_list[k*3];
                    jbin = temp_neigh_list[k*3+1];
                    kbin = temp_neigh_list[k*3+2];
                    no3  = grid_atom_count1[ibin][jbin][kbin];
                    if (no3 != 0){
                        /*xtemp = grid_dr1*(ibin+1);
                        ytemp = grid_dr1*(jbin+1);
                        ztemp = grid_dr1*(kbin+1);*/
                        xtemp = xc;
                        ytemp = yc;
                        ztemp = zc;
                        if (k == 0) {
                            xtemp = xc - grid_dr1[0]*0.5;
                            if (xtemp < lx1[0][0] ) xtemp = lx1[0][0];
                        }
                        if (k == 1) {
                            xtemp = xc + grid_dr1[0]*0.5;
                            if (xtemp > lx1[0][1]) xtemp = lx1[0][1];
                        }
                        if (k == 2){
                            ytemp = yc - grid_dr1[1]*0.5;
                            if (ytemp < ly1[0][0]) ytemp = ly1[0][0];
                        }
                        if (k == 3) {
                            ytemp = yc + grid_dr1[1]*0.5;
                            if (ytemp > ly1[0][1]) ytemp = ly1[0][1];
                        }
                        if (k == 4){
                            ztemp = zc - grid_dr1[2]*0.5;
                            if (ztemp < lz1[0][0]) ztemp = lz1[0][0];
                        }
                        if (k == 5) {
                            ztemp = zc + grid_dr1[2]*0.5;
                            if (ztemp > lz1[0][1]) ztemp = lz1[0][1];
                        }
                        face_flag = 0;
                        //check if the face is already in the list
                        for (l=0;l<void_list1[i].face_xc.size();l++){
                             if (xtemp == void_list1[i].face_xc[l] && ytemp == void_list1[i].face_yc[l] && ztemp == void_list1[i].face_zc[l]){ 
                                 face_flag = 1;
                                 break;
                             }
                        }
                        if (face_flag == 0){
                            void_list1[i].face_xc.push_back(xtemp);
                            void_list1[i].face_yc.push_back(ytemp);
                            void_list1[i].face_zc.push_back(ztemp);
                            total_edge_faces = total_edge_faces + 1;
                        }
                    }
                    at_dist_list.clear();
                    for (l=0;l<no3;l++){
                         temp_list.clear();
                         at_id = grid_list1[ibin][jbin][kbin][l];
                         temp_list.push_back(at_id);
                         xtemp = at_list1[0][at_id].return_xcord();
                         ytemp = at_list1[0][at_id].return_ycord();
                         ztemp = at_list1[0][at_id].return_zcord();
 
                         dx = xtemp-xc;
                         dy = ytemp-yc;
                         dz = ztemp-zc;

			 //check periodicity
                         //check for periodic images in x direction
                         if (k < 2){
                             if (dx > lx1[0][2]/2.0){
                                 dx = dx - lx1[0][2];
                             }
                             if(dx <= (-1.0)*(lx1[0][2]/2.0)){
                                dx = dx + lx1[0][2];
                             }
                             temp_list.push_back(dx);
                         }
                         if (1 < k < 4){
                             //check for periodic images in y direction
                             if (dy > ly1[0][2]/2.0){
                                 dy = dy - ly1[0][2];
                             }
                             if(dy <= (-1.0)*(ly1[0][2]/2.0)){
                                dy = dy + ly1[0][2];
                             }
                             temp_list.push_back(dy);
                         }
                         if (3 < k < 6){
                             //check for periodic images in z direction
                             if (dz > lz1[0][2]/2.0){
                                 dz = dz - lz1[0][2];
                             } 
                             if (dz <= (-1.0)*(lz1[0][2]/2.0)){
                                 dz = dz + lz1[0][2];
                             }
                             temp_list.push_back(dz);
                         }
                      
                         //dr = sqrt(dx*dx + dy*dy + dz*dz);
                         //temp_list.push_back(dr);
                         at_dist_list.push_back(temp_list);
                    }//loop for atoms in neighboring bin ends here
                    if (at_dist_list.size() > 0){
                        sort(at_dist_list.begin(),at_dist_list.end(),sortcol);
                        at_id = at_dist_list[0][0];
                        if (at_tag[at_id] == 0){
                            void_list1[i].edge_atom_list.push_back(at_id);
                            xtemp = at_list1[0][at_id].return_xcord();
                            ytemp = at_list1[0][at_id].return_ycord();
                            ztemp = at_list1[0][at_id].return_zcord();
                            temp_xlist.push_back(xtemp);
                            temp_ylist.push_back(ytemp);
                            temp_zlist.push_back(ztemp);
                            at_tag[at_id] = 1;
                            total_edge_atoms = total_edge_atoms + 1;
                        }
                    }
               }// k-loop ends here
          }
     } //i-loop ends here
     
     cout << "Done finding surface atoms " << endl;
     //maximum size of pore for each void
     find_max_pore_size(no_voids,void_list1,grid_dr1,lx1,ly1,lz1,grid_atom_count1);
     cout << "Done finding maximum pore size " << endl;
     /*myfile.open("voids_size.dat");
     for (i=0;i<no_voids;i++){
          //myfile << i+1 << " " << void_list1[i].max_len << " " << void_list1[i].return_void_vol() << endl;
          myfile << i+1 << " " << void_list1[i].max_len << " " << void_list1[i].return_void_vol(); 
          myfile << " " << void_list1[i].face_xc[void_list1[i].index1] << " " << void_list1[i].face_xc[void_list1[i].index2] ;
          myfile << " " << void_list1[i].face_yc[void_list1[i].index1] << " " << void_list1[i].face_yc[void_list1[i].index2] ;
          myfile << " " << void_list1[i].face_zc[void_list1[i].index1] << " " << void_list1[i].face_zc[void_list1[i].index2] << endl;
     }
     myfile.close();*/
     //create vectors for distribution 
     for (i=0;i<no_voids;i++){
          dist_list1.push_back(void_list1[i].max_len);
          dist_list2.push_back(void_list1[i].return_void_vol());
     }
     double max_dr;
     max_dr = grid_dr1[0];
     for (i=1;i<3;i++){
         if (grid_dr1[i] > max_dr) max_dr = grid_dr1[i]; 
     }
     get_distribution(dist_list1,dist_list2,2*max_dr,"pore_vol");

     cout << "Total number of surface atoms " << total_edge_atoms << endl;
     cout << "done identifying surface atoms" << endl;
      
     //write edge atoms files. One file contains all edge atoms of all voids
     //another file contains edge atoms of each void. 
     /*myfile.open("combined_void_edge_atoms.xyz");
     myfile << total_edge_atoms << endl;
     myfile << endl;
     for (i=0;i<no_voids;i++){
          std::string filename = std::to_string(void_list1[i].return_void_vol()) +"_" + std::to_string(i+1) +"_edge_atoms.xyz";
          myfile1.open(filename);
          no1 = void_list1[i].edge_atom_list.size();
          myfile1 << no1 << endl;
          myfile1 << endl;
          for (j=0;j<no1;j++){
               at_id = void_list1[i].edge_atom_list[j];
               myfile1 << at_list1[0][at_id].return_atomname() << " " << at_list1[0][at_id].return_xcord() << " " << at_list1[0][at_id].return_ycord() << " " << at_list1[0][at_id].return_zcord() << endl;
               myfile << at_list1[0][at_id].return_atomname() << " " << at_list1[0][at_id].return_xcord() << " " << at_list1[0][at_id].return_ycord() << " " << at_list1[0][at_id].return_zcord() << endl;
          }
          myfile1.close();
     }
     myfile.close();*/

     //write surface coordinates files. One file contains all surfaces of all voids
     //another file contains surfaces of each void. 
     myfile.open("combined_void_faces.xyz");
     myfile << total_edge_faces << endl;
     myfile << endl;
     for (i=0;i<no_voids;i++){
          //check to see if void crosses periodic boundary. if yes, remap-coordinates
          no1 = void_list1[i].face_xc.size();
          temp_xlist.clear();
          temp_ylist.clear();
          temp_zlist.clear();
          //check for pbc crossing in x-direction
          for (j=0;j<no1;j++){
               temp_xlist.push_back(void_list1[i].face_xc[j]);
          }
          if (check_pbc_crossing(temp_xlist,lx1[0][0],lx1[0][1])){
              for(j=0;j<no1;j++){
                  if(temp_xlist[j]<(lx1[0][0]+lx1[0][2]/2.0)) temp_xlist[j] = temp_xlist[j]+lx1[0][2];
              }
          }
          //check for pbc crossing in y-direction
          for (j=0;j<no1;j++){
               temp_ylist.push_back(void_list1[i].face_yc[j]);
          }
          if (check_pbc_crossing(temp_ylist,ly1[0][0],ly1[0][1])){
              for(j=0;j<no1;j++){
                  if(temp_ylist[j]<(ly1[0][0]+ly1[0][2]/2.0)) temp_ylist[j] = temp_ylist[j]+ly1[0][2];
              }
          }
          //check for pbc crossing in z-direction
          for (j=0;j<no1;j++){
               temp_zlist.push_back(void_list1[i].face_zc[j]);
          }
          if (check_pbc_crossing(temp_zlist,lz1[0][0],lz1[0][1])){
              for(j=0;j<no1;j++){
                  if(temp_zlist[j]<(lz1[0][0]+lz1[0][2]/2.0)) temp_zlist[j] = temp_zlist[j]+lz1[0][2];
              }
          }
          std::string filename = std::to_string(i+1) + "_"+ std::to_string(void_list1[i].vbin_list.size()) + "_" + std::to_string(void_list1[i].return_void_vol()) + "_surface.xyz";
          //myfile1.open(filename);
          //myfile1 << no1 << endl;
          //myfile1 << endl;
          for (j=0;j<no1;j++){
               myfile << "S " << void_list1[i].face_xc[j] << " " << void_list1[i].face_yc[j] << " " << void_list1[i].face_zc[j] << endl;
               //myfile1 << "O " << void_list1[i].face_xc[j] << " " << void_list1[i].face_yc[j] << " " << void_list1[i].face_zc[j] << endl;
               //myfile1 << "S " << temp_xlist[j] << " " << temp_ylist[j] << " " << temp_zlist[j] << endl;
          }
          //myfile1.close();
     }
     myfile.close();
     
}

bool sortcol( const vector<double>& v1,const vector<double>& v2){
     return v1[1] < v2[1];
}

void surface_atoms(int dirn, double *grid_dr1, double xc, double yc, double zc,atom **at_list1, voids *void_list1,double **lx1, double **ly1, double **lz1){
   int i,j,k;
   int n1,n2;
   double num1,num2;
   double dr1,dr2;
   double small_dr = 1.0;
   

}


void find_max_pore_size(int no_voids,voids *void_list1,double *grid_dr1, double **lx1, double **ly1, double **lz1, int ***grid_atom_count1){
     int i,j,k;
     int l,m,n;
     double dx,dy,dz,dr;
     double x1,y1,z1,x2,y2,z2;
     double xint,yint,zint;
     int xcount,ycount,zcount;
     int xbin,ybin,zbin;
     int n1,n2,empty_flag;
     double xsweep[3],ysweep[3],zsweep[3];
     vector< vector <double> > at_dist_list;
     vector< double > temp_list;
     
     for (i=0;i<no_voids;i++){
          n1 = void_list1[i].face_xc.size();
          void_list1[i].max_len = 0.0;
          for (j=0;j<void_list1[i].face_xc.size();j++){
               x1 = void_list1[i].face_xc[j];
               y1 = void_list1[i].face_yc[j];
               z1 = void_list1[i].face_zc[j];
               for (k=j+1;k<void_list1[i].face_xc.size();k++){
                    x2 = void_list1[i].face_xc[k];
                    y2 = void_list1[i].face_yc[k];
                    z2 = void_list1[i].face_zc[k];

                    dx = x2 - x1;
                    dy = y2 - y1;
                    dz = z2 - z1;

                    check_periodic(dx,lx1);
                    check_periodic(dy,ly1);
                    check_periodic(dz,lz1);
                      
                    dr = sqrt(dx*dx+dy*dy+dz*dz);
                    if (void_list1[i].max_len < dr){
                        //cout << "starting for void " << i+1 << endl;
                        //cout << i+1 << " " << x1 << " " << y1 << " " << z1 << endl;
                        //cout << i+1 << " " << x2 << " " << y2 << " " << z2 << endl;
                        empty_flag = 0;
                        set_array(xsweep,x1,x2,lx1[0][2]);
                        set_array(ysweep,y1,y2,ly1[0][2]);
                        set_array(zsweep,z1,z2,lz1[0][2]);
                        if (xsweep[2] == 1){
                            xcount = int((xsweep[1]-xsweep[0])/(grid_dr1[0]));
                        }else{
                            xcount = int((xsweep[0]+lx1[0][2]-xsweep[1])/(grid_dr1[0]));
                        }
                        if (ysweep[2] == 1){
                            ycount = int((ysweep[1]-ysweep[0])/(grid_dr1[1]));
                        }else{
                            ycount = int((ysweep[0]+ly1[0][2]-ysweep[1])/(grid_dr1[1]));
                        }
                        if (zsweep[2] == 1){
                            zcount = int((zsweep[1]-zsweep[0])/(grid_dr1[2]));
                        }else{
                            zcount = int((zsweep[0]+lz1[0][2]-zsweep[1])/(grid_dr1[2]));
                        }
                        
                        //cout << xcount << " " << ycount << " " << zcount << endl;
                        xint = xsweep[0]+grid_dr1[0]*0.5*xsweep[2];
                        for (l=0;l<xcount;l++){
                             xint = xint + grid_dr1[0]*xsweep[2];
                             if (xint < lx1[0][0]) xint = xint + lx1[0][2];
                             if (xint > lx1[0][1]) xint = lx1[0][1];
                             if ((xint-lx1[0][0]) == lx1[0][2]){
                                 xbin = int((xint-lx1[0][0])/grid_dr1[0])-1;
                             }else{
                                 xbin = int((xint-lx1[0][0])/grid_dr1[0]);
                             }
                             yint = ysweep[0]+grid_dr1[1]*0.5*ysweep[2];
                             for (m=0;m<ycount;m++){
                                  yint = yint + grid_dr1[1]*ysweep[2];
                                  if (yint < ly1[0][0]) yint = yint + ly1[0][2];
                                  if (yint > ly1[0][1]) yint = ly1[0][1];
                                  if ((yint-ly1[0][0]) == ly1[0][2]){
                                      ybin = int((yint-ly1[0][0])/grid_dr1[1])-1;
                                  }else{
                                      ybin = int((yint-ly1[0][0])/grid_dr1[1]);
                                  }    
                                  zint = zsweep[0]+grid_dr1[2]*0.5*zsweep[2];
                                  for (n=0;n<zcount;n++){
                                       zint = zint + grid_dr1[2]*zsweep[2];
                                       if (zint < lz1[0][0]) zint = zint + lz1[0][2];
                                       if (zint > lz1[0][1]) zint = lz1[0][1];
                                       if ((zint-lz1[0][0]) == lz1[0][2]){
                                           zbin = int((zint-lz1[0][0])/grid_dr1[2])-1;
                                       }else{
                                           zbin = int((zint-lz1[0][0])/grid_dr1[2]);
                                       }    
                                       //zbin = int((zint-lz1[0][0])/grid_dr1[2]);
                                       //cout << xint << " " << yint << " " << zint << endl;
                                       //cout << xbin << " " << ybin << " " << zbin << endl;
                                       if (grid_atom_count1[xbin][ybin][zbin] > 1){
                                           empty_flag = 1;
                                           break;
                                       }
                                  }//zsweep
                                  if (empty_flag == 1) break;
                             }//ysweep
                             if (empty_flag == 1) break;
                        }//xsweep
                        if (empty_flag == 0){
                            void_list1[i].max_len = dr;
                            void_list1[i].index1 = j;
                            void_list1[i].index2 = k;
                        }
                    }
               }
          }
     }
}

void check_periodic(double& dr, double **len){
     if (dr > len[0][2]/2.0){
         dr = dr - len[0][2];
     } 
     if (dr <= (-1.0)*(len[0][2]/2.0)){
         dr = dr + len[0][2];
     }
}

void set_array(double array[], double val1, double val2,double len){
     if (val1 < val2){
         array[0] = val1;
         array[1] = val2;
     }else{
         array[0] = val2;
         array[1] = val1;
     }
     
    if (abs(val2-val1)> len/2.0) {
        array[2] = -1;
    }else{
        array[2] = 1;
    }
}

void get_distribution(const vector<double>& dist_list1,const vector<double>& dist_list2,double dr,string filename){
     int i,j,k;
     int nbins,bin_no;
     double min,max;
     ofstream myfile;
     std::string filename1 = filename +"_dist.dat";
     myfile.open(filename1);
     myfile << "void_size void_population total_void_vol vol_fraction" << endl; 
     
     min =*min_element(dist_list1.begin(),dist_list1.end());
     max =*max_element(dist_list1.begin(),dist_list1.end());
     if (dr > 0) {
         nbins = int(max)/dr + 1;
     }else{
         cout << "bin size for distribution cannot be zero or negative" << endl;
         exit(1);
     }

     vector <double> list1(nbins,0.0);
     vector <int> list1_count(nbins,0);
     vector <double> list2(nbins,0.0);

     for (i=0;i<dist_list1.size();i++){
          bin_no = int(dist_list1[i]/dr);
          if (bin_no == nbins) bin_no = nbins - 1;
          if (bin_no > nbins-1){
              cout << "Exceeded number of bins " << endl;
              exit(1);
          }
          list1[bin_no] = list1[bin_no] + dist_list1[i];
          list1_count[bin_no] = list1_count[bin_no] + 1;
          list2[bin_no] = list2[bin_no] + dist_list2[i];
     }

     for (i=0;i<nbins;i++){
       //   if (list1_count[i] != 0){
              //myfile << i+1 << "   " << list1_count[i] << "  " << list1[i] << "   " << list2[i] << "   " << list1[i]/list1_count[i] << "   " << list2[i]/accumulate(list2.begin(),list2.end(),0.0)<< endl;
              //myfile << i+1 << "   " << list1_count[i] << "  " << list1[i] << "   " << list2[i] << "   " << (i+0.5)*dr << "   " << list2[i]/accumulate(list2.begin(),list2.end(),0.0)<< endl;
              myfile << (i+0.5)*dr << "  " << list1_count[i] << "  " << list2[i] << "   " << list2[i]/accumulate(list2.begin(),list2.end(),0.0)<< endl;
     //     }
     }

     myfile.close();
}

bool check_pbc_crossing(const vector< double > temp_vector,double cell_low,double cell_high){
     double min,max;

     min =*min_element(temp_vector.begin(),temp_vector.end());
     max =*max_element(temp_vector.begin(),temp_vector.end());
   
     if (abs(min-cell_low)< 4.0 && abs(max-cell_high)< 4.0) return true;
     
     return false;
}
