#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include "atom.h"
#include "functions.h"
#include "settings.h"

//using namespace std;

void make_connection_table(ofstream &file1,int nframes1, int natoms1, atom **at_list1, 
                           double **lx1, double **ly1, double **lz1,double *alpha1, double *beta1, double *gama1,
                           int ****grid_list1, int ***grid_atom_count1,int nx1,int ny1,int nz1,settings settings_list1){
                                  

          int i,j,k,m,n;
          int p,q,r;
          int plow,phigh,qlow,qhigh,rlow,rhigh;
          int nid,bin_pop;
          int num2,count;
          int i1, i2;
          //int neighbours[15];
          int no_neighbours;
          int check_flag,flag1;
          double dx,dy,dz,dr;
          double x1,y1,z1,x2,y2,z2;
          int xbin1,ybin1,zbin1,xbin2,ybin2,zbin2;
          int bin_indexes[6][3];
          int file_write_flag = 1; //1 means write file1 and file2
          string temp,cond;
          string filename;
          ostringstream convert;
          ofstream file2,file3;

          double cut_off = settings_list1.return_distance_cut_off();
                  
          //cout << "starting connection table " << endl;
          convert << (nframes1+1);
          filename = convert.str();
          //file3.open((filename+"_connection_table.txt").c_str());
          if (file_write_flag == 1) file2.open("temp_connection_table.txt");
          //identify atom type number
          for (i=0;i<1;i++){
              //write connection table
              if (file_write_flag == 1){
                  file1 << "Frame no " << nframes1+1 << endl;
                  file2 << "Frame no " << nframes1+1 << endl;
              }
              
              //cout << "for fram " << i << " last atom has atom no " << at_no[natoms1-1] << endl;
              for (j=0;j<natoms1;j++){
                  //cout << "starting atom " << j+1 << endl;
                  check_flag = 0;
                  
                  //initialize all neighbour and bond order related information
                  no_neighbours = 0;
                  
                  temp.clear();
                  temp = at_list1[i][j].return_atomname();
                  
                  x1 = at_list1[i][j].return_xcord();
                  y1 = at_list1[i][j].return_ycord();
                  z1 = at_list1[i][j].return_zcord();

                  xbin1 = at_list1[i][j].return_xbin();
                  ybin1 = at_list1[i][j].return_ybin();
                  zbin1 = at_list1[i][j].return_zbin();

                  plow = qlow = rlow = -1;
                  phigh = qhigh = rhigh = 2;
                  
                  //following part is needed becuase some cells near pbs can have only 1-2 atoms
                  if (xbin1 == 0 || xbin1 == 1)         plow = -2;
                  if (xbin1 == nx1-1 || xbin1 == nx1-2) phigh = 3;
                  if (ybin1 == 0 || ybin1 == 1)         qlow = -2;
                  if (ybin1 == ny1-1 || ybin1 == ny1-2) qhigh = 3;
                  if (zbin1 == 0 || zbin1 == 1)         rlow = -2;
                  if (zbin1 == nz1-1 || zbin1 == nz1-2) rhigh = 3;
                  count = 0;
                  if (xbin1 == 0 || xbin1 == nx1-1 || ybin1 == 0 || ybin1 == ny1-1 || zbin1 == 0 || zbin1 == nz1-1) check_flag = 1;
                  
                  for (p=plow;p<phigh;p++){
                       xbin2 = xbin1 + p;
                       if (xbin2 < 0) xbin2 = nx1 + xbin2;
                       if (xbin2 > nx1-1) xbin2 = xbin2 - nx1;
                       for (q=qlow;q<qhigh;q++){
                            ybin2 = ybin1 + q;
                            if (ybin2 < 0) ybin2 = ny1 + ybin2;
                            if (ybin2 > ny1-1) ybin2 = ybin2 - ny1;
                           for (r=rlow;r<rhigh;r++){
                                zbin2 = zbin1 + r;
                                if (zbin2 < 0) zbin2 = nz1 + zbin2 ;
                                if (zbin2 > nz1-1) zbin2 = zbin2 - nz1;
                                bin_pop = grid_atom_count1[xbin2][ybin2][zbin2];
                                for (n=0;n<bin_pop;n++){
                                     //if (j == 1841) cout << xbin2 << " " << ybin2 << " " << zbin2 << endl;
                                     nid = grid_list1[xbin2][ybin2][zbin2][n];
                                     temp.clear();
                                     temp = at_list1[i][nid].return_atomname();

                                     x2 = at_list1[i][nid].return_xcord();
                                     y2 = at_list1[i][nid].return_ycord();
                                     z2 = at_list1[i][nid].return_zcord();
                                     //if (j == 3965 && nid == 3964) cout << "found the pair" << endl;
                                     dx = x2-x1;
                                     dy = y2-y1;
                                     dz = z2-z1;
                      
                                     //check periodicity
                                     //check for periodic images in x direction
                                     if (dx > lx1[i][2]/2.0){
                                        dx = dx - lx1[i][2];
                                     }
                                     if(dx <= (-1.0)*(lx1[i][2]/2.0)){
                                        dx = dx + lx1[i][2];
                                     }
                                     //check for periodic images in y direction
                                     if (dy > ly1[i][2]/2.0){
                                         dy = dy - ly1[i][2];
                                     }
                                     if(dy <= (-1.0)*(ly1[i][2]/2.0)){
                                         dy = dy + ly1[i][2];
                                     }
                                     //check for periodic images in z direction
                                     if (dz > lz1[i][2]/2.0){
                                         dz = dz - lz1[i][2];
                                     } 
                                     if(dz <= (-1.0)*(lz1[i][2]/2.0)){
                                        dz = dz + lz1[i][2];
                                     }
                      
                                     dr = sqrt(dx*dx + dy*dy + dz*dz);
                                     num2 = 200;
                            
                                     if(dr < cut_off && j != nid){
                                         //if edge cell, then check to make sure this atom is not already in the neighbor list
                                         // this is required for very small grid like only 2 cells in any direction 
                                         flag1 = 0;
                                         if (check_flag == 1){
                                             for (m=0;m < no_neighbours;m++){
                                                  if (at_list1[i][j].neigh_indexes[m] == nid+1) flag1 = 1;
                                             }
                                         }
                                         if (flag1 == 0){
                                             count = no_neighbours;
                                             at_list1[i][j].set_neigh_index(count,nid+1);
                                             //neighbours[count] = nid+1;
                                             no_neighbours = no_neighbours + 1;
                                             if (no_neighbours > max_neigh){
                                                 cout << "atom " << j+1 << " has more than 16 neighbors. Change neighbors array size " << endl;
                                                 exit(EXIT_FAILURE);
                                             }
                            
                                         }
                                     } //bond order if ends here
                            
                                }
                           }
                       }
                  }
               at_list1[i][j].set_num_neighbors(no_neighbours);   
               if (file_write_flag == 1){
                  file1 << setw(8) << j+1 << " " << setw(3) <<at_list1[i][j].num_neighbors << " " ;
                  file2 << setw(8) << j+1 << " " << setw(3) <<at_list1[i][j].num_neighbors << " " ;
                  //file3 << setw(8) << j+1 << " " << setw(3) <<no_neighbours << " " << endl ;
                  for(k=0;k<no_neighbours;k++){
                                              n = at_list1[i][j].neigh_indexes[k];
                                              //file1 << setw(2) << at_list1[i][n-1].return_atomname() << setw(5) 
                                              file1 << setw(2) << at_list1[i][n-1].return_atomname() << " " << setw(8) 
                                              << n << " " ;
                                              file2 << setw(2) << at_list1[i][n-1].return_atomname() << " " << setw(8) 
                                              << n << " " ;
                                              }
                  file1 << endl;   
                  file2 << endl;
               }
              } // j for loop ends here
          
          
          }//no of frames loop
          
          if (file_write_flag == 1){
              file1.close();
              file2.close();
          }
          //file3.close();
          cout << "finished preparing connection table " << endl;
          
          //cout << at_no[nframes1-1][natoms1-1] << endl;

}
