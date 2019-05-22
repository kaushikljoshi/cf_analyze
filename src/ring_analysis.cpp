#include <iostream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include "functions.h"

using namespace std;

const int max_ring_size = 10; //actual maximum ring allowed will be "max_ring_size -1"
const int max_neighs = 16;
const int ideal_ring_size = 6;

bool check_rings(vector < vector < int > > const &rings, vector < int > const &temp_ring_vertices);
bool check_inter_links(atom **at_list1,vector < int > const &temp_ring_vertices);
bool check_intra_links(atom **at_list1,vector < int > const &temp_ring_vertices);
bool find_cycle(vector < int > &temp_ring_vertices, int curr, atom **at_list1, int bin_dims[3], int bin1[3]);
void write_ring_xyz(int act_frame,int natoms, vector < vector < int > > const &rings, atom **at_list1,double **lx1, double **ly1, double **lz1, settings settings_list1);
void analyze_ring_size(int act_frame, vector< vector < int > > const &rings);
void check_periodic(double &n1,double &n2,double low,double len);

void ring_analysis(int act_frame1, int natoms1, atom **at_list1, int nx1, int ny1, int nz1,
                   double **lx1, double **ly1, double **lz1, settings settings_list1){
   int i,j,k,l,m,n,p,q,r;
   int jhigh,khigh;
   int curr1, curr2, curr3, curr4, curr5, curr6, curr7, curr8;
   int count1,count2,neigh_index;
   int edge_count;
   bool *visited;
   string lines,sub;
   string temp_lines[50];
   int bin_dims[3],bin1[3];
   istringstream iss;
   ifstream ifile1;
   ofstream ofile1;
//   vector < vector < int > > neighbors;
   vector < vector < int > > rings;
//   vector < int > temp_neighbors;
   vector < int > temp_ring_vertices;
   vector < int > sort_ring_vertices;


   cout << "Starting ring analysis" << endl;
   //open output file and write cycle data
   if (act_frame1 == 0){
      ofile1.open("cycles.dat");
   }else{
      ofile1.open("cycles.dat",ios::app);
   }
   bin_dims[0] = nx1;
   bin_dims[1] = ny1;
   bin_dims[2] = nz1;
   
   int ring_count = 0;
   for (i=0;i<natoms1;i++){
        /*if ( i > 1){
           //exit(EXIT_FAILURE);
           break;
        }*/
        if (at_list1[0][i].return_atomname() != "C") continue;
        temp_ring_vertices.clear();
        curr1 = i;
        bin1[0] = at_list1[0][curr1].return_xbin();
        bin1[1] = at_list1[0][curr1].return_ybin();
        bin1[2] = at_list1[0][curr1].return_zbin();
        temp_ring_vertices.push_back(curr1);
        jhigh = at_list1[0][curr1].num_neighbors;
        //cout << "starting ring search for " << i+1 << " atom" << endl;
        for (j = 0;j<jhigh;j++){
             if (temp_ring_vertices.size() > 1)temp_ring_vertices.erase(temp_ring_vertices.begin()+1,temp_ring_vertices.end()); 
             curr2 = at_list1[0][curr1].neigh_indexes[j]-1;
             if (at_list1[0][curr2].return_atomname() != "C") continue;
             if (curr2 > curr1) {
             temp_ring_vertices.push_back(curr2);
             edge_count = 0;
             khigh = at_list1[0][curr2].num_neighbors;
             for (k=0;k<khigh;k++){
                  curr3 = at_list1[0][curr2].neigh_indexes[k]-1; //neighbors[curr2][k];
                  //if (curr3 != curr1 && edge_count < 1){
                  if (curr3 != curr1){
                      if (at_list1[0][curr3].return_atomname() != "C") continue;
                      if (temp_ring_vertices.size() > 2)temp_ring_vertices.erase(temp_ring_vertices.begin()+2,temp_ring_vertices.end());
                     temp_ring_vertices.push_back(curr3);
                     if (find_cycle(temp_ring_vertices, curr3, at_list1, bin_dims,bin1)){
                         //make a deep copy for sorting and comparing
                         sort_ring_vertices.clear();
                         for (l=0;l<temp_ring_vertices.size();l++){
                              sort_ring_vertices.push_back(temp_ring_vertices[l]);
                         }
                         sort(sort_ring_vertices.begin(),sort_ring_vertices.end());
                         if (check_rings(rings,sort_ring_vertices) && check_inter_links(at_list1,temp_ring_vertices)){//make sure newly found ring is really new
                             rings.push_back(sort_ring_vertices);
                             //cout << "found ring of size " << rings[ring_count].size() << endl;
                             ofile1 << ring_count+1 << "  " << rings[ring_count].size() << "  " ;
                             for (int l=0;l<rings[ring_count].size();l++){
                                  //cout << "printing " ;
                                  //cout << rings[ring_count][l] + 1 << " " ;
                                  ofile1 << temp_ring_vertices[l]+1 << " " ;
                             }
                             //cout << "*************************" << endl;
                             ofile1 << endl;
                             edge_count = edge_count + 1;
                             ring_count = ring_count + 1;
                             //Now clear reset temp_ring vertices to first two elements which is original edge
                         }
                     }
                  }
                  //temp_ring_vertices.erase(temp_ring_vertices.begin()+2,temp_ring_vertices.end());
                  //k = k+1;
             }
             //temp_ring_vertices.erase(temp_ring_vertices.begin()+1,temp_ring_vertices.end()); 
           }
        }
        //cout << "Ring search for atom " << i+1 << " done" << endl;
   }

   /*for (i=0;i<rings.size();i++){
        ofile1 << i+1 << "  ";
        for (j=0;j<rings[i].size();j++){
             ofile1 << rings[i][j]+1 << " " ;
        }
        ofile1 << endl;
   }*/
   ofile1.close();
   write_ring_xyz(act_frame1,natoms1,rings,at_list1,lx1,ly1,lz1,settings_list1);
   analyze_ring_size(act_frame1,rings);
   cout << "Finished ring analysis" << endl;
}


bool find_cycle(vector < int > &temp_ring_vertices, int curr, atom **at_list1, int bin_dims[3], int bin1[3]){
   int i,j,k;
   int ihigh,temp_size;
   int next_curr,flag;
   int xbin1,xbin2,ybin1,ybin2,zbin1,zbin2;
   int bin2[3];
   
   int start_vertex = temp_ring_vertices[0];
   temp_size = temp_ring_vertices.size();
   /*for (i=0;i<temp_ring_vertices.size();i++){

        cout << temp_ring_vertices[i]+1 << " ";
   }
   cout << endl;
   cout << "current is " << curr+1 << endl;*/
   ihigh = at_list1[0][curr].num_neighbors;
   for (i=0; i<ihigh;i++){
       next_curr = at_list1[0][curr].neigh_indexes[i]-1; //neighbors[curr][i];
       if (at_list1[0][next_curr].return_atomname() != "C") continue ;
       //cout << "next current is " << next_curr+1 << endl;
       if (next_curr == temp_ring_vertices[0]) return true; //
       if (next_curr != temp_ring_vertices[temp_size-2]){ 
           bin2[0] = at_list1[0][next_curr].return_xbin();
           bin2[1] = at_list1[0][next_curr].return_ybin();
           bin2[2] = at_list1[0][next_curr].return_zbin();
       
           //Following condition will insure that tree search will be limited only
           //for those nodes that are in adjacent bins 
           flag = 0;
           j = 0;
           /*while (flag == 0 && j < 3 ){
                //if (bin1[j] == 0 && bin1[j]!= bin2[j] && bin2[j] < bin_dims[j]-2) return false;
                //if (bin1[j] == bin_dims[j]-1 && bin1[j]!= bin2[j] && bin2[j] > 1) return false;  
                //if (abs(bin1[j] - bin2[j]) > 2) return false;
                if (bin1[j] == 0 || bin1[j] == bin_dims[j]-1){
                    if (bin1[j] == 0 && bin1[j]!= bin2[j] && bin2[j] < bin_dims[j]-2){
                        flag = 1;
                    }
                    if (bin1[j] == bin_dims[j]-1 && bin1[j]!= bin2[j] && bin2[j] > 1){
                        flag = 1;  
                    }
                    j = j+1;
                }else if(abs(bin1[j] - bin2[j]) > 1) {
                    flag = 1;
                    j = j+1;
                }else{
                    j = j + 1;
                }
           }*/ 
          
          if (temp_ring_vertices.size() < max_ring_size-1){;
          //if (flag == 0){
              temp_ring_vertices.push_back(next_curr);
              if (check_inter_links(at_list1,temp_ring_vertices)) {
                  if(find_cycle(temp_ring_vertices,next_curr,at_list1,bin_dims,bin1)) return true;
              }else{
                  //cout << "returning here" << endl;
                  temp_ring_vertices.pop_back();
                  //return false;
              }
           }
        }        
   }  
   temp_ring_vertices.pop_back();
   return false;
}

bool check_inter_links(atom **at_list1, vector < int > const &temp_ring_vertices){
    int i,j,k,l;
    int jhigh,jhigh2;
    int ring_size;
    int vertex_index,neigh_index,neigh_index2;
    int flag;
    
    /*cout << "checking inter links" << endl;
    for (i=0;i<temp_ring_vertices.size();i++){
         cout << temp_ring_vertices[i]+1 << " ";
    } 
    cout << endl;*/
    ring_size = temp_ring_vertices.size();

    //Check interlinks for 1 vertex
    vertex_index = temp_ring_vertices[0];
    jhigh = at_list1[0][vertex_index].num_neighbors;
    for (j= 0; j<jhigh;j++){
         neigh_index = at_list1[0][vertex_index].neigh_indexes[j]-1; //neighbors[i][j];
         for (k=2;k<ring_size-1;k++){
              if (neigh_index == temp_ring_vertices[k]) return false;
         }
    
    }
    //Now do the same for other vertices
    for (i=1; i<ring_size-1;i++){
        vertex_index = temp_ring_vertices[i];
        jhigh = at_list1[0][vertex_index].num_neighbors;
        for (j= 0; j<jhigh;j++){
             neigh_index = at_list1[0][vertex_index].neigh_indexes[j]-1; //neighbors[i][j];
            for (k=i+2;k<ring_size;k++){
                 if (neigh_index == temp_ring_vertices[k]) return false;
            }
        }
    }

    return true;
}

bool check_rings(vector < vector < int > > const &rings, vector < int > const &temp_ring_vertices){
   int i,j,k;
  
   if (rings.size() == 0) return true; 
   for (i=0; i< rings.size(); i++){
       if (rings[i] == temp_ring_vertices){
          return false;
       }
   }
   return true;
}

bool check_intra_links(atom **at_list1, vector < int > const &temp_ring_vertices){

    int i,j,k,l;
    int jhigh,jhigh2;
    int ring_size;
    int vertex_index,neigh_index,neigh_index2;
    int flag;

    cout << "checking inter links" << endl;
    for (i=0;i<temp_ring_vertices.size();i++){
         cout << temp_ring_vertices[i]+1 << " ";
    } 
    cout << endl;
    ring_size = temp_ring_vertices.size();

    for (i=0;i<ring_size-1;i++){
         vertex_index = temp_ring_vertices[i];
         jhigh = at_list1[0][vertex_index].num_neighbors;
         for (j= 0; j<jhigh;j++){
              neigh_index = at_list1[0][vertex_index].neigh_indexes[j]-1; //neighbors[i][j];
              flag = 0;
              for (k=0;k<ring_size;k++){
                   if (neigh_index == temp_ring_vertices[k]){
                       flag = 1;
                       break;
                   }       
              }
              if (flag == 0){
                 //cout << neigh_index+1 << endl;
                 //this neighbor is not a member of the ring
                 jhigh2 = at_list1[0][neigh_index].num_neighbors;
                 for (l=0;l<jhigh2;l++){
                      neigh_index2 = at_list1[0][neigh_index].neigh_indexes[l]-1;
                      if (neigh_index2 != vertex_index){
                          for (k=0;k<ring_size;k++){
                               if (neigh_index2 == temp_ring_vertices[k]) {
                                   cout << "Found problem in 2nd search" << endl;
                                   exit(1);
                                   return false;
                               }
                          }
                      }
                 }
              }
         }
    }

    return true;
}


void write_ring_xyz(int act_frame1,int natoms1, vector < vector < int > > const &rings1, atom **at_list1,double **lx1, double **ly1, double **lz1, settings settings_list1){
   int i,j,k;
   int at_id,flag1;
   double temp_x,temp_y,temp_z;
   double x1,y1,z1,a,b,c,xc,yc,zc;
   double mag_v,vx,vy,vz;
   double xh,yh,zh;
   vector < double > xvector;
   vector < double > yvector;
   vector < double > zvector;
   vector < double > new_x;
   vector < double > new_y;
   vector < double > new_z;
   vector < string > new_name;
   //int at_ring_flags[natoms1][max_ring_size];
   int at_count_per_size[max_ring_size];
   int c2_flag[natoms1];
   int ring_size;
   ofstream ofile1; 
   string ring_names[max_ring_size];

   ring_names[0] = "X";
   ring_names[1] = "X";
   ring_names[2] = "X";
   ring_names[3] = "H";
   ring_names[4] = "Li";
   ring_names[5] = "Be";
   ring_names[6] = "C";
   ring_names[7] = "N";
   ring_names[8] = "O";
   ring_names[9] = "F";


   if (act_frame1 == 0){
      ofile1.open("rings.xyz");
   }else{
      ofile1.open("rings.xyz",ios::app);
   }

  for (i=0;i<natoms1;i++){
       at_list1[0][i].set_ring_defect_flag(0);
       /*for (j=0;j<max_ring_size;j++){
            at_ring_flags[i][j] = 0;
       }*/
  }

  for (i=0;i<max_ring_size;i++){
       at_count_per_size[i] = 0;
  }
  
   //write first two header lines of xyz frame
   ofile1 << rings1.size() << endl;
   ofile1 << endl;
   for (i=0;i<rings1.size();i++){
        temp_x = 0;
        temp_y = 0;
        temp_z = 0;

        for (j=0;j<rings1[i].size();j++){
             at_id = rings1[i][j];
             ring_size = rings1[i].size();
             if (ring_size < 3 || ring_size > max_ring_size){
                 cout << "Ring size exceeding size bounds " << endl;
                 exit(1);
             }
             /*if (at_ring_flags[at_id][ring_size] == 0) {
                 at_ring_flags[at_id][ring_size] = 1;
             }*/
             if (ring_size != ideal_ring_size){
                 if (at_list1[0][at_id].return_ring_defect_flag() == 0){
                     at_list1[0][at_id].set_ring_defect_flag(1);
                 }
             }
             at_count_per_size[ring_size] = at_count_per_size[ring_size] + 1;
             x1 = at_list1[0][at_id].return_xcord();
             y1 = at_list1[0][at_id].return_ycord();
             z1 = at_list1[0][at_id].return_zcord();

             xvector.push_back(x1);
             yvector.push_back(y1);
             zvector.push_back(z1);
        }
        //check for periodic boundary conditions
        if ((*max_element(xvector.begin(),xvector.end())-*min_element(xvector.begin(),xvector.end())) > lx1[0][2]/2.0){
            //cout << "Found periodicity for " << i+1 << " ring in x direction" << endl;
            for (j=0;j<xvector.size();j++){
                 if (xvector[j] < (lx1[0][0]+lx1[0][2]/2.0)) xvector[j] = xvector[j] + lx1[0][2];
            }
        }
        if ((*max_element(yvector.begin(),yvector.end())-*min_element(yvector.begin(),yvector.end())) > ly1[0][2]/2.0){
            for (j=0;j<yvector.size();j++){
                 if (yvector[j] < (ly1[0][0]+ly1[0][2]/2.0)) yvector[j] = yvector[j] + ly1[0][2];
            }
        }
        if ((*max_element(zvector.begin(),zvector.end())-*min_element(zvector.begin(),zvector.end())) > lz1[0][2]/2.0){
            for (j=0;j<zvector.size();j++){
                 if (zvector[j] < (lz1[0][0]+lz1[0][2]/2.0)) zvector[j] = zvector[j] + lz1[0][2];
            }
        }
        
        //Now calculate ring center
        for(j=0;j<xvector.size();j++){
             //cout << xvector[j] << " " ;
             temp_x = temp_x + xvector[j];
             temp_y = temp_y + yvector[j];
             temp_z = temp_z + zvector[j];
        }
        //cout << endl;
        //ofile1 << "C" << " ";
        ofile1 << ring_names[rings1[i].size()] << " ";
        ofile1 << temp_x/xvector.size() << " " ; 
        ofile1 << temp_y/yvector.size() << " " ; 
        ofile1 << temp_z/zvector.size() << " " ; 
        ofile1 << endl;
        
        xvector.clear();
        yvector.clear();
        zvector.clear();
   }

   ofile1.close();

   ofile1.open("6_ring_no_ring.xyz");
   ofile1 << natoms1 << endl;
   ofile1 << endl;
   for (i=0;i<natoms1;i++){
       flag1 = 0;
       if (at_list1[0][i].return_atomname() == "C"){
           /*for (j=0;j<max_ring_size;j++){
                if (j==6) continue;
                if (at_ring_flags[i][j] == 1){
                   flag1 = 1;
                   break;
                }
           }*/
           if (at_list1[0][i].return_ring_defect_flag() == 0){
               ofile1 << "C" << " " << at_list1[0][i].return_xcord() << " " << at_list1[0][i].return_ycord() << " " << at_list1[0][i].return_zcord() << endl;
           }else{
               ofile1 << "O" << " " << at_list1[0][i].return_xcord() << " " << at_list1[0][i].return_ycord() << " " << at_list1[0][i].return_zcord() << endl;
          }
       }else{
               ofile1 << at_list1[0][i].return_atomname() << " " << at_list1[0][i].return_xcord() << " " << at_list1[0][i].return_ycord() << " " << at_list1[0][i].return_zcord() << endl;
       }
   }
   ofile1.close();

  //write coordinates of all atoms that are part of particular ring size
  /*for(i=3;i<max_ring_size;i++){
      std::string filename = std::to_string(i) + "_atoms.xyz";
      ofile1.open(filename);
      ofile1 << at_count_per_size[i] << endl;
      ofile1 << endl;
      for (j=0;j<natoms1;j++){
           if (at_ring_flags[j][i] == 1){
               ofile1 << "C" << " " << at_list1[0][j].return_xcord() << " " << at_list1[0][j].return_ycord() << " " << at_list1[0][j].return_zcord() << endl;
           }
      }
      ofile1.close();
  }*/

}

void analyze_ring_size(int act_frame1, vector< vector < int > > const &rings){
     int i,j,k;
     int size;
     int ring_pop[max_ring_size];    
     ofstream ofile1;


     if (act_frame1 == 0){
         ofile1.open("rings_size_distribution.dat");
         for (i=0;i<max_ring_size;i++){
              ofile1 << i << "_ring " ;
         }
         ofile1 << endl;
     }else{
         ofile1.open("rings_size_distribution.dat",ios::app);
     }

     for (i=0;i<max_ring_size;i++){
          ring_pop[i] = 0;
     }

     for (i=0;i<rings.size();i++){
          size = rings[i].size();
          ring_pop[size] = ring_pop[size] + 1;
     }
 
     for (i=0;i<max_ring_size;i++){
          ofile1 << ring_pop[i] << " ";
     }
     ofile1 << endl;
     ofile1.close();

}

void check_periodic(double &n1,double &n2,double low,double len){
     if (abs(n1-n2) > len/2.0){
         if (n2 > low+len/2.0){
             n2 = n2 - len;
         }else{
             n2 = n2 + len;
         }
     }
}
