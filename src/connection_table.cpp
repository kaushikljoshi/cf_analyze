#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <utility>
#include <chrono>
#include <math.h>
#include "atom.h"
#include "functions.h"
#include "settings.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{   
    int natoms,nframes,nlines;
    int i,j,k,l;
    int header_lines;
    atom **at_list;
    settings settings_list;
    //double **total_bond_order;
    int frames1,frames2;
    string * cell_parameter_lines;
    string *atom_name_types;
    //string periodic;
    double **lx,**ly,**lz;
    double *lows, *highs;
    double *alpha,*beta,*gama;
    int ****grid_list,***grid_atom_count;
    int grid_max_atoms = 100;                    //maximum 100 atoms allowed in each bin
    double grid_dr,act_grid_dr[3]; //cut-off distance for grid spacing in Ang. 
    double xremin,yremin,zremin;
    int nx,ny,nz;
    double remainder;
    int quo;
    int nx_old,ny_old,nz_old;
    ifstream myfile,myfile1,file;
    ofstream myfile2,myfile3,myfile4,myfile5,myfile6;
    
    int no_lines,ans,nlines1;
    int num1;
    string filename,line1,sub,line2,sub1;
    //string *lines;
    string *atom_types;
    istringstream iss,iss1;
    char *temp1;
    string ans1;
    
    string input_atom_types[20];

   settings_list = settings();
    
   get_settings(settings_list);
   cout << "Input structure name is " << settings_list.struct_name << endl;
   find_no_lines(nlines,natoms,nframes,settings_list.struct_name); //************************************
   //cout << nlines << "  " << natoms << "  " << nframes << endl;
    
   if(natoms == 0){
      cout << "No atoms in the system " << endl;
      exit(1);
   }
    
  frames2 = nframes;
  frames1 = 1;
    // allocate memory for atoms
    at_list = new atom *[frames1];
    //total_bond_order = new double *[frames1];
    if (at_list == NULL){
                cout <<"could not allocate memory \n";
                return 1;
                }
    /*if(total_bond_order == NULL){
                        cout << "could not allocate memory \n";
                        return 1;
                        }*/
    
    for (i=0;i<frames1;i++){
        at_list[i] = new atom[natoms];
        //total_bond_order[i] = new double[natoms];
    }
    
    //allocate memory for cell parameter lines
    cell_parameter_lines = new string[frames1];
    
    //allocate memory for cell parameter lines
    lx= new double *[frames1];
    ly= new double *[frames1];
    lz= new double *[frames1];
    for (i=0;i<frames1;i++){
         lx[i] = new double[3];
         ly[i] = new double[3];
         lz[i] = new double[3];
    }
    lows  = new double[3];
    highs = new double[3];
    
    alpha = new double[frames1];
    beta = new double[frames1];
    gama = new double[frames1];
    
    myfile1.open(settings_list.struct_name);
    myfile2.open("connection_table.txt");
    
    if(settings_list.table_flag == 1 || settings_list.table_flag == 2){
            if (settings_list.table_flag == 1){
               file.open("bonds.txt");
               nlines1 = find_no_connt_lines(file);
               file.close();
               //identify number of header lines in connection table
               file.open("bonds.txt");
               if(file){
                     header_lines = identify_header_lines(file);
                     //header_lines = 1;
                     //cout << "Number of header lines are " << header_lines << endl;
                     }else{
                             cout << "bonds.txt file not found" << endl;
                             system("Pause");
                             exit(1);
                             }
               file.clear();
               file.seekg(0, ios::beg);
               file.close();
               //frames2 = nlines1/(natoms + header_lines + 1);
               frames2 = nlines1/(natoms + header_lines + 1);
               cout << "no of frames " << frames2 << endl;
               //cout << "Enter bond order cut_off " << endl;
               //cin >> cut_off;
               file.open("bonds.txt");
             }else{
              file.open("fort.8");
              nlines1 = find_no_connt_lines(file);
              //file.clear();
              //file.seekg(0, ios::beg);
              file.close();
              frames2 = nlines1/(natoms + 1);
              //cout << "Enter bond order cut_off " << endl;
              //cin >> cut_off;
              file.open("fort.8");
             }
            read_input_xmolout(myfile1,frames1,natoms, settings_list.periodic_flag, at_list,cell_parameter_lines,
                              lx,ly,lz,alpha,beta,gama);
                              
            for(i=0;i<frames2;i++){
                if (settings_list.table_flag == 1){
                    read_table(file,i+1,natoms,header_lines,settings_list.bo_cut_off,at_list);
                }
                if (settings_list.table_flag == 2){
                    read_reax_table(file,i+1,natoms,settings_list.bo_cut_off,at_list);
                }
                
                cout << "frame " << i+1 << " is completed" << endl; 
            }
    }else{
	  grid_dr = settings_list.grid_size; 
          for(i=0;i<frames2;i++){
              read_input_xmolout(myfile1,frames1,natoms,settings_list.periodic_flag, at_list,cell_parameter_lines,
                              lx,ly,lz,alpha,beta,gama); //**************************
                              
              //cout << at_list[0][0].return_xcord() << endl; 
    
	      if (i == 0){
		  nx_old = 0;
		  ny_old = 0;
	          nz_old = 0;
	      }else{
		  nx_old = nx;
		  ny_old = ny;
		  nz_old = nz;
	      }
	      //calculate new nx,ny,nz.adjust grid_size in each direction so that each cell will have same size in each direction
	      //This is required for npt simulations because box size can change.
	      //Find nx
              nx = int(lx[frames1-1][2]/grid_dr);
	      xremin = std::fmod(lx[frames1-1][2],grid_dr);
              if (xremin != 0){
                  act_grid_dr[0] = grid_dr+xremin/nx;
              }else{
                 act_grid_dr[0] = grid_dr;
              }
              nx = lx[frames1-1][2]/act_grid_dr[0];
              //cout << "grid_dr after adjustment is " << act_grid_dr[0] << endl;
	      //xremin = std::fmod(lx[frames1-1][2],grid_dr);
              //cout << "new nx is " << nx << endl;
              //exit(1);
              //Find ny
              ny = int(ly[frames1-1][2]/grid_dr);
	      yremin = std::fmod(ly[frames1-1][2],grid_dr);
              if (yremin != 0){
                  act_grid_dr[1] = grid_dr+yremin/ny;
              }else{
                 act_grid_dr[1] = grid_dr;
              }
              ny = ly[frames1-1][2]/act_grid_dr[1];
              //find nz
              nz = int(lz[frames1-1][2]/grid_dr);
	      zremin = std::fmod(lz[frames1-1][2],grid_dr);
              if (zremin != 0){
                  act_grid_dr[2] = grid_dr+zremin/nz;
              }else{
                 act_grid_dr[2] = grid_dr;
              }
              nz = lz[frames1-1][2]/act_grid_dr[2];

              cout << lx[frames1-1][2] << " " << ly[frames1-1][2] << " " <<lz[frames1-1][2] << endl;
              cout << act_grid_dr[0] << " " << act_grid_dr[1] << " " << act_grid_dr[2] << endl;
              cout << nx << " " << ny << " " << nz << endl;
	      if (i == 0){
                  grid_list = new int ***[nx];
                  grid_atom_count = new int **[nx];
                  for (j=0;j<nx;j++){
                       grid_list[j] = new int **[ny];
                       grid_atom_count[j] = new int *[ny];
                       for (k=0;k<ny;k++){
                            grid_list[j][k] = new int *[nz];
                            grid_atom_count[j][k] = new int [nz];
                            for (l=0;l<nz;l++){
                                 grid_list[j][k][l] = new int [grid_max_atoms];
                                 grid_atom_count[j][k][l] = 0;
                            }
                       }
                   }
	      }else{  
                  cout << "starting memoery de-allocation" << endl;
                  //delete and re-assign grid array if change cell grid
		  if (nx != nx_old || ny != ny_old || nz != nz_old){
                      for (j=0;j<nx_old;j++){
                           for (k=0;k<ny_old;k++){
                                for (l=0;l<nz_old;l++){
                                     delete[] grid_list[j][k][l];
                                }
                                delete[] grid_list[j][k];
                                delete[] grid_atom_count[j][k];
                           }
                           delete[] grid_list[j];
                           delete[] grid_atom_count[j];
                      }
                      delete[] grid_list;
                      delete[] grid_atom_count;
                      cout << "Finished memory de-allocation" << endl;
                      grid_list = new int ***[nx];
                      grid_atom_count = new int **[nx];
                      for (j=0;j<nx;j++){
                           grid_list[j] = new int **[ny];
                           grid_atom_count[j] = new int *[ny];
                           for (k=0;k<ny;k++){
                                grid_list[j][k] = new int *[nz];
                                grid_atom_count[j][k] = new int [nz];
                                for (l=0;l<nz;l++){
                                     grid_list[j][k][l] = new int [grid_max_atoms];
                                     grid_atom_count[j][k][l] = 0;
                                }
                           }  
                      }
                  }else{
                      for (j=0;j<nx;j++){
                           for (k=0;k<ny;k++){
                                for(l=0;l<nz;l++){
                                    grid_atom_count[j][k][l] = 0;
                                } 
                           }
                      } 
                  }
	      }
	      bin_atoms(i,natoms,act_grid_dr,grid_max_atoms,at_list,lx,ly,lz,grid_list,grid_atom_count,nx,ny,nz,settings_list);
    
              make_connection_table(myfile2,i,natoms,at_list,lx,ly,lz,alpha,beta,gama,grid_list,grid_atom_count,nx,ny,nz,settings_list);
              
              if (settings_list.ring_flag == 1){
                  auto start_time = high_resolution_clock::now();
                  ring_analysis(i,natoms,at_list,nx,ny,nz,lx,ly,lz,settings_list);
                  auto end_time   = high_resolution_clock::now();
                  auto duration = duration_cast<milliseconds>(end_time-start_time);
                  cout << "Time taken by function: " << duration.count() << endl;
              }
                           
              cout << "frame " << i+1 << " is completed" << endl; 
              cout << "******************************************"<< endl;
	      if ( i == frames2-1){  //last frame. deallocate all memory related to grid lists
                   for (j=0;j<nx;j++){
                        for (k=0;k<ny;k++){
                             for (l=0;l<nz;l++){
                                  delete[] grid_list[j][k][l];
                             }
                            delete[] grid_list[j][k];
                            delete[] grid_atom_count[j][k];
                        }
                        delete[] grid_list[j];
                        delete[] grid_atom_count[j];
                   }
                   delete[] grid_list;
                   delete[] grid_atom_count;
	      }
         } // frame loop ends here
    }
    
    
    
//    cout << total_bond_order[0][4] << endl;
    myfile1.close();
    myfile2.close();
    file.close();
    //deallocate atom_list
    for(i=0;i<frames1;i++){
                           delete[] at_list[i];
                           //at_list[i] = NULL;
                           //delete[] total_bond_order[i];
                           //total_bond_order[i] = NULL;
                           delete[] lx[i];
                           delete[] ly[i];
                           delete[] lz[i];
                           }
    delete[] at_list;
    //delete[] total_bond_order;
    delete[] cell_parameter_lines;
    delete[] lx;
    delete[] ly;
    delete[] lz;
    delete[] lows;
    delete[] highs;
    delete[] alpha;
    delete[] beta;
    delete[] gama;
    
    //delete[] lines;
    
    //system("PAUSE");
    return EXIT_SUCCESS;
}
