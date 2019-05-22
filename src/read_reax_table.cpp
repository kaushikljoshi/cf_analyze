//This file reads reax fort.8 file
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include "atom.h"
#include "functions.h"

using namespace std;

void read_reax_table(istream& myfile,int frame_no, int natoms1,double cut_off, atom **at_list1){
     
     int i,j,k;
     int n1,n2,n3;
     int frames,nbonds;
     int count1,count2,pos,ncount;
     //int neighbours[natoms1];
     //int indexes[natoms1][20];
     double d1;
     //string index_names[natoms1][20];
     string s1,s2,s3,sub;
     string temp[50];
     istringstream iss;
     ofstream myfile2;
     int file_write_flag = 0 ; //1 means write files
     
     
     cout << "Reading connection table " << frame_no << endl;
     if(file_write_flag == 1){
        myfile2.open("temp_connection_table.txt");
        myfile2 << "Frame no " << frame_no << endl;
     }
     count1 = 0;
     frames = 0;
     nbonds = 0;
     cout << cut_off << endl;
     for(i=0;i<natoms1+1;i++){
        s2.clear();
        getline (myfile,s2);
        //cout << s2 << endl;
        if(i>0){
           iss.clear();
           iss.str(s2);
           count2 = 0;
           while (iss >> sub){
                 //cout << sub << endl;
                 temp[count2] = sub;
                 count2 = count2 + 1;
           }
           /*if (i== 2930){
                for (j=0;j<count2;j++){
                     cout << temp[j] << endl;
                }
           }*/
           n1 = atoi(temp[0].c_str()); //atom no+1
           n1 = n1-1;
           n2 = nbonds;  //no of neighbours
           //cout << n1 << "    " << n2 << endl;
           //neighbours[n1] = n2;
           ncount = 0;
           for(j=0;j<n2;j++){
               n3 = atoi(temp[j+2].c_str());
               d1 = atof(temp[j+3+nbonds].c_str());
               if (d1 >= cut_off){
                   s3 = at_list1[frames][n3-1].return_atomname();
                   
                   //indexes[n1][ncount] = n3;
                   //index_names[n1][ncount] = s3;
                   at_list1[frames][n1].neigh_indexes[ncount] = n3;
                   
                   //count2 = count2 + 1;
                   ncount = ncount + 1;
                   if (ncount > max_neigh){
                       cout << "Increase maximum number of neighbors " << endl;
                       exit(EXIT_FAILURE);
                   }
               }
           }
           //neighbours[n1] = ncount;
           at_list1[frames][n1].set_num_neighbors(ncount);
        }else{//First header line
           //cout << "starting header line" << endl;
           //cout << s2 << endl;
           iss.clear();
           iss.str(s2);
	   count2 = 0;
           pos = 0;
           //cout << s2 << endl;
           while(iss >> sub){
                 //cout << sub << endl;
                 temp[count2] = sub;
                 count2 = count2 + 1;
                 if (sub == "#Bonds:"){
                     pos = count2;
                 }
           }
           if (pos == 0){
               cout << "Could not find bonds section in fort.8" << endl;
               exit (EXIT_FAILURE);
           }
           nbonds = atoi(temp[pos].c_str());
        }//if for header lines ends here
    }//main for loop ends here
    

    if (file_write_flag == 1){
        for(i=0;i<natoms1;i++){
            n1 = at_list1[frames][i].num_neighbors; //neighbours[i];
            myfile2 << setw(10) << i+1 << " " << setw(3) << n1 << " ";
            for(j=0;j<n1;j++){
                n2 = at_list1[frames][i].neigh_indexes[j];
                s2 = at_list1[frames][n2-1].return_atomname();
                myfile2<< setw(2) << s2 << " " << setw(10) << n2 << " " ;
            }
            myfile2 << endl;
        }
        myfile2.close();
    }
    
}
