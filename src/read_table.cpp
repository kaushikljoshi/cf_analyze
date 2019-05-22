#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include "atom.h"
#include "functions.h"

using namespace std;

void read_table(istream& myfile,int frame_no, int natoms1, int hlines, double cut_off,atom **at_list1){
     
     cout << "Reading connection table " << frame_no << endl;
     int i,j,k;
     int n1,n2,n3;
     int frames;
     int count1,count2;
     //int neighbours[natoms1];
     //int indexes[natoms1][max_neigh];
     double dn1,dn2;
     //string index_names[natoms1][14];
     string s1,s2,s3,sub;
     string temp[50];
     istringstream iss;
     ofstream myfile2;
     int file_write_flag = 0; //1 means write files
     
     
     //myfile2.open("temp_connection_table.txt");
     
     //myfile2 << "Frame no " << frame_no << endl;
     count1 = 0;
     frames = 0;
     cout << "Bond order cut_off is " << cut_off << endl;
     for(i=0;i<natoms1+hlines;i++){
        s2.clear();
        getline (myfile,s2);
        if(i>hlines-1){
           iss.clear();
           iss.str(s2);
           count2 = 0;
           while(iss >> sub){
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
               n2 = atoi(temp[2].c_str()); //no of neighbours
               //cout << i << "  " << n1 << "    " << n2 << endl;
               //neighbours[n1] = n2;
               count1 = 0;
               for(j=0;j<n2;j++){
                   n3 = atoi(temp[j+3].c_str());
                   dn1 = atof(temp[j+4+n2].c_str());
                   if (dn1 > cut_off){
                       s3 = at_list1[frames][n3-1].return_atomname();
                       //indexes[n1][count1] = n3;
                       //index_names[n1][count1] = s3;
                       at_list1[frames][n1].neigh_indexes[count1] = n3;
                       count2 = count2 + 1;
                       count1 = count1 + 1;
                       if (count1 > max_neigh){
                           cout << "Increase maximum number of neighbors\n";
                           exit(EXIT_FAILURE);
                       }
                   }
               }
               //neighbours[n1] = count1;
               at_list1[frames][n1].set_num_neighbors(count1);
        }//if for header lines ends here
     }//main for loop ends here
     getline(myfile,s2);

     if (file_write_flag == 1){
         for(i=0;i<natoms1;i++){
             n1 = at_list1[frames][i].num_neighbors;
             myfile2 << setw(10) << i+1 << " " << setw(3) << n1 << " ";
             for(j=0;j<n1;j++){
                 n3 = at_list1[frames][i].neigh_indexes[j];
                 s3 = at_list1[frames][n3-1].return_atomname();
                 myfile2<< setw(2) << s3 << " " << setw(10) << n3 << " " ;
             }
             myfile2 << endl;
         }
     }

     //myfile2.close();
    
}
