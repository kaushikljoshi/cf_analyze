#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include "atom.h"
#include "functions.h"

using namespace std;

void connection_table(istream myfile1, int frame_no, int natoms1, int hlines, atom **at_list1){

     int i,j,k;
     int n1,n2,n3;
     int frames;
     int count1,count2;
     //int neighbours[natoms1];
     //int indexes[natoms1][20];
     //string index_names[natoms1][20];
     string s1,s2,s3,sub;
     string temp[50];
     int file_write_flag = 0 // 1 means write myfile2
     istringstream iss;
     ofstream myfile2;
    
     if (file_write_flag == 1){ 
         myfile2.open("temp_connection_table.txt");
         myfile2 << "Frame no " << frame_no << endl;
     }
     count1 = 0;
     frames = 0;
     
     for(i=0;i<natoms1+hlines;i++){
                                   getline(myfile1,s2);
                                   if(i>hlines-1){
                                           iss.clear();
                                           iss.str(s2);
                                           count2 = 0;
                                           while(iss >> sub){
                                                     temp[count2] = sub;
                                                     count2 = count2 + 1;
                                                     }
                                           if(i==2931){
                                            for (j=0;j<count2;j++){
                                                cout << temp[count2]<< endl;
                                              }
                                           } 
                                           n1 = atoi(temp[0].c_str()); //atom no+1
                                           n1 = n1-1;
                                           n2 = atoi(temp[2].c_str()); //no of neighbours
                                           //neighbours[n1] = n2;
                                           if (n2 > max_neigh){
                                              cout << "Atom " << i+1 << " has more than " << max_neigh << " so increase max_neigh " << endl;
                                              exit(EXIT_FAILURE);
                                           }
                                           at_list1[frames][i].set_num_neighbors(n2);
                                           count2 = 0;
                                           for(j=0;j<n2;j++){
                                                             n3 = atoi(temp[j+3].c_str());
                                                             s3 = at_list1[frames][n3-1].return_atomname();
                                                             
                                                             //indexes[n1][count2] = n3;
                                                             //index_names[n1][count2] = s3;
                                                             at_list1[frames][i].neigh_indexes[j] = n3;
                                                             
                                                             count2 = count2 + 1;
                                                             }
                                           }//if for header lines ends here
                                   }//main for loop ends here
     if (file_write_flag == 1){
        for(i=0;i<natoms1;i++){
                            myfile2 << setw(10) << i << " " << setw(3) << at_list1[frames][i].num_neighbors << " " ;
                            
                            n1 = at_list1[frames][i].num_neighbors;
                            for(j=0;j<n1;j++){
                                              n2 = at_list1[frames][i].neigh_indexes[j]
                                              myfile2<< setw(2) << at_list1[frames][j].return_atomname() << " " << setw(10) << n2+1 << " " ;
                                              }
                            myfile2 << endl;
                            }
        myfile2.close();
     }
}
