#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "atom.h"
#include "functions.h"

using namespace std;

void input_connection_table(istream myfile1, int frame_no, int natoms1, int hlines, atom **at_list1){

     int i,j,k;
     int n1,n2,n3;
     int frames;
     int count1,count2;
     int neighbours[natoms1];
     int indexes[natoms1][20];
     string index_names[natoms1][20];
     string s1,s2,s3,sub;
     string temp[50];
     istringstream iss;
     ofstream myfile2;
     
     myfile2.open("temp_connection_table.txt");
     
     myfile2 << "Frame no " << frame_no << endl;
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
                                           n1 = atoi(temp[0].c_str()); //atom no+1
                                           n1 = n1-1;
                                           n2 = atoi(temp[2].c_str()); //no of neighbours
                                           neighbours[n1] = n2;
                                           count2 = 0;
                                           for(j=0;j<n2;j++){
                                                             n3 = atoi(temp[j+3].c_str());
                                                             s3 = at_list1[frames][n3-1].return_atomname();
                                                             
                                                             indexes[n1][count2] = n3;
                                                             index_names[n1][count2] = s3;
                                                             
                                                             count2 = count2 + 1;
                                                             }
                                           }//if for header lines ends here
                                   }//main for loop ends here

     for(i=0;i<natoms1;i++){
                            myfile2 << setw(10) << i << " " << setw(3) << neighbours[i] << " " ;
                            
                            n1 = neighbours[i];
                            for(j=0;j<n1;j++){
                                              myfile2<< setw(2) << index_names[i][j] << " " << setw(10) << indexes[i][j] << " " ;
                                              }
                            myfile2 << endl;
                            }
     myfile2.close();
}
