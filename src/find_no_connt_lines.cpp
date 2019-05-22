#include <iostream>
#include <fstream>
#include "functions.h"

using namespace std;

int find_no_connt_lines(istream &myfile){
     using namespace std;
     
     int i,j,k;
     int a,nlines1;
     string line;
     string filename;
     istringstream instream;
     
     cout << "Counting number of lines in file" << endl;
     
     nlines1 =0;
     //myfile.open("bonds.txt");
     if(myfile){
                while (getline(myfile, line))
                      {
                       nlines1 = nlines1 + 1;
                       }
                 }else{
                       cout << "connection table file not found" << endl;
                       //system("Pause");
                       //exit(1);
                       }
     //myfile.clear();
     //myfile.seekg(0, ios::beg);
     //myfile.close();
     
     cout <<"Number of lines in connection table are " << nlines1 << endl;
     
     return nlines1;
}
