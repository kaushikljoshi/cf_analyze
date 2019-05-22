#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "atom.h"
#include "functions.h"

using namespace std;

int identify_header_lines(istream &myfile1){
    
    int i,j,k;
    int count1;
    string s1,s2,s3;
    string sub;
    istringstream iss;
    
    s1 = "N";
    count1 = 0;
    while (s1 == "N"){
          getline(myfile1,s2);
          iss.clear();
          iss.str(s2);
          while(iss >> sub){
                    if(sub == "#"){
                           count1 = count1 + 1;
                           break;
                           }else{
                                 s1 = "Y";
                                 }
                    }
          }
    
    //cout << "Number of header lines are " << count1 << endl;
    
    return count1;
    
}
