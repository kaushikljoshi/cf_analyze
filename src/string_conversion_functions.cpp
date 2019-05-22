#include <iostream>
#include <cstdlib>
#include <fstream>
#include "functions.h"

using namespace std;

int string_to_integer(char  *s1){
       char *end =0;
       int d=0;
       
       d = strtol(s1, &end, 10);
       
       if(*end != 0){
               cout << "Error in converting string " << s1<< endl;
               system("Pause");
               exit(1);
               }
       return d;
       }

double string_to_double(char  *s1){
       char *end =0;
       double d=0.0;
       
       d = strtod(s1, &end);
       
       if(*end != 0){
               cout << "Error in converting string " << s1<< endl;
               system("Pause");
               exit(1);
               }
       return d;
       }

string float_to_string(int a){
       std::string s;
       std::stringstream out;
       out << a;
       s = out.str();
       return s;
}
