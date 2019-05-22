#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include "atom.h"
#include "functions.h"
#include "settings.h"

using namespace std;

void get_settings(settings& settings_list1){
     int i,j,k;
     int count;
     string line,sub;
     string temp[50];
     ifstream myfile;
     istringstream iss;

     myfile.open("settings.txt");
     if(myfile){
        while (getline(myfile,line)){
               iss.clear();
               iss.str(line);
               count = 0;
               while (iss >> sub){
                      temp[count] = sub;
                      count = count + 1;
                      if (count > 50){
                          cout << "Cannot have more than 50 lines in settings file" << endl;
                          exit(1);
                      }
               }
               if (temp[0] == "input_file"){
                   settings_list1.struct_name = temp[1];
               }
               if (temp[0] == "table_flag"){
                   settings_list1.table_flag = atoi(temp[1].c_str());
               }
               if (temp[0] == "ring_flag"){
                   settings_list1.ring_flag = atoi(temp[1].c_str());
               }
               if (temp[0] == "void_flag"){
                   settings_list1.void_flag = atoi(temp[1].c_str());
               }
               if (temp[0] == "periodic_flag"){
                   settings_list1.periodic_flag = atoi(temp[1].c_str());
               }
               if (temp[0] == "grid_size"){
                   settings_list1.periodic_flag = atof(temp[1].c_str());
               }
               if (temp[0] == "bo_cut_off"){
                   settings_list1.bo_cut_off = atof(temp[1].c_str());
               }
               if (temp[0] == "distance_cut_off"){
                   settings_list1.distance_cut_off = atof(temp[1].c_str());
               }
               if (temp[0] == "empty_bin_cut_off"){
                   settings_list1.empty_bin_cut_off = atoi(temp[1].c_str());
               }
        }
     }else{
        cout << "Settings file not found" << endl;
        exit(1);
     }
     myfile.close();
}
