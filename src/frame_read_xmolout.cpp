#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "atom.h"
#include "functions.h"
#include "settings.h"

using namespace std;
//********************************************************************************************
void find_no_lines(int &nlines1,int &natoms1,int &nframes1, string filename){
     using namespace std;
     
     int i,j,k;
     int a;
     string line;
     ifstream myfile;
     istringstream instream;

     //cout << "Counting number of lines in file" << endl;
     
     nlines1 =0;
     myfile.open(filename);
     if(myfile){
        while (getline(myfile, line))
            {
             nlines1 = nlines1 + 1;
            }
     }else{
      cout << "Input structure file not found" << endl;
      exit(1);
     }
     myfile.clear();
     myfile.seekg(0, ios::beg);
     myfile.close();
     
     cout <<"Number of lines in file are " << nlines1 << endl;
     
    myfile.open (filename);
    getline (myfile,line);
    instream.clear();
    instream.str(line);
     if (instream >> a ) {
            instream >> ws;        // Skip white space, if any.
            if (instream.eof()) {  // true if we're at end of string.
                //cout << "OK." << endl;
            } else {
                cout << "BAD. Too much on the line." << endl;
                //system("Pause");
                //exit(1);
            }
        } else {
            cout << "BAD: Didn't find the three items." << endl;
            //system("Pause");
            //exit(1);
        }
    natoms1 = a;
    myfile.clear();
    myfile.seekg(0, ios::beg);
    myfile.close();
    
    
    nframes1 = nlines1/(natoms1 + 2);
    cout << "number of atoms is " << natoms1 <<endl;
    cout << "number of frames is " << nframes1 << endl;

}
//********end of find no of lines function ************************************************

//************begining of read input xmolout function**************************************
int read_input_xmolout(istream &myfile, int nframes1, int natoms1, int periodic, atom **at_list1, string *cell_parameter_lines1,
                         double **lx1, double **ly1, double **lz1, double *alpha1, double *beta1, double *gama1){
    int i,j,k;
    int a,count;
    double xc,yc,zc;
    double xlow,xhigh,ylow,yhigh,zlow,zhigh;
    double dx,dy,dz;
    double ax,ay,az,ang1,ang2,ang3;
    //string ans;
    string line,line1,at,sub; 
    //ifstream myfile;
    istringstream instream,iss;
    
          
    // read in the xmolout file
    
    cout << "reading xyz file " << endl;
    for (i=0;i<nframes1;i++){
    //i=0;
          getline (myfile,line);
          instream.clear();
          instream.str(line);
          if (instream >> a ) {
             instream >> ws;        // Skip white space, if any.
             if (instream.eof()) {  // true if we're at end of string.
                //cout << "OK." << endl;
                } else {
                  cout << "BAD. Too much on the line." << endl;
                  }
                  } else {
                  //cout << "BAD: Didn't find the required number of items." << endl;
                  }
          //cout << "no of atoms is in frame" << i << "is " << a << endl;
          getline (myfile,line);
          
          //cell_parameter_lines1[i] = line;
          //cout << "cell parameter lines " << line << endl;
    
    
          if(periodic == 0){
                 lx1[i][0] = -100.0;
                 lx1[i][1] = 100.0;
                 lx1[i][2] = 200.0;

                 ly1[i][0] = -100.0;
                 ly1[i][1] = 100.0;
                 ly1[i][2] = 200.0;

                 lz1[i][0] = -100.0;
                 lz1[i][1] = 100.0;
                 lz1[i][2] = 200.0;
               
                 alpha1[i] = 90.0;
                 beta1[i]  = 90.0;
                 gama1[i]  = 90.0;
          }else{
                find_cell_parameters(line,i,lx1,ly1,lz1,alpha1,beta1,gama1);
          }
          cout << lx1[0][0] << " " << ly1[0][0] << " " << lz1[0][0] << endl;
          //in the following part, if value of any coordinate is negative,
          //the the cell parameter of corresponding direction is added
          //so that all coordinates are positive.
          xlow  = lx1[0][0];
          xhigh = lx1[0][1];
          ylow  = ly1[0][0];
          yhigh = ly1[0][1];
          zlow  = lz1[0][0];
          zhigh = lz1[0][1];
          for (j=0;j<natoms1;j++){
                                line.clear();
                                getline (myfile,line);
                                iss.clear();
                                iss.str(line);
                                count =0;
                                //cout << "starting atom " << j << endl;
                                while(iss >> sub){
                                          if(count == 0){
                                             at = sub;
                                          }
                                          if(count == 1){
                                             xc = atof(sub.c_str());
                                             if (xc < xlow){
                                                 cout << "Warning: Atom is out of cell bounds in x-direction. Resetting cell bounds" << endl;
                                                 xlow = xc;
                                                 lx1[0][0] = xc;
                                             }
                                             if (xc > xhigh){
                                                 cout << "Warning: Atom is out of cell bounds in x-direction. Resetting cell bounds" << endl;
                                                 xhigh = xc;
                                                 lx1[0][1] = xc;
                                             }
                                             /*if (xc < lx1[0][0]){
                                                 xc = xc + lx1[0][2];
                                                 if (xc > lx1[0][1]){
                                                     cout << "Atom " << j+1 << " out of box bounds in x direction" << endl;
                                                     exit(1);
                                                 }*/
                                          }
                                          if(count == 2){
                                             yc = atof(sub.c_str());
                                             if (yc < ylow){
                                                 cout << "Warning: Atom is out of cell bounds in y-direction. Resetting cell bounds" << endl;
                                                 ylow = yc;
                                                 ly1[0][0] = yc;
                                             }
                                             if (yc > yhigh){
                                                 cout << "Warning: Atom is out of cell bounds in y-direction. Resetting cell bounds" << endl;
                                                 yhigh = yc;
                                                 ly1[0][1] = yc;
                                             }
                                          }
                                          if(count == 3){
                                             zc = atof(sub.c_str());
                                             if (zc < zlow){
                                                 cout << "Warning: Atom is out of cell bounds in z-direction. Resetting cell bounds" << endl;
                                                 zlow = zc;
                                                 lz1[0][0] = zc;
                                             }
                                             if (zc > zhigh){
                                                 cout << "Warning: Atom is out of cell bounds in y-direction. Resetting cell bounds" << endl;
                                                 zhigh = zc;
                                                 lz1[0][1] = zc;
                                             }
                                          }
                                          count = count + 1;
                                          }
                                   //cout << "line is "<< j << " " << xc <<" "<< yc <<" "<< zc <<" " << at << endl;
                                   at_list1[i][j].set_coordinates (xc,yc,zc,at);   
                                   }
        }  
        lx1[0][2] = lx1[0][1] - lx1[0][0];  
        ly1[0][2] = ly1[0][1] - ly1[0][0];  
        lz1[0][2] = lz1[0][1] - lz1[0][0];  
    //cout << at_list1[0][natoms1-1].return_atomname()<<endl;
    
    //myfile.close();
      cout << "finished reading xyz file " << endl;
     }
//*****************end of read input xmolout function*************************************

//***********begining of find cell parameters function************************************
void find_cell_parameters(string s1, int frame_no, double **lx1, double **ly1, double **lz1, double *alpha1, double *beta1, double *gama1){
     int i,j,k;
     int count;
     string temp1;
     string sub;
     string *temp_strs;
     istringstream iss(s1);
     istringstream iss1(s1);
     
     cout << s1 << endl;
     count= 0;
     
     while(iss >> sub){
               count = count + 1;
               }

     if (count < 9){
         cout << "Error:Cell parameter line does not contain all parameters" << endl;
         cout << "Exiting the code" << endl;
         exit(1);
     }
     cout << "Count is " << count << endl;
     
     temp_strs = new string[count];
     
     i=0;
     while (iss1 >> sub){
           temp_strs[i] = sub;
           //cout << temp_strs[i] << endl;
           i= i +1;
           }
     
     gama1[frame_no] = atof(temp_strs[count-1].c_str());
     beta1[frame_no]  = atof(temp_strs[count-2].c_str());
     alpha1[frame_no]  = atof(temp_strs[count-3].c_str());
     
     lz1[frame_no][0] = atof(temp_strs[count-5].c_str()); //zlo
     lz1[frame_no][1] = atof(temp_strs[count-4].c_str()); //zhi
     ly1[frame_no][0] = atof(temp_strs[count-7].c_str()); //ylo
     ly1[frame_no][1] = atof(temp_strs[count-6].c_str()); //yhi
     lx1[frame_no][0] = atof(temp_strs[count-9].c_str()); //xlo
     lx1[frame_no][1] = atof(temp_strs[count-8].c_str()); //xhi

     lx1[frame_no][2] = lx1[frame_no][1] - lx1[frame_no][0]; //lx/a
     ly1[frame_no][2] = ly1[frame_no][1] - ly1[frame_no][0]; //ly/b
     lz1[frame_no][2] = lz1[frame_no][1] - lz1[frame_no][0]; //lz/c
     
     cout << lx1[frame_no][0] << " " << ly1[frame_no][0] << " " << lz1[frame_no][0] << endl ; 
     
     delete[] temp_strs;
     }
//**********end of find cell parameters function******************************************

