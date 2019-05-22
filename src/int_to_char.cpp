#include <sstream>
#include <string>
#include <iostream>

using namespace std;
int main ( int argc, char **argv) {
    char *ch;
    char *temp1;
    char temp2;
    int i,k,num1,size,j;
    int *int_num;

    cout << "Enter length of number " << endl;
    cin >> size;
    cout << "Enter the number " << endl;
    cin >> num1;
    //cout << num1 << endl ;
    std::ostringstream oss;
    
    ch = new char[size];
    int_num = new int[size];
    oss << num1;
    strcpy(ch, oss.str().c_str());
    
    for(k=0;k<size;k++){
                        temp2 = ch[k];
                        temp1 = &temp2;
                        int_num[k] = atoi(temp1);
                      }
    
    for(k=0;k<size;k++){
                        cout << int_num[k] << endl;
                        }
    
    delete[] ch;
    delete[] int_num;
    system("Pause");
}
