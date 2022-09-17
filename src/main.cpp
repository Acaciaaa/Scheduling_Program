#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include "my.h"
using namespace std;

extern void read(char* name);
extern void read_();
extern void cal();
extern void write(char* name);
extern void write_();

int main(int argc, char* argv[]){
    read(argv[1]);
    cal();
    if(argc == 3)
        write(argv[2]);
    else
        write_();
   
    delete[] OpTool;
    return 0;
}