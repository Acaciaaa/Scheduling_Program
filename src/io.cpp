#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "my.h"
using namespace std;

int NResourceType, NOpType, AreaLimit;
float Target;
vector<ResourceInf> ResourceLibrary; // NResourceType个
int NBlock, NOperation;
vector<int> OpCategory; // NOpType个
vector<BlockInf> Block; // NBlock个
vector<OpInf> Op; // NOperation个

void read(char* name){
    fstream file;
    file.open(name, ios::in|ios::binary);
    if(!file.is_open())
        cout << "open error" << endl;

    string line;

    // head
    getline(file, line);
    istringstream sin(line);
    sin >> NResourceType >> NOpType >> Target >> AreaLimit;
    
    // Resource Library
    for(int i = 0; i < NResourceType; ++i){
        ResourceInf tmp;
        getline(file, line); sin.clear(); sin.str(line);
        int is_seq, area, is_pipe = 0;
        sin >> is_seq >> area;
        if(is_seq == 1){
            tmp.isSequential = true;
            tmp.area = area;
            getline(file, line); sin.clear(); sin.str(line);
            sin >> tmp.latency >> tmp.delay >> is_pipe;
        }
        else{
            tmp.isSequential = false;
            tmp.area = area;
            tmp.latency = 0;
            getline(file, line); sin.clear(); sin.str(line);
            sin >> tmp.delay;
        }
        if(is_pipe == 1)
            tmp.isPipelined = true;
        else
            tmp.isPipelined = false;
        getline(file, line); sin.clear(); sin.str(line);
        sin >> tmp.nr;
        int optype;
        for(int j = 0; j < tmp.nr; ++j){
            sin >> optype;
            tmp.optypes.push_back(optype);
        }
        ResourceLibrary.push_back(tmp);
    }

    // CDFG
    getline(file, line); sin.clear(); sin.str(line);
    sin >> NBlock >> NOperation;

    // Category
    getline(file, line); sin.clear(); sin.str(line);
    int cate;
    for(int i = 0; i < NOpType; ++i){
        sin >> cate;
        OpCategory.push_back(cate);
    }

    // Block
    for(int i = 0; i < NBlock; ++i){
        BlockInf tmp;
        getline(file, line); sin.clear(); sin.str(line);
        sin >> tmp.nop >> tmp.npred >> tmp.nsucc >> tmp.times;
        int control_n[3] = {tmp.nop, tmp.npred, tmp.nsucc};
        int t;
        for(int j = 0; j < 3; ++j){
            getline(file, line);
            if(!line.empty()){
                sin.clear(); sin.str(line);
                for(int k = 0; k < control_n[j]; ++k){
                    sin >> t;
                    switch(j){
                        case 0: tmp.ops.push_back(t); break;
                        case 1: tmp.pred.push_back(t); break;
                        default: tmp.succ.push_back(t); break;
                    }
                }
            }
        }
        Block.push_back(tmp);
    }

    // Data flow
    for(int i = 0; i < NOperation; ++i){
        OpInf tmp;
        getline(file, line); sin.clear(); sin.str(line);
        sin >> tmp.optype;
        getline(file, line); sin.clear(); sin.str(line);
        sin >> tmp.ninputs;
        int str;
        for(int j = 0; j < tmp.ninputs; ++j){
            sin >> str;
            tmp.input.push_back(str);
        }
        Op.push_back(tmp);
    }

    file.close();
}

void write(char* name){
    fstream file;
    file.open(name, ios::out|ios::binary);
    if(!file.is_open())
        cout << "open output.txt error" << endl;
    
    for(int i = 0; i < NOperation; ++i){
        if(NeedSchedule(OpCategory[Op[i].optype]))
            file << Op[i].cycle << " ";
        else
            file << "-1 ";
    }
    file << "\n";
    for(int i = 0; i < NResourceType; ++i)
        file << ResourceLibrary[i].num << " ";
    file << "\n";
    for(int i = 0; i < NOperation; ++i){
        if(OpCategory[Op[i].optype] >= 6)
            file << OpTool[Op[i].optype] << " " << Op[i].instance << "\n";
        else 
            file << "-1\n";
    }

    file.close();
}

void write_(){
    for(int i = 0; i < NOperation; ++i){
        if(NeedSchedule(OpCategory[Op[i].optype]))
            cout << Op[i].cycle << " ";
        else
            cout << "-1 ";
    }
    cout << "\n";
    for(int i = 0; i < NResourceType; ++i)
        cout << ResourceLibrary[i].num << " ";
    cout << "\n";
    for(int i = 0; i < NOperation; ++i){
        if(OpCategory[Op[i].optype] >= 6)
            cout << OpTool[Op[i].optype] << " " << Op[i].instance << "\n";
        else 
            cout << "-1\n";
    }
}