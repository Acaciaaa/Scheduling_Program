#include <vector>
#include <string>
using namespace std;
#define NeedSchedule(a) (!(a == 1 || a == 2 || a == 5))

typedef struct{
    bool isSequential;
    int area;
    int latency;
    float delay;
    bool isPipelined;
    int nr;
    vector<int> optypes; // 这个resource可以执行的运算类型：nr个

    int num = 0; // 这个resource的数量
}ResourceInf;

typedef struct{
    int nop, npred, nsucc;
    float times;
    vector<int> ops; // 这个block内的运算：nop个
    vector<int> pred, succ;

    int end = 0; // 代表下一个block开始的cycle
}BlockInf;

typedef struct{
    int optype;
    int ninputs;
    vector<int> input;

    int cycle; // 非相对于block，是绝对的
    int instance = -1; // 实例
    int cycle_h;
    int instance_h = -1;
}OpInf;

extern int NResourceType, NOpType, AreaLimit;
extern float Target;
extern vector<ResourceInf> ResourceLibrary;
extern int NBlock, NOperation;
extern vector<int> OpCategory;
extern vector<BlockInf> Block;
extern vector<OpInf> Op;

extern int* OpTool;
extern bool Test;

typedef struct{
    int nsucc = 0;
    vector<int> succ;
    int npred = 0;
    vector<int> pred;
}tmpInf;

typedef struct{
    int num;
    vector<int> borrow;
}record;