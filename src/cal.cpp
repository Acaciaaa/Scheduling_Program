#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <iterator>
#include <map>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include "my.h"
extern "C"{
    #include "glpk.h"
}
using namespace std;

int lat(int op_resource){
    return ResourceLibrary[op_resource].latency + 1;
}

int* OpTool;
int MaxOp = 15;
vector<int> block_order;
bool** dominants;

void dominant(vector<int> &path){
    int nowid = path[path.size()-1];
    vector<bool> come_by(NBlock);
    for(int i = 0; i < NBlock; ++i) come_by[i] = false;
    for(int i = 0; i < path.size(); ++i)
        come_by[path[i]] = true;

    // cout << "path:" << endl;
    // for(int i = 0; i < path.size(); ++i)
    //     cout << path[i] << " ";
    // cout << endl;
    // for(int i = 0; i < NBlock; ++i)
    //     cout << come_by[i] << " ";
    // cout << endl;

    for(int i = 0; i < NBlock; ++i)
        if(dominants[nowid][i] == true && come_by[i] == false)
            dominants[nowid][i] = false;
    for(int i = 0; i < Block[nowid].nsucc; ++i){
        int nextid = Block[nowid].succ[i];
        if(find(path.begin(), path.end(), nextid) == path.end()){
            path.push_back(nextid);
            dominant(path);
            path.pop_back();
        }
    }
}

// 处理基本块的顺序
void def_order(){
    dominants = new bool*[NBlock];
    for(int i = 0; i < NBlock; ++i){
        dominants[i] = new bool[NBlock];
        for(int j = 0; j < NBlock; ++j)
            dominants[i][j] = true;
    }
    if(Block[0].npred != 0)
        cout << "entry is not block0" << endl;
    else{
        vector<int> path = {0};
        dominant(path);
    }
    // cout << "dominants:" << endl;
    // for(int i = 0; i < NBlock; ++i){
    //     for(int j = 0; j <  NBlock; ++j)
    //         cout << dominants[i][j] << " ";
    //     cout << endl;
    // }
    // cout << "original:" << endl;
    // for(int i = 0; i < NBlock; ++i){
    //     cout << Block[i].npred << ":";
    //     for(int j = 0; j < Block[i].npred; ++j)
    //         cout << Block[i].pred[j] << " ";
    //     cout << endl << Block[i].nsucc << ":";
    //     for(int j = 0; j < Block[i].nsucc; ++j)
    //         cout << Block[i].succ[j] << " ";
    //     cout << endl;
    // }
    
    // 消回边
    vector<int> new_pred;
    for(int i = 0; i < NBlock; ++i){
        new_pred.clear();
        int minus = 0;
        for(int j = 0; j < Block[i].npred; ++j){
            int blockid = Block[i].pred[j];
            if(dominants[blockid][i]){ // 是回边
                //cout << blockid << " to " << i << endl;
                minus++;
                Block[blockid].nsucc--;
                for(auto p = Block[blockid].succ.begin(); p != Block[blockid].succ.end(); ++p)
                    if(*p == i){
                        Block[blockid].succ.erase(p);
                        break;
                    }
            }
            else
                new_pred.push_back(blockid);
        }
        Block[i].npred -= minus;
        Block[i].pred.assign(new_pred.begin(), new_pred.end());
    }

    // 拓扑排序
    int* degree = new int[NBlock];
    bool* if_order = new bool[NBlock];
    for(int i = 0; i < NBlock; ++i){
        if_order[i] = false;
        degree[i] = Block[i].npred;
    }

    while(block_order.size() != NBlock){
        for(int i = 0; i < NBlock; ++i)
            if(degree[i] == 0 && !if_order[i]){
                if_order[i] = true;
                block_order.push_back(i);
                for(int j = 0; j < Block[i].nsucc; ++j)
                    degree[Block[i].succ[j]]--;
            }
    }

    for(int i = 0; i < NBlock; ++i)
        delete[] dominants[i];
    delete[] dominants;
    delete[] if_order;
    delete[] degree;

    // for(int i = 0; i < NBlock; ++i)
    //     cout << block_order[i] << " ";
    // cout << endl;
}

vector<int> begin_time;
bool comp(const int & a, const int & b){
    return begin_time[a] < begin_time[b];
}

int heuristic(vector<int> & schedule_ops, int block_begin){
    vector<int> ready_list;

    // ALAP + 顺便初始化ready_list
    int op_num = schedule_ops.size();
    vector<tmpInf> node(op_num);
    for(int i = 0; i < op_num; ++i){ // 初始化node
        int myop = schedule_ops[i];
        bool ready = true;
        for(int j = 0; j < Op[myop].ninputs; ++j){
            int predop = Op[myop].input[j];
            if(predop == -1)
                continue;
            auto p = find(schedule_ops.begin(), schedule_ops.end(), predop);
            if(p != schedule_ops.end()){
                ready = false;
                int opid = distance(schedule_ops.begin(), p);
                node[opid].nsucc++;
                node[opid].succ.push_back(i);
                node[i].npred++;
                node[i].pred.push_back(opid);
            }
        }
        if(ready)
            ready_list.push_back(i); // 无前驱点的op
    }
    // for(int i = 0; i < op_num; ++i){
    //     cout << i << " pred:";
    //     for(int j = 0; j < node[i].npred; ++j)
    //         cout << " " << node[i].pred[j];
    //     cout << endl;
    // }
    
    int max_t = 0;
    int t = 0;
    vector<int> end_time(op_num);
    begin_time.clear();
    for(int i = 0; i < op_num; ++i) {end_time[i] = 1; begin_time.push_back(0);}
    for(int i = 0; i < op_num; ++i) // 最后一个cycle的node
        if(node[i].nsucc == 0){
            end_time[i] = 0;
            begin_time[i] = -lat(OpTool[Op[schedule_ops[i]].optype]);
            max_t = min(max_t, begin_time[i]);
        }
    bool flag = true;
    while(flag){
        t--;
        for(int i = 0; i < op_num; ++i)
            if(end_time[i] == 1){
                bool ready = true;
                for(int j = 0; j < node[i].nsucc; ++j)
                    if(end_time[node[i].succ[j]] == 1 || t > begin_time[node[i].succ[j]])
                        ready = false;
                if(ready){
                    end_time[i] = t;
                    begin_time[i] = t - lat(OpTool[Op[schedule_ops[i]].optype]);
                    max_t = min(max_t, begin_time[i]);
                }
            }
        flag = false;
        for(int i = 0; i < op_num; ++i)
            if(end_time[i] == 1){
                flag = true;
                break;
            }
    }
    max_t = -max_t;
    for(int i = 0; i < op_num; ++i)
        begin_time[i] += max_t;

    // cout << endl << "新block" << endl;
    // cout << "max_t:" << max_t << endl;
    // for(int i = 0; i < op_num; ++i)
    //     cout << i << ":" << begin_time[i] << " " << end_time[i] << endl;

    // List Scheduling
    // 工具卡片---------------------------------------------------------------------------
    map<int, record> m;
    for(int i = 0; i < NResourceType; ++i){
        int num = ResourceLibrary[i].num;
        if(num != 0){
            record tmp;
            tmp.num = num;
            for(int i = 0; i < num; ++i) tmp.borrow.push_back(-1);
            m.insert(make_pair(i, tmp));
        }
    }
    // 工具卡片---------------------------------------------------------------------------
    int max_real = 0;
    t = 0;
    vector<int> begin_real(op_num), end_real(op_num), instan(op_num); for(int i = 0; i < op_num; ++i) {begin_real[i] = -1; instan[i] = -1;}
    flag = true;
    while(flag){
        // 更新一下工具卡片：先不考虑紧接着执行的情况
        for(auto p = m.begin(); p != m.end(); ++p){
            int tool = p->first;
            if(ResourceLibrary[tool].isPipelined || lat(tool) == 1){
                for(int i = 0; i < p->second.num; ++i)
                    p->second.borrow[i] = -1;
            }
            else{
                for(int i = 0; i < p->second.num; ++i){
                    int borrower = p->second.borrow[i];
                    if(borrower == -1)
                        continue;
                    if(t >= end_real[borrower])
                        p->second.borrow[i] = -1;
                }
            }
        }
        //cout << "更新完工具了" << endl;
        sort(ready_list.begin(), ready_list.end(), comp);
        // cout << "ready_list:" << endl;
        // for(int i = 0; i < ready_list.size(); ++i)
        //     cout << ready_list[i] << " ";
        // cout << endl;

        // 调度ready_list
        for(int i = 0; i < ready_list.size(); ++i){
            int schedule_id = ready_list[i];
            int tool = OpTool[Op[schedule_ops[schedule_id]].optype];
            // 如果是store/load就不管工具直接怼上去
            int cate = OpCategory[Op[schedule_ops[schedule_id]].optype];
            if(cate == 3 || cate == 4){
                begin_real[schedule_id] = t;
                end_real[schedule_id] = t + lat(tool);
                max_real = max(max_real, end_real[schedule_id]);
                continue;
            }
            if(m.find(tool) == m.end()) cout << "这个工具不存在" << endl;
            for(int j = 0; j < m[tool].num; ++j){
                if(m[tool].borrow[j] == -1){
                    m[tool].borrow[j] = schedule_id;
                    begin_real[schedule_id] = t;
                    end_real[schedule_id] = t + lat(tool);
                    max_real = max(max_real, end_real[schedule_id]);
                    instan[schedule_id] = j;
                    break;
                }
            }
        }

        t++;

        // 更新ready_list
        ready_list.clear();
        for(int i = 0; i < op_num; ++i){
            //cout << i << ":" << begin_real[i] << " " << end_real[i] << endl;
            if(begin_real[i] == -1){
                bool ready = true;
                for(int j = 0; j < node[i].npred; ++j)
                    if(begin_real[node[i].pred[j]] == -1 || t < end_real[node[i].pred[j]])
                        ready = false;
                if(ready)
                    ready_list.push_back(i);
            }
        }
        // cout << "后ready_list:" << endl;
        // for(int i = 0; i < ready_list.size(); ++i)
        //     cout << ready_list[i] << " ";
        // cout << endl;

        flag = false;
        for(int i = 0; i < op_num; ++i)
            if(begin_real[i] == -1){
                flag = true;
                break;
            }
        // if(t == 10)
        //     break;
    }
    // cout << "最后：" << max_real << endl;
    // for(int i = 0; i < op_num; ++i)
    //     cout << i << ":" << begin_real[i] << " " << end_real[i] << endl;
    // cout << endl << endl;
    for(int i = 0; i < op_num; ++i){
        int opid = schedule_ops[i];
        Op[opid].cycle_h = block_begin + begin_real[i];
        Op[opid].instance_h = instan[i];
    }

    return max_real + block_begin;
}

int inner_schedule(vector<int> & raw_ops, int split, int begin_cycle){
    int max_latency = 0;
    vector<int> schedule_ops;
    for(int i = split; i < raw_ops.size() && i < (split + MaxOp); ++i){
        schedule_ops.push_back(raw_ops[i]);
        max_latency += lat(OpTool[Op[raw_ops[i]].optype]);
    }
    //cout << "!!!!!!max_latency:" << max_latency << endl;
    int op_num = schedule_ops.size();
    
    // if(schedule_ops[0] == 17){
    // cout << "更小的块:" << endl;
    // for(int i = 0; i < op_num; ++i){
    //     cout << schedule_ops[i] << ":";
    //     for(int j = 0; j < Op[schedule_ops[i]].ninputs; ++j)
    //         cout << " " << Op[schedule_ops[i]].input[j];
    //     cout << endl;
    // }
    // }

    int m = op_num * max_latency + 1;
    //cout << "m:" << m << endl;

    /*
    有max_latency * op_num + 1个原变量
    辅助变量[1, op_num]是∑t(x_it)
    */  

    // 初始化
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_cols(lp, m);
    // x_it是binary + m>=0 && 是整数
    for(int i = 1; i <= op_num; ++i)
        for(int j = 1; j <= max_latency; ++j){
            glp_set_col_kind(lp, (i - 1) * max_latency + j, GLP_BV);
            glp_set_col_bnds(lp, (i - 1) * max_latency + j, GLP_DB, 0, 1);
        }
    glp_set_col_bnds(lp, m, GLP_DB, 0.0, max_latency+1);
    glp_set_col_kind(lp, m, GLP_IV);
    // min m
    glp_set_obj_coef(lp, m, 1.0);
    // 辅助变量处理：辅助变量数*原变量数的二维数组
    vector<vector<int> > auxi;

    // ∑(x_it) = 1
    for(int i = 1; i <= op_num; ++i){
        // 声明辅助变量和界限
        glp_add_rows(lp, 1);
        glp_set_row_bnds(lp, i, GLP_FX, 1.0, 1.0);

        // 声明辅助变量和原变量关系
        vector<int> tmp(m);
        for(int j = 0; j < m; ++j)
            tmp[j] = 0;
        int beginning = (i - 1) * max_latency;
        for(int j = 0; j < max_latency; ++j)
            tmp[beginning + j] = 1;
        auxi.push_back(tmp);
    }
    
    int serial_number = op_num + 1; // 编号
    unordered_map<int, vector<int> > hash; // 存的是schedule_ops的下标
    
    // dependence constraint
    for(int i = 0; i < op_num; ++i){
        int myop = schedule_ops[i]; // 真实的op标号 可以用作Op下标
        // cout << "i:" << i << " op:" << myop << endl;
        int mytool = OpTool[Op[myop].optype];
        if(mytool == -1){
            cout << "mytool err: op:" << myop << endl;
            cout << "this optype:" << Op[myop].optype << " has no resource TAT" << endl;
            exit(1);
        }

        // 为resource constraint做准备
        if(OpCategory[Op[myop].optype] != 3 && OpCategory[Op[myop].optype] != 4){
            if(hash.find(mytool) == hash.end())
                hash.insert(make_pair(mytool, vector<int> {i}));
            else
                hash[mytool].push_back(i);
        }

        // m constraint---------------------------------------
        // 声明辅助变量和原变量关系
        vector<int> tmp(m);
        for(int j = 0; j < m; ++j)
            tmp[j] = 0;
        int beginning = i * max_latency;
        for(int j = 0; j < max_latency; ++j)
            tmp[beginning + j] = -(j + 1);
        tmp[m-1] = 1;
        auxi.push_back(tmp);

        // 声明辅助变量和界限
        glp_add_rows(lp, 1);
        glp_set_row_bnds(lp, serial_number, GLP_LO, lat(mytool), 0);

        // 辅助变量多一个：编号++
        ++serial_number;
        // m constraint---------------------------------------

        auto & input = Op[myop].input;
        for(int in = 0; in < input.size(); ++in){
            if(input[in] == -1)
                continue;
            // 这个input在我的原变量里排第几
            auto num_which = find(schedule_ops.begin(), schedule_ops.end(), input[in]);
            if(num_which == schedule_ops.end())
                continue;
            //cout << "distance:" << distance(num_which, schedule_ops.begin()) << endl;
            int tool = OpTool[Op[input[in]].optype];
            if(tool == -1){
                cout << "inputtool err: op:" << input[in] << endl;
                cout << "this optype:" << Op[input[in]].optype << " has no resource TAT" << endl;
                exit(1);
            }

            // 声明辅助变量和原变量关系
            vector<int> tmp(m);
            for(int j = 0; j < m; ++j)
                tmp[j] = 0;
            int beginning = i * max_latency;
            for(int j = 0; j < max_latency; ++j)
                tmp[beginning + j] = j + 1;
            for(int j = 0; j < max_latency; ++j)
                tmp[distance(schedule_ops.begin(), num_which) * max_latency + j] = -(j + 1);
            auxi.push_back(tmp);

            // 声明辅助变量和界限
            glp_add_rows(lp, 1);
            if(ResourceLibrary[tool].isSequential){
                if(ResourceLibrary[mytool].isSequential)
                    glp_set_row_bnds(lp, serial_number, GLP_LO, ResourceLibrary[tool].latency+1, 0);
                else{
                    if((ResourceLibrary[tool].delay + ResourceLibrary[mytool].delay) <= Target)
                        glp_set_row_bnds(lp, serial_number, GLP_LO, ResourceLibrary[tool].latency, 0);
                    else
                        glp_set_row_bnds(lp, serial_number, GLP_LO, ResourceLibrary[tool].latency+1, 0);
                }
            }
            else{
                glp_set_row_bnds(lp, serial_number, GLP_LO, 1, 0);
                // if(ResourceLibrary[mytool].isSequential)
                //     glp_set_row_bnds(lp, serial_number, GLP_LO, 1, 0);
                // else{
                //     if((ResourceLibrary[tool].delay + ResourceLibrary[mytool].delay) <= Target)
                //         glp_set_row_bnds(lp, serial_number, GLP_LO, 0, 0);
                //     else
                //         glp_set_row_bnds(lp, serial_number, GLP_LO, 1, 0);
                // }
            }

            // 辅助变量多一个：编号++
            ++serial_number;
        }
    }
    
    // if(schedule_ops[0] == 7){
    // for(int i = 0; i < auxi.size(); ++i){
    //     for(int j = 0; j < auxi[i].size(); ++j)
    //         cout << auxi[i][j] << " ";
    //     cout << endl;
    // }
    // cout << "serial_number: " << serial_number << endl;
    // for(int i = 1; i < serial_number; ++i)
    //     cout << "lb:" << glp_get_row_lb(lp, i) << " ub:" << glp_get_row_ub(lp, i) << endl;
    // }


    // resource constraint
    for(auto p = hash.begin(); p != hash.end(); ++p){
        ResourceInf & inf = ResourceLibrary[p->first];
        if(inf.isPipelined){
            for(int t = 0; t < max_latency; ++t){
                // 声明辅助变量和原变量关系
                vector<int> tmp(m);
                for(int j = 0; j < m; ++j)
                    tmp[j] = 0;
                for(auto q = p->second.begin(); q != p->second.end(); ++q)
                    tmp[(*q) * max_latency + t] = 1;
                auxi.push_back(tmp);

                // 声明辅助变量和界限
                glp_add_rows(lp, 1);
                glp_set_row_bnds(lp, serial_number, GLP_DB, 0.0, double(inf.num));

                // 辅助变量多一个：编号++
                ++serial_number;
            }
        }
        else{
            int tmp_lat = inf.latency;
            for(int t = 0; (t + tmp_lat) < max_latency; ++t){
                // 声明辅助变量和原变量关系
                vector<int> tmp(m);
                for(int j = 0; j < m; ++j)
                    tmp[j] = 0;
                for(auto q = p->second.begin(); q != p->second.end(); ++q){
                    for(int l = 0; l <= tmp_lat; ++l)
                        tmp[(*q) * max_latency + t + l] = 1;
                }
                auxi.push_back(tmp);

                // 声明辅助变量和界限
                glp_add_rows(lp, 1);
                glp_set_row_bnds(lp, serial_number, GLP_DB, 0.0, double(inf.num));

                // 辅助变量多一个：编号++
                ++serial_number;
            }
        }
    }
    
    // for(int i = 0; i < auxi.size(); ++i){
    //     for(int j = 0; j < auxi[i].size(); ++j)
    //         cout << auxi[i][j] << " ";
    //     cout << endl;
    // }
    // cout << "serial_number: " << serial_number << endl;
    // for(int i = 1; i < serial_number; ++i)
    //     cout << "lb:" << glp_get_row_lb(lp, i) << " ub:" << glp_get_row_ub(lp, i) << endl;

    // 加入原变量和辅助变量关系
    int total_num = (serial_number-1) * m;
    int *ia = new int[total_num + 1], *ja = new int[total_num + 1];
    double *ar = new double[total_num + 1];
    int n = 0;
    for(int i = 1; i < serial_number; ++i)
        for(int j = 1; j <= m; ++j){
            ++n;
            ia[n] = i;
            ja[n] = j;
            ar[n] = auxi[i-1][j-1];
            // cout << i << " " << j << " " << auxi[i-1][j-1] << endl;
        }
    glp_load_matrix(lp, total_num, ia, ja, ar);

    // 让glpk算
    glp_smcp p1;
    glp_init_smcp(&p1);
    p1.msg_lev = GLP_MSG_OFF;
    glp_simplex(lp, &p1);
    glp_iocp p2;
    glp_init_iocp(&p2);
    p2.msg_lev = GLP_MSG_OFF;
    glp_intopt(lp, &p2);
    // cout << "z:" << glp_mip_obj_val(lp) << endl;
    // for(int i = 0; i < op_num; ++i){
    //     for(int t = 0; t < max_latency; ++t)
    //         cout << glp_mip_col_val(lp, i * max_latency + t + 1) << " ";
    //     cout << endl;
    // }

    // 填上每个op的cycle和end
    int ans = begin_cycle + glp_mip_obj_val(lp);
    
    for(int i = 0; i < op_num; ++i)
        for(int t = 0; t < max_latency; ++t)
            if(glp_mip_col_val(lp, i * max_latency + t + 1) == 1){
                Op[schedule_ops[i]].cycle = begin_cycle + t;
                // if(schedule_ops[0] == 7){
                //     cout << schedule_ops[i] << ":" << t + 1 << endl;
                // }
            }
    // 填上每个op的instance 刚好用一下hash
    vector<int> free_instance(max_latency); // 给pipeline用的
    vector<vector<int> >op_start(max_latency); // 给非pipeline用的
    for(auto p = hash.begin(); p != hash.end(); ++p){
        ResourceInf & inf = ResourceLibrary[p->first];
        if(inf.isPipelined){ // 只有同一个cycle开始的会冲突
            for(int i = 0; i < max_latency; ++i) free_instance[i] = 0;
            for(auto q = p->second.begin(); q != p->second.end(); ++q){
                auto & tmpop = Op[schedule_ops[*q]];
                tmpop.instance = free_instance[tmpop.cycle-begin_cycle]++;
            }
        }
        else{ // 只要周期有重合就会冲突
            int interval = inf.latency + 1;
            for(int i = 0; i < max_latency; ++i) op_start[i] = vector<int>{};
            for(auto q = p->second.begin(); q != p->second.end(); ++q){
                int tmpp = schedule_ops[*q];
                op_start[Op[tmpp].cycle-begin_cycle].push_back(tmpp);
            }
            int ins = 0;
            while(ins < inf.num){
                int i = 0;
                while(i < max_latency){
                    if(op_start[i].size() > ins){
                        Op[op_start[i][ins]].instance = ins;
                        i += interval;
                    }
                    else
                        i++;
                }
                ins++;
            }
        }
    }
    
    // cout << "更小的块:" << endl;
    // for(int i = 0; i < op_num; ++i){
    //     int id = schedule_ops[i];
    //     cout << id << ":";
    //     cout << " cycle:" << Op[id].cycle << " tool&lat:" << OpTool[Op[id].optype] << " " << lat(OpTool[Op[id].optype]) << " ins:" << Op[id].instance;
    //     cout << endl;
    // }
    // cout << "ans:" << ans << endl << endl;
    

    delete[] ia, ja, ar;
    glp_delete_prob(lp);
    return ans;
}
// 调度
void schedule(int block){
    vector<int> schedule_ops;
    auto & ops = Block[block].ops;
    for(int i = 0; i < ops.size(); ++i){
        //cout << ops[i] << " ";
        int op_type = Op[ops[i]].optype;
        //cout << "op_type:" << op_type << endl;
        int op_category = OpCategory[op_type];
        //cout << "op_category:" << op_category << endl;
        if(NeedSchedule(op_category))
            schedule_ops.push_back(ops[i]);
    }
    int op_num = schedule_ops.size();

    int block_begin = 1;
    for(int i = 0; i < Block[block].npred; ++i)
        block_begin = max(block_begin, Block[Block[block].pred[i]].end);
    
    int split = 0, tmp = block_begin;
    while(split < op_num){
        tmp = inner_schedule(schedule_ops, split, tmp);
        split += MaxOp;
    }
    Block[block].end = tmp;

    int tmp_h = heuristic(schedule_ops, block_begin);
    //cout << "ilp:" << tmp << "  heu:" << tmp_h << endl;
    if(tmp_h < tmp){// 启发式的更好
        Block[block].end = tmp_h;
        for(int i = 0; i < op_num; ++i){
            int opid = schedule_ops[i];
            Op[opid].cycle = Op[opid].cycle_h;
            Op[opid].instance = Op[opid].instance_h;
        }
    }

    return;
}

// 遍历求最小覆盖
bool auxi_area(vector<vector<int> > &opTool, vector<bool> &opHaveTool, int &area_sum){
    if(area_sum > AreaLimit)
        return false;
    for(int i = 0; i < NOpType; ++i)
        if(opHaveTool[i] == false){
            for(int j = 0; j < opTool[i].size(); ++j){
            // 尝试用这个tool
            int tool = opTool[i][j]; auto p = &ResourceLibrary[tool];
            area_sum += p->area;
            vector<int> tmpop;
            for(int j = 0; j < p->nr; ++j){
                int tmp = p->optypes[j];
                if(opHaveTool[tmp] == false){
                    opHaveTool[tmp] = true;
                    OpTool[tmp] = tool;
                    tmpop.push_back(tmp);
                }
            }
            if(auxi_area(opTool, opHaveTool, area_sum))
                return true;
            for(int j = 0; j < tmpop.size(); ++j){
                opHaveTool[tmpop[j]] = false;
                OpTool[tmpop[j]] = -1;
            }
            area_sum -= p->area;
            }
        }
    return true;
}

void cal(){
    OpTool = new int[NOpType];
    for(int i = 0; i < NOpType; ++i) OpTool[i] = -1;

    // allocation
    vector<vector<int> > opTool; // 表示每个optype的工具resource
    vector<int> opFrequency; // 表示每个optype的频率：即有多少操作是这个optype
    vector<bool> opHaveTool; // 表示每个optype有没有工具用
    for(int i = 0; i < NOpType; ++i){ // 初始化
        opTool.push_back(vector<int>{});
        opFrequency.push_back(0);
        opHaveTool.push_back(false);
    }
    for(int i = 0; i < ResourceLibrary.size(); ++i){ // 算opTool
        auto p = &ResourceLibrary[i];
        for(int j = 0; j < p->nr; ++j){
            int optype = p->optypes[j];
            opTool[optype].push_back(i);
        }
    }
    for(int i = 0; i < NOperation; ++i){ // 算opFrequency
        opFrequency[Op[i].optype]++;
    }
    // 最低要求
    int area_sum = 0;
    for(int i = 0; i < NOpType; ++i){
        if(opTool[i].size() == 0 || opFrequency[i] == 0) // 说明是个不用schedule的type/没有op是这个type
            opHaveTool[i] = true;
        else if(opTool[i].size() == 1 && opHaveTool[i] == false){ // 说明这个type只有一种工具
            int tool = opTool[i][0]; auto p = &ResourceLibrary[tool];
            area_sum += p->area; // load和store保证是0
            for(int j = 0; j < p->nr; ++j){
                int tmp = p->optypes[j];
                if(opHaveTool[tmp] == false){ // 已经有tool的type不去打扰它
                    opHaveTool[tmp] = true;
                    OpTool[tmp] = tool;
                }
            }
        }
    }
    if(!auxi_area(opTool, opHaveTool, area_sum))
        cout << "无最小覆盖 arealimit有问题" << endl;
    // cout << "area_sum:" << area_sum << " area_limit:" << AreaLimit << endl;
    // for(int i = 0; i < NOpType; ++i)
    //     cout << opFrequency[i] << " ";
    // cout << endl;
    // for(int i = 0; i < NOpType; ++i)
    //     cout << opHaveTool[i] << " ";
    // cout << endl;
    // for(int i = 0; i < NOpType; ++i){
    //     cout << "I:" << i << endl;
    //     for(int j = 0; j < opTool[i].size(); ++j)
    //         cout << opTool[i][j] << " ";
    //     cout << endl;
    // }
    // cout << endl;
    for(int i = 0; i < NOpType; ++i) // 用最低要求初始化每个resource的num
        if(OpCategory[i] >= 6)
            ResourceLibrary[OpTool[i]].num = 1;
    
    // cout << "area_sum:" << area_sum << " area_limit:" << AreaLimit << endl;
    // for(int i = 0; i < NOpType; ++i)
    //     cout << "optype:" << i << " tool:" << OpTool[i] << " freq:" << opFrequency[i] << endl;
    // 按照arealimit启发式增加一些resource
    vector<int> freq(opFrequency);
    sort(freq.rbegin(), freq.rend()); auto p = freq.begin();
    while(true){
        if(*p < 4){
            p = freq.begin();
            continue;
        }
        auto q = find(opFrequency.begin(), opFrequency.end(), *p);
        int type = distance(opFrequency.begin(), q);
        if(OpCategory[type] < 6){
            // category是1 2 3 4 5
            p++;
            continue;
        }
        auto & inf = ResourceLibrary[OpTool[type]];
        if(area_sum + inf.area > AreaLimit)
            break;
        inf.num++;
        area_sum += inf.area;
        p++;
    }
    // cout << "area_sum:" << area_sum << " area_limit:" << AreaLimit << endl;
    // for(int i = 0; i < NResourceType; ++i)
    //     cout << "resource:" << i << " num:" << ResourceLibrary[i].num << endl;

    // 确定基本块的schedule顺序
    def_order();

    // binding包含在schedule里面了
    for(int i = 0; i < NBlock; ++i)
        schedule(block_order[i]);
}