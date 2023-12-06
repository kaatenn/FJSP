//
// Created by 86137 on 2023/12/3.
//

#ifndef FJSP_JOB_H
#define FJSP_JOB_H

#include "fstream"
#include <map>

using namespace std;

struct Job {
    int item_id = 0;
    int job_id = 0;
    map<int, int> machine_time_map;
    Job() = default;
    Job(const Job& other) = default;
};

#endif //FJSP_JOB_H
