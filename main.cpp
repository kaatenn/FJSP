#include <string>
#include <vector>

#include "json.hpp"
#include "Job.h"
#include "Chromosome.h"
#include "chrono"
#include "iostream"

using std::ifstream;
using std::string;
using std::vector;

using nlohmann::json;

string read_json(const string& filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        exit(-1);
    }
    string json_str;
    string line;
    while (getline(fin, line)) {
        json_str += line;
    }
    return json_str;
}

void analyse_json(vector<Job>& jobs, const string& jobs_json_str) {
    json j_jobs = json::parse(jobs_json_str);

    for (const auto& job: j_jobs["job_list"]) {
        Job j;
        j.item_id = job["item_id"];
        j.job_id = job["job_id"];
        for (const auto& machine_time: job["machine_map"]) {
            j.machine_time_map[machine_time["machine_id"]] = machine_time["work_time"];
        }
        jobs.push_back(j);
    }
}

int main() {
    string file_name = "object.json";
    string job_json_str = read_json(file_name);
    vector<Job> jobs;
    analyse_json(jobs, job_json_str);
    auto start_time = std::chrono::high_resolution_clock::now();
    Population population(jobs);
    population.iteration(1500);
    population.output_result();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    std::cout << "Duration: " << duration.count() << " seconds." << std::endl;
    return 0;
}
