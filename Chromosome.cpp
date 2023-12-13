#include "Chromosome.h"

#include <random>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <set>
#include <nlohmann/json.hpp>

// #define DEBUG
#define POPULATION_CAPACITY 3000
#define SURVIVAL_CAPACITY 100
#define MUTATION_OFFSET 0.1
#define CROSS_OFFSET 0.1

using std::mt19937;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::random_device;
using std::cout;
using std::set;
using nlohmann::json;

using std::max_element;
using std::min_element;
using std::accumulate;
using std::shuffle;
using std::for_each;
using std::sort;
using std::find_if;

//-------------------- Chromosome --------------------

Chromosome::Chromosome(const vector<Job> &jobs) {
    this->segment_jobs = jobs;
    // Sort jobs by hash code f(x) = item_id * 1000 + job_id
    sort(segment_jobs.begin(), segment_jobs.end(), [](const Job &left, const Job &right) -> bool {
        return left.item_id * 100000 + left.job_id < right.item_id * 1000 + right.job_id;
    });
    // Save the index of the first job
    for (int i = 0; i < segment_jobs.size(); ++i) {
        if (job_index_map.count(segment_jobs[i].item_id) == 0) {
            job_index_map[segment_jobs[i].item_id] = i;
        }
    }
    init_chromosome();
}

int Chromosome::make_span() const {
    map<int, int> machine_working_time_map;
    map<int, int> machine_item_map;
    int span = 0;
    for (int i = 0; i < this->chromosome_jobs.size(); ++i) {

        int current_machine = this->chromosome_machine_id[i];

        if (machine_working_time_map.count(current_machine) != 0) {
            // Calculate the time until the current job begin.
            int last_time = machine_working_time_map[current_machine]; // The last time of the job on the machine
            span += last_time;
            vector<int> available_machine;
            for_each(machine_working_time_map.begin(), machine_working_time_map.end(), [&](pair<const int, int> p) {
                p.second -= last_time;
                if (p.second <= 0) {
                    // The time of the job on the machine of course will be 0
                    available_machine.push_back(p.first);

                }
            });
            for_each(available_machine.begin(), available_machine.end(), [&machine_working_time_map, &machine_item_map]
            (const int
                                                                                                     machine_id) {
                machine_item_map.erase(machine_id);
                machine_working_time_map.erase(machine_id);
            });
        }

        // If item isn't available, calculate the time until the current item available.
        auto it = find_if(machine_item_map.begin(), machine_item_map.end(), [&](pair<const int, int> p) {
            return p.second == this->chromosome_jobs[i].item_id;
        });
        if (it != machine_item_map.end()) {
            int last_time = machine_working_time_map[it->first];
            span += last_time;
            vector<int> available_machine;
            for_each(machine_working_time_map.begin(), machine_working_time_map.end(), [&](pair<const int, int> p) {
                p.second -= last_time;
                if (p.second <= 0) {
                    available_machine.push_back(p.first);
                }
            });
            for_each(available_machine.begin(), available_machine.end(), [&machine_working_time_map, &machine_item_map]
            (const int
                                                                                                     machine_id) {
                machine_item_map.erase(machine_id);
                machine_working_time_map.erase(machine_id);
            });
        }

        machine_working_time_map[current_machine] = this->chromosome_jobs[i].machine_time_map.at(current_machine);
        machine_item_map[current_machine] = this->chromosome_jobs[i].item_id;
    }

    if (!machine_working_time_map.empty()) {
        // Add the max time to the total span
        span += max_element(machine_working_time_map.begin(), machine_working_time_map.end(),
                            [](const pair<const int, int> &left, const pair<const int, int> &right) -> bool {
                                return left.second < right.second;
                            })->second;
    }
    return span;
}

void Chromosome::init_chromosome() {
    //----- Generate the chromosome -----
    int chromosome_size = static_cast<int>(this->segment_jobs.size());
    // The random engine for randomly selecting
    mt19937 mt_rand(random_device{}());

    // Shuffle the element in the chromosome to generate a new chromosome under rules
    vector<int> item_repeated; // Generate a list of items repeating by the number of jobs
    for_each(this->segment_jobs.begin(), this->segment_jobs.end(), [&item_repeated](const Job &job) {
        item_repeated.push_back(job.item_id);
    });
    shuffle(item_repeated.begin(), item_repeated.end(), mt_rand);
    map<int, int> current_job_index;
    for (const int item_id: item_repeated) {
        // Create a pair if the item haven't been arranged
        if (current_job_index.count(item_id) == 0) {
            current_job_index[item_id] = 0;
        }
        this->chromosome_jobs.push_back(this->segment_jobs[job_index_map[item_id] + (current_job_index[item_id]++)]);
    }

    for (int i = 0; i < chromosome_size; ++i) {
        // Randomly select a machine
        uniform_int_distribution<int> dist(0, static_cast<int>(chromosome_jobs[i].machine_time_map.size() - 1));
        vector<int> machine_ids;
        for_each(chromosome_jobs[i].machine_time_map.begin(), chromosome_jobs[i].machine_time_map.end(),
                 [&machine_ids](const
                                pair<const int, int> &p) {
                     machine_ids.push_back(p.first);
                 });
        int selected_machine_id = machine_ids[dist(mt_rand)];
        this->chromosome_machine_id.push_back(selected_machine_id);
    }
}

void Chromosome::cross(Chromosome &other, double probability) {
    mt19937 mt_rand(random_device{}());
    uniform_real_distribution<double> dist(0, 1);

    // Machine cross
    for (int i = 0; i < this->chromosome_jobs.size(); ++i) {
        if (dist(mt_rand) < probability) {
            // Swap the machine id
            int temp = this->chromosome_machine_id[i];
            this->chromosome_machine_id[i] = other.chromosome_machine_id[i];
            other.chromosome_machine_id[i] = temp;

            // Fix the machine id if it is illegal
            this->fix_machine_id(i);
            other.fix_machine_id(i);
        }
    }

    // Job cross
    for (int i = 0; i < this->chromosome_jobs.size(); ++i) {
        if (dist(mt_rand) < probability) {
            // Swap the job
            if (this->chromosome_jobs[i].item_id != other.chromosome_jobs[i].item_id) {
                this->fix_job(i, other.chromosome_jobs[i].item_id);
                this->fix_job(i, this->chromosome_jobs[i].item_id);
            }
        }
    }
}

void Chromosome::fix_machine_id(int origin_index) {
    if (this->chromosome_jobs[origin_index].machine_time_map.count(this->chromosome_machine_id[origin_index]) == 0) {
        vector<int> machine_ids;
        // If the machine id is illegal, find a legal machine id
        mt19937 mt_rand(random_device{}());
        for_each(this->chromosome_jobs[origin_index].machine_time_map.begin(),
                 this->chromosome_jobs[origin_index].machine_time_map.end(),
                 [&machine_ids](const
                                pair<const int, int> &p) {
                     machine_ids.push_back(p.first);
                 });
        uniform_int_distribution<int> dist(0, static_cast<int>(machine_ids.size() - 1));
        this->chromosome_machine_id[origin_index] = machine_ids[dist(mt_rand)];
    }
}

void Chromosome::fix_job(int origin_index, int target_item_id) {
    int j = origin_index;
    while (j < this->chromosome_jobs.size() && this->chromosome_jobs[j].item_id != target_item_id) {
        ++j;
    }
    if (j < this->chromosome_jobs.size()) {
        // Fix the job id between the exchanged jobs
        int last_job_id = this->chromosome_jobs[origin_index].job_id;
        for (int k = origin_index; k < j; k++) {
            if (this->chromosome_jobs[k].item_id == this->chromosome_jobs[origin_index].item_id) {
                last_job_id = this->chromosome_jobs[k].job_id;
                this->chromosome_jobs[k].job_id--;
                this->fix_machine_id(k);
            }
        }
        if (last_job_id != this->chromosome_jobs[origin_index].job_id) {
            this->chromosome_jobs[origin_index].job_id = last_job_id + 1;
            this->fix_machine_id(origin_index);
        }
        Job temp = this->chromosome_jobs[origin_index];
        this->chromosome_jobs[origin_index] = this->chromosome_jobs[j];
        this->chromosome_jobs[j] = temp;
        int temp_machine_id = this->chromosome_machine_id[origin_index];
        this->chromosome_machine_id[origin_index] = this->chromosome_machine_id[j];
        this->chromosome_machine_id[j] = temp_machine_id;
    } else {
        j = origin_index;
        while (j >= 0 && this->chromosome_jobs[j].item_id != target_item_id) {
            --j;
        }
        if (j >= 0) {
            int last_job_id = this->chromosome_jobs[origin_index].job_id;
            for (int k = origin_index; k > j; k--) {
                if (this->chromosome_jobs[k].item_id == this->chromosome_jobs[origin_index].item_id) {
                    last_job_id = this->chromosome_jobs[k].job_id;
                    this->chromosome_jobs[k].job_id++;
                    this->fix_machine_id(k);
                }
            }
            if (last_job_id != this->chromosome_jobs[origin_index].job_id) {
                this->chromosome_jobs[origin_index].job_id = last_job_id - 1;
                this->fix_machine_id(origin_index);
            }
            Job temp = this->chromosome_jobs[origin_index];
            this->chromosome_jobs[origin_index] = this->chromosome_jobs[j];
            this->chromosome_jobs[j] = temp;
            int temp_machine_id = this->chromosome_machine_id[origin_index];
            this->chromosome_machine_id[origin_index] = this->chromosome_machine_id[j];
            this->chromosome_machine_id[j] = temp_machine_id;
        }
    }
}

void Chromosome::mutation(double probability) {
    mt19937 mt_rand(random_device{}());
    uniform_real_distribution<double> dist(0, 1);
    for (int i = 0; i < this->chromosome_jobs.size(); ++i) {
        if (dist(mt_rand) < probability) {
            // Reselect the job id
            set<int> job_ids;
            for_each(this->chromosome_jobs.begin(), this->chromosome_jobs.end(),
                     [&job_ids](const Job &job) {
                         job_ids.insert(job.job_id);
                     });
            uniform_int_distribution<int> dist_job(0, static_cast<int>(job_ids.size() - 1));
            int selected_index = dist_job(mt_rand);
            auto it = job_ids.begin();
            for (int j = 0; j < selected_index; ++j) {
                it++;
            }
            int selected_job_id = *it;
            if (selected_job_id != this->chromosome_jobs[i].job_id)
                this->fix_job(i, selected_job_id);
            // Reselect machine
            vector<int> machine_ids;
            for_each(this->chromosome_jobs[i].machine_time_map.begin(),
                     this->chromosome_jobs[i].machine_time_map.end(),
                     [&machine_ids](const
                                    pair<const int, int> &p) {
                         machine_ids.push_back(p.first);
                     });
            uniform_int_distribution<int> dist_machine(0, static_cast<int>(machine_ids.size() - 1));
            this->chromosome_machine_id[i] = machine_ids[dist_machine(mt_rand)];
        }

        if (dist(mt_rand) < probability) {
            // Reselect the machine id
            vector<int> machine_ids;
            for_each(this->chromosome_jobs[i].machine_time_map.begin(),
                     this->chromosome_jobs[i].machine_time_map.end(),
                     [&machine_ids](const
                                    pair<const int, int> &p) {
                         machine_ids.push_back(p.first);
                     });
            uniform_int_distribution<int> dist_machine(0, static_cast<int>(machine_ids.size() - 1));
            this->chromosome_machine_id[i] = machine_ids[dist_machine(mt_rand)];
        }
    }
}

//-------------------- Population --------------------

Population::Population(const vector<Job> &jobs) {

    for (int i = 0; i < SURVIVAL_CAPACITY; ++i) {
        Chromosome root_chromosome(jobs);
        this->chromosomes.push_back(root_chromosome);
    }

    refresh_fitness();
}

void Population::refresh_fitness() {
    this->fitness.clear();
    for (const auto &chromosome: this->chromosomes) {
        this->fitness.push_back(1.0 / chromosome.make_span());
    }
    this->max_fitness = *max_element(this->fitness.begin(), this->fitness.end());
    this->min_fitness = *min_element(this->fitness.begin(), this->fitness.end());
    this->sum_fitness = accumulate(this->fitness.begin(), this->fitness.end(), 0.0);
    this->average_fitness = sum_fitness / static_cast<int>(this->chromosomes.size());
    this->filial_generation.clear();
}

void Population::iteration(int iteration_time) {
    for (int i = 0; i < iteration_time; ++i) {
#ifdef DEBUG
        cout << "Iteration " << i << ": " << endl;
#endif
        this->select_individuals();
#ifdef DEBUG
        this->print_population("Selected individuals");
#endif
        this->cross();
        clear_repeat();
#ifdef DEBUG
        this->print_population("Crossed individuals");
#endif
        this->mutation();
        clear_repeat();
#ifdef DEBUG
        this->print_population("Mutated individuals");
#endif
        int min_span = INT32_MAX;
        for_each(this->chromosomes.begin(), this->chromosomes.end(), [&min_span](Chromosome &chromosome) {
            if (chromosome.make_span() < min_span) {
                min_span = chromosome.make_span();
            }
        });
        this->iteration_record.push_back(min_span);

        // Debug the min span
#ifdef DEBUG
        cout << endl;
        cout << min_span << endl;
#endif
    }
#ifdef DEBUG
    for_each(this->iteration_record.begin(), this->iteration_record.end(), [](int span) {
        cout << span << " ";
    });
#endif
}

double Population::probability_cross(int index_individual1, int index_individual2) const {
    mt19937 mt_rand(random_device{}());
    uniform_real_distribution<double> dist(0, 1);
    double k1 = dist(mt_rand);
    double k2 = dist(mt_rand);
    double f_cross_max = max(this->fitness[index_individual1], this->fitness[index_individual2]);
    if (f_cross_max > this->average_fitness) {
        return k1 * (max_fitness - f_cross_max) / (max_fitness - min_fitness) + CROSS_OFFSET;
    } else {
        return k2;
    }
}

double Population::probability_mutation(int index_individual) const {
    mt19937 mt_rand(random_device{}());
    uniform_real_distribution<double> dist(0, 1);
    double k3 = dist(mt_rand);
    double k4 = dist(mt_rand);
    if (this->fitness[index_individual] > this->average_fitness) {
        return k3 * (max_fitness - this->fitness[index_individual]) / (max_fitness - min_fitness) + MUTATION_OFFSET;
    } else {
        return k4;
    }
}

double Population::probability_accept(int index_individual) const {
    return pow(this->fitness[index_individual], 3) / this->sum_fitness * pow
            (SURVIVAL_CAPACITY, 5);
}

void Population::select_individuals() {
    sort(this->chromosomes.begin(), this->chromosomes.end(), [](const Chromosome &left, const Chromosome &right) -> bool {
        return left.make_span() < right.make_span();
    });
    mt19937 mt_rand(random_device{}());
    uniform_real_distribution<double> dist(0, 1);
    for (int i = 1; i < chromosomes.size(); ++i) {
        // The first chromosome is the best chromosome
        if (dist(mt_rand) < probability_accept(i)) {
            this->filial_generation.push_back(this->chromosomes[i]);
        }
    }
    this->filial_generation.push_back(this->chromosomes[0]);
    this->chromosomes = this->filial_generation;
    refresh_fitness();
}

void Population::cross() {
    // Shuffle the chromosomes in the population
    sort(this->chromosomes.begin(), this->chromosomes.end(), [](const Chromosome &left, const Chromosome &right) -> bool {
        return left.make_span() < right.make_span();
    });
    // Cross the chromosomes
    for (int i = 0; i < this->chromosomes.size(); ++i) {
        int next_index = (i + 1) % static_cast<int>(this->chromosomes.size());
        this->filial_generation.push_back(this->chromosomes[i]);
        this->chromosomes[i].cross(this->chromosomes[next_index], probability_cross(i, next_index));
        this->filial_generation.push_back(this->chromosomes[i]);
        this->filial_generation.push_back(this->chromosomes[next_index]);
        if (filial_generation.size() >= POPULATION_CAPACITY) {
            break;
        }
    }

    this->chromosomes = this->filial_generation;
    refresh_fitness();
}

void Population::mutation() {
    sort(this->chromosomes.begin(), this->chromosomes.end(), [](const Chromosome &left, const Chromosome &right) -> bool {
        return left.make_span() < right.make_span();
    });
    for (int i = 3; i < this->chromosomes.size(); ++i) {
        // The first three chromosomes are the best chromosomes
        this->chromosomes[i].mutation(probability_mutation(i));
    }
}

void Population::print_population(const char* title) const {
    cout << title << ":" << endl;
    for (const auto &chromosome: this->chromosomes) {
        cout << chromosome.make_span() << " ";
    }
    cout << endl;
}

void Population::clear_repeat() {
    // If the same span exists, remove the chromosome with the smaller index
    sort(this->chromosomes.begin(), this->chromosomes.end(), [](const Chromosome &left, const Chromosome &right) -> bool {
        return left.make_span() < right.make_span();
    });
    for (int i = 0; i < this->chromosomes.size(); ++i) {
        for (int j = i + 1; j < this->chromosomes.size(); ++j) {
            if (this->chromosomes[i].make_span() == this->chromosomes[j].make_span()) {
                this->chromosomes.erase(this->chromosomes.begin() + j);
                --j;
            }
        }
    }

}

void Population::output_result() {
    json json_output = this->iteration_record;
    std::ofstream fout("result.json");
    fout << json_output;
}
