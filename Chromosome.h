#ifndef FJSP_CHROMOSOME_H
#define FJSP_CHROMOSOME_H

#include <vector>
#include "Job.h"

using std::vector;

class Chromosome {
    vector<Job> segment_jobs;
    map<int, int> job_index_map; // Key means the item id, Value mean the index of first job of the item
    vector<Job> chromosome_jobs;
    vector<int> chromosome_machine_id;
public:
    Chromosome() = default;
    Chromosome(const Chromosome& other) = default;
    explicit Chromosome(const vector<Job>& jobs);
    void cross(Chromosome& other, double probability);
    void mutation(double probability);
    [[nodiscard]] int make_span() const;
private:
    void init_chromosome();
    void fix_machine_id(int origin_index);
    void fix_job(int origin_index, int target_item_id);
};

class Population {
    vector<int> iteration_record;
    vector<Chromosome> chromosomes;
    vector<Chromosome> filial_generation; // The children of the chromosomes
    vector<double> fitness;
    double average_fitness{};
    double max_fitness{};
    double min_fitness{};
    double sum_fitness{};

public:
    explicit Population(const vector<Job>& jobs);
    Population(const Population& other) = default;
    void iteration(int iteration_time);

    void output_result();

private:
    [[nodiscard]] double probability_cross(int index_individual1, int index_individual2) const;
    [[nodiscard]] double probability_mutation(int index_individual) const;
    [[nodiscard]] double probability_accept(int index_individual) const;
    void refresh_fitness();
    void select_individuals();
    void cross();
    void mutation();
    void print_population(const char* title) const;
    void clear_repeat();
};


#endif //FJSP_CHROMOSOME_H
