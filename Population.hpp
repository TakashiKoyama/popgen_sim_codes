#pragma once
#ifndef POPULATION_HPP
#define POPULATION_HPP

#include "Individual.hpp"
#include <vector>
#include <random>
#include <iostream>
#include <string>

class Population{
private:
  int popSize;
  int loci;
  std::vector<Individual> individuals;
  std::mt19937 mt;
  double lambda_mu;
  double lambda_rho;
  std::uniform_int_distribution<unsigned> uniformDist_ind;
  //std::vector<double> fitness;
    //int generationToFix = 1;
public:
  //Constructor
  Population( const int popSize_input,
              const int loci_input,
              const double maf_input,
              const double mu_input,
              const double r_input);

  //Functions
  //std::vector<int> getIndividualID();
  std::vector<int> getIndividualGT();
  void getSamples(std::vector<Individual>& samples_out);
  void printIndividualAlleles();
  void setMutations();
  void setRecombinations();
  void nextGeneration(const double maFitness_input);
  void iterateToFix(const double maFitness_input,
                    int& generationToFix_out,
                    bool& fixedAllele_out);
  std::vector<double> getMAFs();
  double getHeterozygosity();
  void getStatistics(const int samples_in, double& pi_out, int& s_out);
};

#endif
