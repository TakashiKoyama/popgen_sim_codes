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
  std::vector<Individual> individuals;
  std::mt19937 mt;
  //std::vector<double> fitness;
    //int generationToFix = 1;
  public:
  //Constructor宣言
  Population(const int popSize_input, const double maf_input);

  //関数宣言
  //std::vector<int> getIndividualID();
  std::vector<int> getIndividualGT();
  void nextGeneration(const double maFitness_input);
  void iterateToFix(const double maFitness_input, int& generationToFix_out, bool& fixedAllele_out);
  double getIndividualS(const Individual individual_input);
};

#endif
