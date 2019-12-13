#pragma once
#ifndef SAMPLES_HPP
#define SAMPLES_HPP

#include "Population.hpp"

class Samples{
private:
  std::vector<Individual> samples;
  int sampleSize;
  int loci;

public:
  //Constructor
  Samples(const int sampleSize, const Population pop_in);
  //Functions
  void printIndividualAlleles();
  std::vector<double> getMAFs();
  double getHeterozygosity();
  void getPiAndS(double& pi_out, int& s_out);
};
#endif
