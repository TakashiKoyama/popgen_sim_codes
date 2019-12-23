#pragma once
#ifndef SAMPLES_HPP
#define SAMPLES_HPP

#include "Population.hpp"

class Samples{
private:
  std::vector<Individual> samples;
  std::vector<std::vector<bool>> sampleSegAlleles;
  unsigned int sampleSize;
  unsigned int numOfSegSites;

public:
  //Constructor
  Samples(){}
  Samples(const int sampleSize, const Population pop_in);
  //Functions
  void printIndividualAlleles();
  std::vector<double> getMAFs();
  double getHeterozygosity();
  unsigned int getS();
  double getPi();
};
#endif
