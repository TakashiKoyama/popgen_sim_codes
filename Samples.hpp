#pragma once
#ifndef SAMPLES_HPP
#define SAMPLES_HPP

#include "Population.hpp"
#include <cmath>
#include <map>
#include <unordered_map>

class Samples{
private:
  std::vector<Individual> samples;
  std::vector<std::vector<bool>> sampleSegAlleles;
  unsigned int sampleSize;
  unsigned int numOfSegSites;
  std::vector<unsigned int> segSite_positions;

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
  std::vector<std::vector<double>> getHFs();
  //std::map<int, double> getLD();
  void getLD(std::map<unsigned int, double>& r2Table, std::map<unsigned int, unsigned long>& r2Redundancy);
};
#endif
