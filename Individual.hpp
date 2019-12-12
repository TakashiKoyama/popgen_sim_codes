#pragma once
#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <vector>

class Individual{
private:
  int id;
  bool gt;
  std::vector<bool> alleles;
  //double s;
  //std::vector alleles;
public:
  // Constructor
  Individual( const int id_input,
              const int numOfLoci_input,
              const bool allele_input);
  //Function
  //getter
  int getID(){ return id; }//get id from individual object
  bool getGT(){ return gt; }
  const std::vector<bool> getAlleles(){ return alleles; }// get allele vector from individual object
  //setter
  void toggleAlleles(const int locus_input){ alleles.at(locus_input) = !alleles.at(locus_input); }
  void replaceAlleles(const std::vector<bool> alleles_input){ alleles = alleles_input; }
};

#endif
