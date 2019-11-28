#pragma once
#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

class Individual{
private:
  int id;
  bool gt;
  double s;
  //std::vector alleles;
public:
  //Constructor
  Individual(const int id_input, const bool gt_input) : id(id_input), gt(gt_input) {}//Constructor
  //Function
  //getter
  const int getID(){ return id; }//get id from individual
  const bool getGT(){ return gt; }
  };

#endif
