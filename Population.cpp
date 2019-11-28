#include "Population.hpp"
#include "Individual.hpp"
#include <algorithm>

//Constructor
//1st generation
Population::Population(const int popSize_input, const double maf_input){
  popSize = popSize_input;
  individuals.reserve(popSize_input);
  // check inputs
  //if (maf_input*size_input < 1)

  for (int i = 1; i <= popSize_input; ++i){
    Individual individual(0,0);
    if (i <= maf_input*popSize_input) individuals.emplace_back(i, 1);
    else individuals.emplace_back(i, 0);
  }
  //random generator
  std::random_device engine;
  std::mt19937 mt_gen(engine());
  mt = mt_gen;
}

//Function
//getter
//individual ID or GT
//std::vector<int> Population::getIndividualID(){
std::vector<int> Population::getIndividualGT(){
  std::vector<int> id_vector;
  id_vector.reserve(popSize);
  for(int i = 0; i < popSize; ++i){
    //id_vector.push_back(individuals[i].getID());
    id_vector.push_back(individuals[i].getGT());
  }
  return id_vector;
}

//Construct next generation
void Population::nextGeneration(const double maFitness_input){
  std::vector<Individual> pop_out;
  pop_out.reserve(popSize);
  //fitness.clear();
  std::vector<double> fitness;
  fitness.reserve(popSize);
  //Construct fitness vector
  for (int i = 0; i < popSize; ++i){
    if (individuals[i].getGT()) fitness.push_back(maFitness_input);
    else fitness.push_back(1.0);
  }
  //pickup individuals with no selection
  //std::uniform_int_distribution<unsigned> dist(0,popSize - 1);
  std::discrete_distribution<std::size_t> dist(fitness.begin(), fitness.end());
  for (int i = 0; i < popSize; ++i){
    pop_out.push_back(individuals[dist(mt)]);
  }
  individuals = pop_out;
}

//Iterator until fixation
void Population::iterateToFix(const double maFitness_input, int& generationToFix_out, bool& fixedAllele_out){
  nextGeneration(maFitness_input);
  //int generationToFix;
  for(generationToFix_out = 2; ; ++generationToFix_out){
    //std::vector<int> id_vector = getIndividualID();
    std::vector<int> id_vector = getIndividualGT();
    //sort(id_vector.begin(),id_vector.end());
    int diversity = 0;
    for (std::size_t i = 0; i < id_vector.size(); ++i){
      diversity += id_vector[i];
    }
    //for(int i = 0; i != (size-1); ++i){
      //if (id_vector[i] != id_vector[(i+1)]) ++diversity;
    //}
    //std::cout << id_vector[(size - 1)] << std::endl;
    if (diversity == 0 || diversity == popSize){
      fixedAllele_out = diversity/popSize;
      break;
    }
    else nextGeneration(maFitness_input);
  }
}
