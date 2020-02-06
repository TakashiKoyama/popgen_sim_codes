#include "Population.hpp"
#include "Individual.hpp"
#include <algorithm>

//Constructor
//1st generation
Population::Population( const int popSize_input,
                        const int loci_input,
                        const double maf_input,
                        const double mu_input,
                        const double r_input){
  popSize = popSize_input;
  loci = loci_input;
  individuals.reserve(popSize);
  // check inputs
  //if (maf_input*size_input < 1)

  for (int i = 1; i <= popSize; ++i){
    Individual individual(0,0,0);
    if (i <= maf_input*popSize) individuals.emplace_back(i, loci, 1);//minor allele is 1
    else individuals.emplace_back(i, loci, 0);//major allele is 0
  }
  //random generator
  std::random_device engine;
  std::mt19937 mt_gen(engine());
  mt = mt_gen;
  std::uniform_int_distribution<unsigned> uniformDist_ind_gen(0, popSize - 1);
  uniformDist_ind = uniformDist_ind_gen;
  lambda_mu = mu_input * loci * popSize;
  lambda_rho = r_input * loci * popSize / 2;
}

//Function
//getter
//individual ID or GT
//std::vector<int> Population::getIndividualID(){
std::vector<int> Population::getIndividualGT(){
  std::vector<int> id_vector;
  id_vector.reserve(popSize);
  for(int i = 0; i < popSize; ++i){
    //id_vector.push_back(individuals.at(i).getID());
    id_vector.push_back(individuals.at(i).getGT());
  }
  return id_vector;
}

std::vector<Individual> Population::getSamples(unsigned int sampleSize_in){
  std::vector<Individual> samples_out;
  samples_out.reserve(sampleSize_in);
  for (unsigned int i = 0; i < sampleSize_in; ++i){
    samples_out.push_back(individuals.at(uniformDist_ind(mt)));
  }
  return samples_out;
}

void Population::printIndividualAlleles(){
  for (auto individualobj : individuals){
    std::vector<bool> alleles = individualobj.getAlleles();
    for (auto allele : alleles){
      std::cout << allele;
    }
    std::cout << ", ";
  }
  std::cout << std::endl;
}

//setter
void Population::setMutations(){
  //set the number of mutations per generation
  std::poisson_distribution<> dist(lambda_mu);
  int mutations = dist(mt);
  //asign the mutations
  std::uniform_int_distribution<unsigned> uniformDist_locus(0, loci - 1);
  for (int i = 0; i < mutations; ++i){
    int mutIndividual = uniformDist_ind(mt);
    int mutLocus = uniformDist_locus(mt);
    individuals.at(mutIndividual).toggleAlleles(mutLocus);
  }
}

void Population::setRecombinations(){
  //set the number of recombination per generation
  std::poisson_distribution<> dist(lambda_rho);
  int recombinations = dist(mt);
  //assign the recombinations
  std::uniform_int_distribution<unsigned> uniformDist_locus(1, loci-1);
  for (int i = 0; i < recombinations; ++i){
    int recIndividual1 = uniformDist_ind(mt);
    int recIndividual2 = uniformDist_ind(mt);
    if (recIndividual1 == recIndividual2){
      --i;
      continue;
    }
    std::vector<bool> alleles1 = individuals.at(recIndividual1).getAlleles();
    std::vector<bool> alleles2 = individuals.at(recIndividual2).getAlleles();
    int recLocus = uniformDist_locus(mt);
    //swap alleles
    std::swap_ranges(alleles1.begin(), alleles1.begin()+recLocus, alleles2.begin());
    //set new alleles to the individuals
    individuals.at(recIndividual1).replaceAlleles(alleles1);
    individuals.at(recIndividual2).replaceAlleles(alleles2);
  }
}

//Construct next generation
void Population::nextGeneration(const double maFitness_input){
  std::vector<Individual> pop_out;
  pop_out.reserve(popSize);
  std::vector<double> fitness;
  fitness.reserve(popSize);
  //Construct fitness vector
  for (int i = 0; i < popSize; ++i){
    std::vector<bool> alleles = individuals.at(i).getAlleles();
    if (alleles.at(0)) fitness.push_back(maFitness_input);
    else fitness.push_back(1.0);
  }
  //pickup individuals with selection
  std::discrete_distribution<std::size_t> discreteDist_ind(fitness.begin(), fitness.end());
  for (int i = 0; i < popSize; ++i){
    std::size_t ind = discreteDist_ind(mt);
    pop_out.push_back(individuals.at(ind));
  }
  individuals = pop_out;
  //set mutations and recombinations on next generation
  setMutations();
  setRecombinations();
}

//Iterator until fixation
void Population::iterateToFix(const double maFitness_input,
                              int& generationToFix_out,
                              bool& fixedAllele_out){
  nextGeneration(maFitness_input);
  //int generationToFix;
  for(generationToFix_out = 2; ; ++generationToFix_out){
    //printIndividualAlleles();
    std::vector<bool> id_vector;
    id_vector.reserve(individuals.size());
    for (auto individualobj : individuals){
      //std::vector<int> id_vector = getIndividualID();
      //std::vector<int> id_vector = getIndividualGT();
      id_vector.push_back(individualobj.getAlleles().at(0));
    }
    //sort(id_vector.begin(),id_vector.end());
    int diversity = 0;
    for (auto id : id_vector){
      diversity += id;
    }
    //for(int i = 0; i != (size-1); ++i){
      //if (id_vector.at(i) != id_vector.at(i+1)) ++diversity;
    //}
    //std::cout << id_vector.at(size - 1) << std::endl;
    if (diversity == 0 || diversity == popSize){
      fixedAllele_out = diversity/popSize;
      break;
    }
    else nextGeneration(maFitness_input);
  }
}
