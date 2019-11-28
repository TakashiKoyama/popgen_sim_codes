#include <iostream>
#include <vector>
#include <fstream>
#include "Individual.hpp"
#include "Population.hpp"
//using namespace std;

int main(){
  int n = 100;//population size
  double p = 0.01;//minor allele frequency
  double fitness = 1.06;//minor allele fitness
  int trial = 1000000;//number of trials
  std::cout << "pop size: " << n << ", MAF: " << p << ", MA Fitness: " << fitness << ", Num of Trials: " << trial << std::endl;

  int generationSum = 0;
  int MAfixedEvent = 0;
  std::ofstream outfileGeneration("fixedGenerationList.txt", std::ios_base::app);

  for (int i = 1; i <= trial; ++i){
    Population population(n, p);//populationクラス定義
    //std::vector<int> id_vectorInitial = population.getIndividualID();
    //std::vector<int> id_vectorInitial = population.getIndividualGT();
    int generationToFix;
    bool fixedAllele;
    population.iterateToFix(fitness, generationToFix, fixedAllele);
    //std::vector<int> id_vectorFixed = population.getIndividualID();
    //std::vector<int> id_vectorFixed = population.getIndividualGT();
    //generationSum += generationToFix;

    if (fixedAllele == 0){
      //--i;
    } else if (fixedAllele == 1) {
      ++MAfixedEvent;
      generationSum += generationToFix;
      outfileGeneration << generationToFix << "\n";
    } else std::cerr << "data exception\n";
  }
  outfileGeneration.close();
  std::cout << "Num. of fixed events of mutant: " << MAfixedEvent << std::endl;
  std::cout << "Fixation probability: " << MAfixedEvent/trial << std::endl;
  std::cout << "Ave. fixed generation of mutant: " << (generationSum/MAfixedEvent) << std::endl;
}
