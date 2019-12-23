#include <iostream>
#include <vector>
#include <fstream>
#include "Individual.hpp"
#include "Population.hpp"
#include "Samples.hpp"

int main(){
  int n = 1000;//population size
  double p = 0.0;//minor allele frequency
  int loci = 100;//number of loci
  double mu = 0.00001;//mutation rate(/site/generation)
  double r = 0.0;//recombination rate(/site/generation)
  double fitness = 1.00;//minor allele fitness
  int trial = 1000;//number of trials
  int numOfSamples = 10;
  std::ofstream outfile("out.txt", std::ios_base::app);
  outfile << "Pop size: " << n
  << ", Number of loci: " << loci
  << ", MAF: " << p
  << ", Mutation rate: " << mu
  << ", Recombination rate: " << r
  << ", MA Fitness: " << fitness
  << ", Number of E(H) calculations: " << trial << "\n";
  outfile << "Generation (*2n)\tHeterozygosity\tSegregatingSites\tPi" << std::endl;
  Population population(n, loci, p, mu, r);
  for (int i = 1; i <= (2 * n * trial); ++i){
    population.nextGeneration(fitness);
    //get heterozygosity at each 2n generations
    if (i % (2 * n) == 0){
      Samples samples(numOfSamples, population);
      samples.printIndividualAlleles();
      double h = samples.getHeterozygosity();
      double pi = samples.getPi();
      unsigned int s = samples.getS();
      //double r2;
      outfile << i/(2 * n) << "\t" << h << "\t" << s << "\t" << pi << "\n";
    }
  }
  outfile.close();

  //std::cout << "10th generation" << std::endl;
  //population.printIndividualAlleles();

  /*int generationSum = 0;
  int MAfixedEvent = 0;
  std::ofstream outfileGeneration("fixedGenerationList.txt", std::ios_base::app);

  for (int i = 1; i <= trial; ++i){
    Population population(n, loci, p, mu);//populationクラス定義
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
  std::cout << "Fixation probability: " << (MAfixedEvent * 1.0)/trial << std::endl;
  if (MAfixedEvent) std::cout << "Ave. fixed generation of mutant: " << (generationSum/MAfixedEvent) << std::endl;
*/}
