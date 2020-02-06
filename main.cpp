#include <iostream>
#include <vector>
#include <fstream>
#include "Individual.hpp"
#include "Population.hpp"
#include "Samples.hpp"

int main(){
  int n = 1000;//population size (nreps in ms)
  double p = 0.0;//minor allele frequency
  int loci = 1001;//number of loci
  double mu = 0.000001;//mutation rate(/site/generation)
  double r = 0.000001;//recombination rate(/site/generation)
  double fitness = 1.00;//minor allele fitness
  int sampling = 100000;
  int sampleSize = 10;//sample size (nsam in ms)
  std::ofstream outfile_parameters("Parameters.txt", std::ios_base::app);
  std::ofstream outfile_summary("SummaryStatistics.txt", std::ios_base::app);
  std::ofstream outfile_r2("PairwiseR2.txt", std::ios_base::app);
  outfile_parameters << "Parameter\tValue"
  << "\nPop_size\t" << n
  << "\nNumber_of_loci\t" << loci
  << "\nMAF\t" << p
  << "\nMutation_rate\t" << mu
  << "\nRecombination_rate\t" << r
  << "\nMA_Fitness\t" << fitness
  << "\nSample_size\t" << sampleSize
  << "\nNumber_of_samplings\t" << sampling << std::endl;
  outfile_summary << "Generation (*2n)\tHeterozygosity\tSegregatingSites\tPi" << std::endl;
  //set objects
  Population population(n, loci, p, mu, r);
  std::map<unsigned int, double> r2Table;
  std::map<unsigned int, unsigned long> r2Redundancy;
  for (int i = 1; i <= (2 * n * sampling); ++i){
    population.nextGeneration(fitness);
    //get summary statistics at each 2n generations
    if (i % (2 * n) == 0){
      Samples samples(sampleSize, population);
      double h = samples.getHeterozygosity();
      double pi = samples.getPi();
      unsigned int s = samples.getS();
      outfile_summary << i/(2 * n) << "\t" << h << "\t" << s << "\t" << pi << std::endl;
      std::cout << "//\nsegsites: " << s << "\n";
      samples.printIndividualAlleles();
      //std::map<int, double> r2Table_tmp = samples.getLD();
      samples.getLD(r2Table, r2Redundancy);
      //sum up each sample R2 table
      /*for (auto i_ptr = r2Table_tmp.begin(); i_ptr != r2Table_tmp.end(); ++i_ptr){
        if(i_ptr->first != -1){
          r2Table[i_ptr->first] += r2Table_tmp.at(i_ptr->first);
          ++r2Redundancy[i_ptr->first];
        }
      }*/
    }
  }
  //out average R2
  for (auto i_ptr = r2Table.begin(); i_ptr != r2Table.end(); ++i_ptr){
    double r2 = r2Table.at(i_ptr->first)/(r2Redundancy.at(i_ptr->first) * 1.0);
    outfile_r2 << i_ptr->first << "\t" << r2;
    outfile_r2 << std::endl;
  }

  outfile_parameters.close();
  outfile_summary.close();
  outfile_r2.close();

  //std::cout << "10th generation" << std::endl;
  //population.printIndividualAlleles();

  /*int generationSum = 0;
  int MAfixedEvent = 0;
  std::ofstream outfileGeneration("fixedGenerationList.txt", std::ios_base::app);

  for (int i = 1; i <= sampling; ++i){
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
  std::cout << "Fixation probability: " << (MAfixedEvent * 1.0)/sampling << std::endl;
  if (MAfixedEvent) std::cout << "Ave. fixed generation of mutant: " << (generationSum/MAfixedEvent) << std::endl;
*/}
