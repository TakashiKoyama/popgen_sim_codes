#include "Population.hpp"
#include "Individual.hpp"
#include "Samples.hpp"
//Constructor
//Pickup samples
Samples::Samples(const int sampleSize_in, Population pop_in){
  sampleSize = sampleSize_in;
  samples = pop_in.getSamples(sampleSize);
  loci = samples.at(0).getAlleles().size();
}

//Functions
void Samples::printIndividualAlleles(){
  for (auto sampleobj : samples){
    std::vector<bool> alleles = sampleobj.getAlleles();
    for (auto allele : alleles){
      std::cout << allele;
    }
    std::cout << ", ";
  }
  std::cout << std::endl;
}

std::vector<double> Samples::getMAFs(){
  std::vector<unsigned int> alleleCounts(loci);
  std::vector<double> mafs_out(loci);
  //get minor allele count at each locus
  for (int i = 0; i < sampleSize; ++i){
    for (int j = 0; j < loci; ++j){
      alleleCounts.at(j) += samples.at(i).getAlleles().at(j);
    }
  }
  //calculate allele frequency at each locus
  for (unsigned int i = 0; i < alleleCounts.size(); ++i){
    mafs_out.at(i) = (alleleCounts.at(i)*1.0)/sampleSize;
  }
  return mafs_out;
}

double Samples::getHeterozygosity(){
  double heterozygosity_sum = 0.0;
  std::vector<double> mafs = getMAFs();
  for (auto maf : mafs){
    double heterozygosity = 2 * maf * (1 - maf);
    heterozygosity_sum += heterozygosity;
  }
  return heterozygosity_sum/loci;
}

void Samples::getPiAndS(double& pi_out, int& s_out){
  int k = 0;//total number of differences
  pi_out = 0;
  s_out = 0;
  std::vector<bool> s_vec(loci);
  int nC2 = sampleSize * (sampleSize - 1) / 2;

  for (unsigned int i = 0; i < (samples.size() - 1); ++i){
    for (unsigned int j = i + 1; j < samples.size(); ++j){
      for (int m = 0; m < loci; ++m){
        if (samples.at(i).getAlleles().at(m) != samples.at(j).getAlleles().at(m)){
          ++k;
          s_vec.at(m) = 1;
        }
      }
    }
  }
  pi_out = (k * 1.0) / nC2;
  for (auto s : s_vec) s_out += s;
}
