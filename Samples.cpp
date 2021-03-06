#include "Population.hpp"
#include "Individual.hpp"
#include "Samples.hpp"
//Constructor
//Pickup samples
Samples::Samples(const int sampleSize_in, Population pop_in){
  sampleSize = sampleSize_in;
  samples = pop_in.getSamples(sampleSize);
  int loci = samples.at(0).getAlleles().size();
  //anotate polymorphic loci
  std::vector<bool> segSites_vec(loci);
  for (unsigned int i = 0; i < (samples.size() - 1); ++i){
    for (unsigned int j = i + 1; j < samples.size(); ++j){
      for (int m = 0; m < loci; ++m){
        if (samples.at(i).getAlleles().at(m) != samples.at(j).getAlleles().at(m)){
          segSites_vec.at(m) = 1;
        }
      }
    }
  }
  numOfSegSites = 0;
  for (auto segSite : segSites_vec) numOfSegSites += segSite;

  //record positions of the polymorphic sites
  segSite_positions.reserve(numOfSegSites);
  for (unsigned int i = 0; i < segSites_vec.size(); ++i){
    if(segSites_vec.at(i)) segSite_positions.push_back(i);
  }

  //pickup polymorphic loci
  if (numOfSegSites > 0){//check polymorphism
    sampleSegAlleles.reserve(sampleSize);
    for (unsigned int i = 0; i < samples.size(); ++i){
      std::vector<bool> segAlleles;
      segAlleles.reserve(numOfSegSites);
      for (unsigned int j = 0; j < segSites_vec.size(); ++j){
        if (segSites_vec.at(j)) segAlleles.push_back(samples.at(i).getAlleles().at(j));
      }
      sampleSegAlleles.push_back(segAlleles);
    }
  } //else{
    //throw std::runtime_error("no polymorphic site in samples");
  //}
}

//Functions
void Samples::printIndividualAlleles(){
  if (numOfSegSites > 0){
    //print positions
    std::cout << "positions:";
    for (auto pos : segSite_positions){
      std::cout << " " << pos;
    }
    std::cout << std::endl;
    //print alleles
    for (auto segAlleles : sampleSegAlleles){
      for (auto segAllele : segAlleles) std::cout << segAllele;
      std::cout << std::endl;
    }
  }
}

std::vector<double> Samples::getMAFs(){
  std::vector<unsigned int> alleleCounts(numOfSegSites);
  std::vector<double> mafs_out(numOfSegSites);
  //get minor allele count at each locus
  for (unsigned int i = 0; i < sampleSize; ++i){
    for (unsigned int j = 0; j < numOfSegSites; ++j){
      alleleCounts.at(j) += sampleSegAlleles.at(i).at(j);
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
  if (numOfSegSites > 0){
    std::vector<double> mafs = getMAFs();
    for (auto maf : mafs){
      double heterozygosity = 2 * maf * (1 - maf);
      heterozygosity_sum += heterozygosity;
    }
    return heterozygosity_sum/numOfSegSites;
  } else{//when no polymorphic site
    return heterozygosity_sum;
  }
}

unsigned int Samples::getS(){ return numOfSegSites; }

double Samples::getPi(){
  int k = 0;//total number of differences
  double pi_out = 0;
  int nC2 = sampleSize * (sampleSize - 1) / 2;

  if (numOfSegSites > 0){
    for (unsigned int i = 0; i < (sampleSegAlleles.size() - 1); ++i){
      for (unsigned int j = i + 1; j < sampleSegAlleles.size(); ++j){
        for (unsigned int m = 0; m < numOfSegSites; ++m){
          if (sampleSegAlleles.at(i).at(m) != sampleSegAlleles.at(j).at(m)) ++k;
        }
      }
    }
  }
  return pi_out = (k * 1.0) / nC2;
}

std::vector<std::vector<double>> Samples::getHFs(){
  std::vector<std::vector<double>> DDhf_vec; //vector of derived-derived haplotype frequencies
  std::vector<double> mafs = getMAFs(); //use to check haplotype frequency calculation
  //calculate haplotype frequency
  DDhf_vec.reserve(numOfSegSites - 1);
  for (unsigned int i = 0; i < (numOfSegSites - 1); ++i){
    std::vector<double> DDhfs_tmp;
    DDhfs_tmp.reserve(numOfSegSites - 1 - i);
    for (unsigned int j = i + 1; j < numOfSegSites; ++j){
      unsigned int nDD = 0; //haplotype count of derived-derived
      unsigned int nDA = 0; //haplotype count of derived-ancestral
      for (unsigned int n = 0; n < sampleSegAlleles.size(); ++n){
        //count haplotypes
        if (sampleSegAlleles.at(n).at(i) && sampleSegAlleles.at(n).at(j)) ++nDD;
        else if (sampleSegAlleles.at(n).at(i) && !sampleSegAlleles.at(n).at(j)) ++nDA;
      }
      // check concistency between hf and maf
      if ((unsigned int)(mafs.at(i) * sampleSize) != (nDD + nDA)) std::cerr << "inconsistency between haplotype frequency and allele frequency" << std::endl;
      DDhfs_tmp.push_back((nDD*1.0)/sampleSize);
    }
    DDhf_vec.push_back(DDhfs_tmp);
  }
  return DDhf_vec;
}

void Samples::getLD(std::map<unsigned int, double>& r2Table, std::map<unsigned int, unsigned long>& r2Redundancy){
  if (numOfSegSites > 1){
    std::vector<double> mafs = getMAFs();
    std::vector<std::vector<double>> hfs = getHFs();
    for (unsigned int i = 0; i < (numOfSegSites - 1); ++i){//first locus
      for (unsigned int j = i + 1; j < numOfSegSites; ++j){//second locus
        unsigned int interval = segSite_positions.at(j) - segSite_positions.at(i);
        double r2 = std::pow((hfs.at(i).at(j - (i + 1)) - mafs.at(i) * mafs.at(j)), 2.0)/(mafs.at(i) * (1 - mafs.at(i)) * mafs.at(j) * (1 - mafs.at(j)));
        r2Table[interval] += r2;
        ++r2Redundancy[interval];
      }
    }
  }}
