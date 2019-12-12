#include "Individual.hpp"
#include "Population.hpp"

//Constructor
Individual::Individual( const int id_input,
                        const int numOfLoci_input,
                        const bool allele_input){
    id = id_input;
    alleles.reserve(numOfLoci_input);
    for (int i = 1; i <= numOfLoci_input; ++i){
      alleles.push_back(allele_input);
    }
}
//Function
