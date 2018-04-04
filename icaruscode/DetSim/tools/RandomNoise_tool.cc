////////////////////////////////////////////////////////////////////////
/// \file   RandomNoise.cc
/// \author F. Varanini
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IGenNoise.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/Utilities/LArFFT.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include <fstream>

namespace icarus_tool
{

class RandomNoise : IGenNoise
{
public:
    explicit RandomNoise(const fhicl::ParameterSet& pset);
    
    ~RandomNoise();
    
    void configure(const fhicl::ParameterSet& pset)                          override;

    void GenerateNoise(std::vector<float> &noise, double noise_factor) const override;
    
private:
    // Member variables from the fhicl file
    std::string         fInputNoiseHistFileName;
    std::string         fHistogramName;
    
    // We'll recover the bin contents and store in a vector
    // with the likely false hope this will be faster...
    std::vector<double> fNoiseHistVec;
};
    
//----------------------------------------------------------------------
// Constructor.
RandomNoise::RandomNoise(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
RandomNoise::~RandomNoise()
{
}
    
void RandomNoise::configure(const fhicl::ParameterSet& pset)
{
   
    return;
}

void RandomNoise::GenerateNoise(std::vector<float> &noise, double noise_factor) const
{
    //ART random number service
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine&                        engine = rng->getEngine("noise");
    CLHEP::RandGaussQ                              rGauss(engine, 0.0, noise_factor);
    
    //In this case noise_factor is a value in ADC counts
    //It is going to be the Noise RMS
    //loop over all bins in "noise" vector
    //and insert random noise value
     for (unsigned int i=0; i<noise.size(); i++)
         noise.at(i) = rGauss.fire();

    return;
}
    
DEFINE_ART_CLASS_TOOL(RandomNoise)
}
