#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <functional>

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Grid.h>
#include <crpropa/GridTools.h>
#include <crpropa/Candidate.h>
#include <crpropa/Vector3.h>

#include "ppInteractions/Constants.h"
#include "ppInteractions/ParticleDecay.h"

#ifdef _OPENMP
	#include "omp.h"
#endif


/**
 @class NucleusNucleusInteraction
 @brief Handles pp interactions following the parametrisations from:
 
 Kelner et al. Phys. Rev. D 74 (2006) 034018.
 Kafexhiu et al. Phys. Rev. D 90 (2014) 123014.

 Note that this is a hybrid approach: Kafexhiu is used only to model the gamma-ray
 spectrum due to the pi-0 decay, whereas Kelner+ is used to model everything else.
 */
class NucleusNucleusInteraction : public crpropa::Module
{
public:
	std::vector<double> logFraction;
	std::vector<double> logIncidentEnergy;
	std::vector<double> logProbabilities;
	std::vector<double> neutralPionFraction;
	std::vector<double> chargedPionFraction;
	std::vector<double> etaMesonFraction;
	std::vector<double> probabilities;

	crpropa::ref_ptr<crpropa::Grid1f> densityGrid;

	bool haveElectrons;
	bool havePhotons;
	bool haveNeutrinos;
	bool havePions;
	bool haveMuons;
	bool doDecayChargedPion;
	bool doDecayNeutralPion;
	bool doDecayMuon;
	double limit;
	double thinning;
	double normMatterField;
	bool isDensityConstant;

	NucleusNucleusInteraction(double normMatterField = 1., double thinning = 0, double limit = 0.1);
	NucleusNucleusInteraction(crpropa::ref_ptr<crpropa::Grid1f> densityGrid, double normMatterField = 1., double thinning = 0, double limit = 0.1);
	void setLimit(double limit);
	void setThinning(double thinning);
	void setFieldNorm(double normMatterField);
	void setIsDensityConstant(bool densityConstant);
	void initMesonSpectra();
	double crossSection(double energy) const;
	double energyFractionNeutralPion(double energy, double xmin = 0, double xmax = 1) const;
	double energyFractionChargedPion(double energy, double xmin = 0, double xmax = 1) const;
	double energyFractionEtaMeson(double energy, double xmin = 0, double xmax = 1) const;
	double lossLength(int id, double E) const; 
	double lossLength(int id, double E, crpropa::Vector3d position) const; 
	void process(crpropa::Candidate *candidate) const;
	void performInteraction(crpropa::Candidate *candidate) const;
};







