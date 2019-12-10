#include "ProtonProtonInteraction.h"

using namespace crpropa;

ProtonProtonInteraction::ProtonProtonInteraction(double normBaryonField, bool photons, bool electrons, bool neutrinos, double thinning, double limit) : Module() {

	initSpectra();
	setFieldNorm(normBaryonField);
	setHaveElectrons(electrons);
	setHaveNeutrinos(neutrinos);
	setHavePhotons(photons);
	setLimit(limit);
	setIsDensityConstant(true);
	setDescription("ProtonProtonInteraction");
}

ProtonProtonInteraction::ProtonProtonInteraction(ref_ptr<ScalarGrid> grid, double normBaryonField, bool photons, bool electrons, bool neutrinos, double limit) : Module() {
	
	initSpectra();
	setFieldNorm(normBaryonField);
	setHaveElectrons(electrons);
	setHaveNeutrinos(neutrinos);
	setHavePhotons(photons);
	setLimit(limit);
	setDescription("ProtonProtonInteraction");
	setIsDensityConstant(false);
	if (normBaryonField != 1) 
		scaleGrid(grid, normBaryonField);
	densityGrid = grid;
}

void ProtonProtonInteraction::setIsDensityConstant(bool b) {
	isDensityConstant = b;
}

void ProtonProtonInteraction::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void ProtonProtonInteraction::setHavePhotons(bool b) {
	havePhotons = b;
}

void ProtonProtonInteraction::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void ProtonProtonInteraction::setFieldNorm(double x) {
	if (x == 0) 
		x = std::numeric_limits<double>::min();
	normBaryonField = x;
}

void ProtonProtonInteraction::setLimit(double l) {
	limit = l;
}

void ProtonProtonInteraction::setThinning(double thinning) {
	thinning = thinning;
}

void ProtonProtonInteraction::initSpectra() {
	positronSpec = new LeptonSpectrum(-11);
	muonNuSpec = new LeptonSpectrum(14);
	muonAntiNuSpec = new LeptonSpectrum(-14);
	electronNuSpec = new LeptonSpectrum(12);
	gammaSpec = new GammaSpectrum();
	pionSpec = new PionSpectrum();
	tabFracEnergy = pionSpec->energy;
	tabFracProb = pionSpec->prob;
	tabFrac = pionSpec->frac;

	// for (int i = 0; i < gammaSpec->prob.size(); i++){
	// 	std::cout << gammaSpec->prob[i] / eV << std::endl;//" " << gammaSpec->frac << " " << gammaSpec->prob << std::endl;
	// }

	/* TEST
  	std::ofstream myfile("test.txt");
  	int k = 0;
  	for (int i = 0; i < tabFracEnergy.size(); i++) {
  		for (int j = 0; j < tabFracProb.size(); j++) {
  			myfile << tabFracEnergy[i] / eV << "\t" << tabFracProb[j] << "\t" << tabFracEnergy[k] << std::endl;
  			k++;
  		}
  	}
	myfile.close();
	*/
}

void ProtonProtonInteraction::process(Candidate *candidate) const {

	// check if nucleus
	int id = candidate->current.getId();
	if (not (isNucleus(id)))
		return;

	// the loop should be processed at least once for limiting the next step
	double step = candidate->getCurrentStep();
	do {
		double E = candidate->current.getEnergy();
		double z = candidate->getRedshift();
		double rate;
		if (isDensityConstant == false) {
			Vector3d pos = candidate->current.getPosition();
			rate = 1 / lossLength(id, E, pos);
		} else
			rate = 1 / lossLength(id, E);

		// find interaction mode with minimum random decay distance
		Random &random = Random::instance();
		double randDistance = std::numeric_limits<double>::max();       
		double d = -log(random.rand()) / rate;
		if (d > randDistance)
			continue;
		randDistance = d;

		// check if interaction doesn't happen
		if (step < randDistance) {
			// limit next step to a fraction of the mean free path
			candidate->limitNextStep(limit / rate);
			return;
		}
		// interact and repeat with remaining step
		performInteraction(candidate);
		step -= randDistance;
	} while (step > 0);
}

void ProtonProtonInteraction::performInteraction(Candidate *candidate) const {
	
	double en = candidate->current.getEnergy();
	double w0 = candidate->getWeight();

	Random &random = Random::instance();

	// select position of secondary within step
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// check if in tabulated energy range
	if (en < tabFracEnergy.front() or (en > tabFracEnergy.back()))
		return;

	// Use Kafexhiu et al. 2014 approach to compute pion multiplicity
	// Use Kelner et al. 2006 to scale lepton spectra wrt gamma rays.
	int mult_g = pionSpec->neutralPionMultiplicity(en);
	int mult_pos = positronSpec->computeMultiplicity(pionSpec, mult_g);
	int mult_nuEl = electronNuSpec->computeMultiplicity(pionSpec, mult_g);
	int mult_nuMu = muonNuSpec->computeMultiplicity(pionSpec, mult_g);
	int mult_nuAMu = muonAntiNuSpec->computeMultiplicity(pionSpec, mult_g);

	// if (en / eV > 1e13)
	// std::cout << en / eV << " " << mult_g << " " << mult_pos << " " << mult_nuMu <<  " " << mult_nuAMu << " " << mult_nuEl << std::endl;

	double fmin = 1e-6;
	double fmax = 1;
	int i = 0;
	double xtot = 0;
	double dEtot_g = 0;
	while(i < mult_g and xtot < 1 and fmax > fmin) {
		// double x_pi = interpolate2d(en, random.randUniform(fmin, fmax), tabFracEnergy, tabFracProb, tabFrac);
		double x_pi = interpolate2d(en, random.rand(), tabFracEnergy, tabFracProb, tabFrac);
		double en_pi = x_pi * en;
		double x0 = interpolate2d(en, random.rand(), gammaSpec->energy, gammaSpec->prob, gammaSpec->frac);
		if (havePhotons) {
			double f = x0 * en_pi / en;
			if (random.rand() < pow(f, thinning)) {
				double w = w0 / pow(f, thinning);
				candidate->addSecondary(22, x0 * en_pi, pos, w);
			}
			if (random.rand() < pow(1 - f, thinning)) {
				double w = w0 / pow(1 - f, thinning);
				candidate->addSecondary(22, (1 - x0) * en_pi, pos, w);
			}			
		}
		fmax -= x_pi;
		xtot += x_pi;
		dEtot_g += x0 * en_pi;
		i++;
	} 
	xtot = 0;
	fmax = 1;
	i = 0;
	double dEtot_pos = 0;
	while(i < mult_pos and xtot < 1 and fmax > fmin) {
		// double x_pi = interpolate2d(en, random.randUniform(fmin, fmax), tabFracEnergy, tabFracProb, tabFrac);
		double x_pi = interpolate2d(en, random.rand(), tabFracEnergy, tabFracProb, tabFrac);
		double en_pi = x_pi * en;
		double x0 = interpolate(random.randUniform(fmin, fmax), positronSpec->prob, positronSpec->frac);
		if (haveElectrons) {
			double f = x0 * en_pi / en;
			if (random.rand() < pow(f, thinning)) {
				double w = w0 / pow(f, thinning);
				candidate->addSecondary(-11, x0 * en_pi, pos, w);
			}
		}
		fmax -= x_pi;
		xtot += x_pi;
		dEtot_pos += x0 * en_pi;
		i++;
	} 
	xtot = 0;
	fmax = 1;
	i = 0;
	double dEtot_nuMu = 0;
	while(i < mult_nuMu and xtot < 1 and fmax > fmin) {
		// double x_pi = interpolate2d(en, random.randUniform(fmin, fmax), tabFracEnergy, tabFracProb, tabFrac);
		double x_pi = interpolate2d(en, random.rand(), tabFracEnergy, tabFracProb, tabFrac);
		double en_pi = x_pi * en;
		double x0 = interpolate(random.randUniform(fmin, fmax), muonNuSpec->prob, muonNuSpec->frac);
		if (haveNeutrinos) {
			double f = x0 * en_pi / en;
			if (random.rand() < pow(f, thinning)) {
				double w = w0 / pow(f, thinning);
				candidate->addSecondary(14, x0 * en_pi, pos, w);
			}
		}
		fmax -= x_pi;
		xtot += x_pi;
		dEtot_nuMu += x0 * en_pi;
		i++;
	} 
	xtot = 0;
	fmax = 1;
	i = 0;
	double dEtot_nuAMu = 0;
	while(i < mult_nuAMu and xtot < 1 and fmax > fmin) {
		// double x_pi = interpolate2d(en, random.randUniform(fmin, fmax), tabFracEnergy, tabFracProb, tabFrac);
		double x_pi = interpolate2d(en, random.rand(), tabFracEnergy, tabFracProb, tabFrac);
		double en_pi = x_pi * en;
		double x0 = interpolate(random.randUniform(fmin, fmax), muonAntiNuSpec->prob, muonAntiNuSpec->frac);
		if (haveNeutrinos) {
			double f = x0 * en_pi / en;
			if (random.rand() < pow(f, thinning)) {
				double w = w0 / pow(f, thinning);
				candidate->addSecondary(-14, x0 * en_pi, pos, w);
			}
		}
		fmax -= x_pi;
		xtot += x_pi;
		dEtot_nuMu += x0 * en_pi;
		i++;
	} 
	xtot = 0;
	fmax = 1;
	i = 0;
	double dEtot_nuEl = 0;
	while(i < mult_nuEl and xtot < 1 and fmax > fmin) {
		// double x_pi = interpolate2d(en, random.randUniform(fmin, fmax), tabFracEnergy, tabFracProb, tabFrac);
		double x_pi = interpolate2d(en, random.rand(), tabFracEnergy, tabFracProb, tabFrac);
		double en_pi = x_pi * en;
		double x0 = interpolate(random.randUniform(fmin, fmax), electronNuSpec->prob, electronNuSpec->frac);
		if (haveNeutrinos) {
			double f = x0 * en_pi / en;
			if (random.rand() < pow(f, thinning)) {
				double w = w0 / pow(f, thinning);
				candidate->addSecondary(12, x0 * en_pi, pos, w);
			}
		}
		fmax -= x_pi;
		xtot += x_pi;
		dEtot_nuEl += x0 * en_pi;
		i++;
	} 

	double dEtot = dEtot_g - dEtot_pos - dEtot_nuMu - dEtot_nuAMu - dEtot_nuEl;
	double newE = std::min(en - dEtot, en); // redundant way to prevent energy gain


	// Guarantees that remaining proton has a lorentz factor of at least 10.
	// Note that this is redundant with the module MinimumEnergy, but prevents problems.
	if (newE < 1e10 * eV) {
		candidate->setActive(false);
		return;
	}

	candidate->current.setEnergy(newE);
	// candidate->setActive(false);
}

double ProtonProtonInteraction::lossLength(int id, double E) const {
	double rate = crossSection(E) * normBaryonField;
	return 1. / rate;
}

double ProtonProtonInteraction::lossLength(int id, double E, Vector3d position) const {
	double rate = crossSection(E) * densityGrid->interpolate(position);
	return 1. / rate;
}

double ProtonProtonInteraction::crossSection(double en) const {
	// Parametrisation from Kafexhiu et al. 2014
	// Note: only works for en >> Ethr
	// values are given in mb, hence the 1e-31 factor
	double ethr = 2.797e8 * eV;
	double x = (en - mass_proton * c_squared) / ethr;
    double res = (30.7 - 0.96 * log(x) + 0.18 * log(x) * log(x)) * pow(1. - pow(1 / x, 1.9), 3);
    return res * 1e-31; 
}

LeptonSpectrum::LeptonSpectrum() {
}

LeptonSpectrum::LeptonSpectrum(int id) {
	setLepton(id);
	init();
}

void LeptonSpectrum::init() {
	prob.clear();
	ratio.clear();
	frac.clear();
	switch(leptonId) {
		case -11:
			positronDistribution();
			break;
		case 12:
			electronNeutrinoDistribution();
			break;
		case 14:
			muonNeutrinoDistribution();
			break;
		case -14:
			muonAntiNeutrinoDistribution();
			break;
		default:
			throw std::invalid_argument("Unknown id for lepton.");
	}	
}

void LeptonSpectrum::setLepton(int id) {
	leptonId = id;
}

void LeptonSpectrum::setFraction(std::vector<double> fraction) {
	frac = fraction;
}

void LeptonSpectrum::setRatio(std::vector<double> ratio) {
	ratio = ratio;
}

void LeptonSpectrum::setProbability(std::vector<double> probability) {
	prob = probability;
}

int LeptonSpectrum::computeMultiplicity(PionSpectrum *ps, int mult0) const {
	// Takes the gamma-ray multiplicity and used the scaling relations by
	// Kelner et al. to obtain the multiplicity of the other particles.
	// Gamma-ray multiplicities are from Kafexhiu et al. 2014.
	double m = (double) mult0;
	int nxspec = 10;
	int nintegral = 100;
	double r = 0;
	if (leptonId == 14) {
		double l = 0.427;
		double res = 0;
		for (int i = 0; i < nxspec; i++) {
			double xmin = (double) i / nxspec;
			double xmax = (double) (i + 1) / nxspec;
			double a = ps->computeSlopeInInterval(xmin, xmax);
			r += pow(l, a);
		}
		r = (double) r / nxspec;
	} else {
		double r = 0;
		double m = (double) mult0;
		int nxspec = 10;
		int nintegral = 100;
		for (int i = 0; i < nxspec; i++) {
			double xmin = (double) i / nxspec;
			double xmax = (double) (i + 1) / nxspec;
			double a = ps->computeSlopeInInterval(xmin, xmax);
			double res = 0;
			for (int j = 0; j < nintegral; j++) {
				double dx = (double) (xmax - xmin) / nintegral;
				double x = (double) (2 * j + 1) * dx / 2;
				double x0 = interpolate(x, frac, ratio);
				res = res + x0 * pow(x, a - 1) * dx;
			}
			r += a * res;
		}
	}
	Random &random = Random::instance();
	if (random.rand() > 0.5) 	// correction for floor function.
		return std::floor(r * m) + 1; 
	else
		return std::floor(r * m);
} 

double LeptonSpectrum::energyFraction(double pmin, double pmax) const {
	Random &random = Random::instance();
	return interpolate(random.randUniform(pmin, pmax), prob, frac);
}

double LeptonSpectrum::ratioToPhoton(double x) const {
	return interpolate(x, frac, ratio);
}

void LeptonSpectrum::muonNeutrinoDistribution() {
	// Muon neutrino due to pion+ decay (pi+ -> mu+ + nu_mu)
	// See Eqs. 23-26 from Kelner et al. 2006.
	double r = pow(mMuon / mChargedPion, 2);
	double res = 0;
	for (int i = 0; i < nsamples; i++) {
		double x = (double) (i + 1) / nsamples;
		frac.push_back(x);
		if (x <= r) {
			res = res + 2.367;
			prob.push_back(res);
			ratio.push_back(2.367);
		} else {
			prob.push_back(res);
			ratio.push_back(0.);
		}
	}
	for (int i = 0; i < prob.size(); i++)
		prob[i] /= res;
}

void LeptonSpectrum::positronDistribution() { 
	// Electron spectrum; follows Eq. 36 from Kelner et al. 2006.
	double r = pow(mMuon / mChargedPion, 2);
	double res = 0;
	for (int i = 0; i < nsamples; i++) {
		double x = (double) (i + 1) / nsamples;
		double g_1 = (3 - 2 * r) / (9 * (1 - r) * (1 - r));
		double g_2 = 9 * x * x - 6 * log(x) - 4 * x * x * x - 5;
		double g = g_1 * g_2;
		double h1_1 = (3 - 2 * r) / (9 * (1 - r) * (1 - r));
		double h1_2 = 9 * r * r - 6 * log(r) - 4 * r * r * r - 5;
		double h1 = h1_1 * h1_2;
		double h2_1 = (1 + 2 * r) * (r - x) / (9 * r * r);
		double h2_2 = 9 * (r + x) - 4 * (r * r + r * x + x * x);
		double h2 = h2_1 * h2_2;
		if (x <= r) {
			res = res + h1 + h2;
			prob.push_back(res);
			ratio.push_back(h1 + h2);
		} else {
			res += g;
			prob.push_back(res);
			ratio.push_back(g);
		}
		frac.push_back(x);
	}
	for (int i = 0; i < prob.size(); i++)
		prob[i] /= res;
}

void LeptonSpectrum::electronNeutrinoDistribution() {
	// Electron neutrino spectrum; follows Eq. 40 from Kelner et al. 2006.
	// Note that there is a typo in this reference. The correct version is presented
	// in the erratum: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.79.039901
	double r = pow(mMuon / mChargedPion, 2);
	double res = 0;
	for (int i = 0; i < nsamples; i++) {
		double x = (double) (i + 1) / nsamples;
		double g_1 = 2 / (3 * (1 - r) * (1 - r)); 
		double g_2 = (6 * (1 - x) * (1 - x) + r * (5 + 5 * x - 4 * x * x)) * (1 - x) + 6 * r * log(x);
		double g = g_1 * g_2;
		double h1_1 = 2 / (3 * (1 - r) * (1 - r));
		double h1_2 = (1 - r) * (6 - 7 * r + 11 * r * r - 4 * r * r * r) + 6 * r * log(r);
		double h1 = h1_1 * h1_2;
		double h2_1 = 2 * (r - x) / (3 * r * r);
		double h2_2 = (7 * r * r - 4 * r * r * r + 7 * x * r - 4 * x * r * r - 2 * x * x - 4 * r * x * x);
		double h2 = h2_1 * h2_2;
		if (x >= r) {
			res = res + h1 + h2;
			prob.push_back(res);
			ratio.push_back(h1 + h2);
		} else {
			res += g;
			prob.push_back(res);
			ratio.push_back(g);
		}
		frac.push_back(x);
	}
	for (int i = 0; i < prob.size(); i++)
		prob[i] /= res;
}

void LeptonSpectrum::muonAntiNeutrinoDistribution() {
	// Initiates the table computed following Eq. 36 from Kelner et al. 2006.
	positronDistribution();
}

PionSpectrum::PionSpectrum() {
	initSpectrum();

	// std::cout << "pi frac: " << frac.size() << "  " << frac[0] << " " << frac[frac.size()-1] << std::endl;
	// std::cout << "pi prob: " << prob.size() << "  " << prob[0] << " " << prob[prob.size()-1] << std::endl;
	// std::cout << "pi ener: " << energy.size() << "  " << energy[0] / eV << " " << energy[energy.size()-1] / eV << std::endl;
}

void PionSpectrum::initSpectrum() {
	// Initiates the table computed following Eq. 11 of Kafexhiu et al. 2014.

	// clear previously loaded tables
	frac.clear();
	prob.clear();
	energy.clear();

	const static double mp = 938.27231 * 1e6; 
	const static double mpi = 134.9767 * 1e6;
	const int nx = 600;
	const int ne = 110; 
	const int np = 600;
	double frac_tmp[nx][ne];
	double prob_tmp[nx];
	double energy_tmp[ne];
	const double emin = 1e10 * eV;
	const double emax = 1e21 * eV;
	const double fmin = 1e-6;
	const double fmax = 1;
	std::vector<double> f_tmp;

	for (int i = 0; i < np; i++) {
		double dp = log10(fmax / fmin) / np;
		double p = log10(fmin) + (double) (i + 1) * dp;
		p = pow(10, p);
		prob.push_back(p);
	}

	for (int j = 0; j < ne; j++ ) {	
		double de = log10(emax / emin) / ne;
		double Ep = log10(emin) + (double) (j + 1) * de;
		Ep = pow(10, Ep);
		energy.push_back(Ep);
		double p = 0;
		double pmax = 0;
		std::vector<double> p_tmp;
		std::vector<double> x_tmp;
		for (int i = 0; i < nx; i++) {
			double dx = log10(fmax / fmin) / nx;
			double x = log10(fmin) + (double) (i + 1) * dx;
			x = pow(10, x);
		    double L = log((Ep / eV - mp) / 1e12);
		    double Bpi = 5.58 + 0.78 * L + 0.1 * L * L;
		    double r = 3.1 / pow(Bpi, 1.5);
		    double alpha = 0.89 / (pow(Bpi, 0.5) * (1 - exp(-0.33 * Bpi)));
		    double f1 = 4 * alpha * Bpi * pow(x, alpha - 1);
		    double f2 = pow(1 - pow(x, alpha), 4) / pow(pow(1 + r * pow(x, alpha), 3), 4);
		    double f3 = 1 / (1 - pow(x, alpha)) + 3 * r / (1 + r * pow(x, alpha));
		    double f4 = sqrt(1 - mpi / (x * Ep / eV));
		    double f = p + f1 * f2 * f3 * f4;
		    if (f < std::numeric_limits<double>::max() and f > std::numeric_limits<double>::min())
			    p += f1 * f2 * f3 * f4;
		    x_tmp.push_back(x);
		    p_tmp.push_back(p);
		    if (i == nx - 1)
		    	pmax = p;
		}

		for (int k = 0; k < p_tmp.size(); k++) {
			p_tmp[k] /= pmax;
		}

		for (int k = 0; k < prob.size(); k++) {
			double y = interpolate(prob[k], p_tmp, x_tmp);	
			frac.push_back(y);
		}
		
	}
}

double PionSpectrum::energyFraction(double en) const {
	Random &random = Random::instance();
	return interpolate2d(random.rand(), en, prob, energy, frac);
}

int PionSpectrum::neutralPionMultiplicity(double en) const {
	// Parametrisation from Kafexhiu et al. 2014 for the pi0 multiplicity
	// Experimental data up to 5 GeV, GEANT 4.10.0 up to 100 TeV, and SYBILL 2.1 above.
	double Tp = en - mass_proton * c_squared;
	if (Tp > 1e9 * eV and Tp < 5e9 * eV) { 
		double q = (Tp - 2.797e8 * eV) / (mass_proton * c_squared);
		double res = -6e-3 + 0.237 * q - 0.023 * q * q;
		return (int) std::floor(res);
	} else if (Tp >= 5e9 * eV and Tp < 1e14 * eV) { 
		double a1 = 0.728;
		double a2 = 0.596;
		double a3 = 0.491;
		double a4 = 0.2503;
		double a5 = 0.117;
		double csi = (Tp - 3e9 * eV) / (mass_proton * c_squared);
		double res = a1 * pow(csi, a4) * (1 + exp(-a2 * pow(csi, a5))) * (1 - exp(-a3 * pow(csi, 0.25)));
		return (int) std::floor(res);
	} else {
		double a1 = 5.436;
		double a2 = 0.254;
		double a3 = 0.072;
		double a4 = 0.075;
		double a5 = 0.166;
		double csi = (Tp - 3e9 * eV) / (mass_proton * c_squared);
		double res = a1 * pow(csi, a4) * (1 + exp(-a2 * pow(csi, a5))) * (1 - exp(-a3 * pow(csi, 0.25)));
		return (int) std::floor(res);
	}
}

double PionSpectrum::computeSlopeInInterval(double xmin, double xmax) const {
	// In future releases this will be changed to read from a tabulated pion spectrum
	// Note that this value requires simulations to be run with spectrum E^-1
	return 1;

	// // minimum
	// double x = xmin;
	// double L = log((en / eV - mp) / 1e12);
	// double Bpi = 5.58 + 0.78 * L + 0.1 * L * L;
	// double r = 3.1 / pow(Bpi, 1.5);
	// double alpha = 0.89 / (pow(Bpi, 0.5) * (1 - exp(-0.33 * Bpi)));
	// double f1 = 4 * alpha * Bpi * pow(x, alpha - 1);
	// double f2 = pow(1 - pow(x, alpha), 4) / pow(pow(1 + r * pow(x, alpha), 3), 4);
	// double f3 = 1 / (1 - pow(x, alpha)) + 3 * r / (1 + r * pow(x, alpha));
	// double f4 = sqrt(1 - mpi / (x * Ep / eV));
	// double fmin = f1 * f2 * f3 * f4;

	// double x = xmax;
	// L = log((en / eV - mp) / 1e12);
	// Bpi = 5.58 + 0.78 * L + 0.1 * L * L;
	// r = 3.1 / pow(Bpi, 1.5);
	// alpha = 0.89 / (pow(Bpi, 0.5) * (1 - exp(-0.33 * Bpi)));
	// f1 = 4 * alpha * Bpi * pow(x, alpha - 1);
	// f2 = pow(1 - pow(x, alpha), 4) / pow(pow(1 + r * pow(x, alpha), 3), 4);
	// f3 = 1 / (1 - pow(x, alpha)) + 3 * r / (1 + r * pow(x, alpha));
	// f4 = sqrt(1 - mpi / (x * Ep / eV));
	// fmax = f1 * f2 * f3 * f4;

	// double slope = (fmax - fmin) / (xmax - xmin);
	// std::cout << slope << std::endl;
	// return slope;
}

int PionSpectrum::chargedPionMultiplicity(double en) const {
	// Parametrisation from Kelner et al. 2006 (Eq. 4).
	// Parameters are for QGSJet.
	double x = en / mass_proton;
	double L = log(en / (1e12 * eV));
	double Bpi = 5.58 + 0.78 * L + 0.10 * L * L;
	double r = 3.1 / pow(Bpi, 1.5);
	double alpha = 0.89 / (sqrt(Bpi) * (1 - exp(-0.33 * Bpi)));
	double a = 1 - pow(x, alpha);
	double b = pow(1 + r * pow(x, alpha), 3);
	// return (int) std::floor(Bpi * pow(a / b, 4));
	double pf = 4 * alpha * Bpi * pow(x, alpha - 1);
	double f1 = pow((1 - pow(x, alpha)) / pow(1 + r * pow(x, alpha), 3), 4);
	double f2 = 1 / (1 - pow(x, alpha)) + 3 * r / (1 + r * pow(x, alpha));
	double f3 = sqrt(1 - mChargedPion / (x * en));
	return (int) std::floor(pf * f1 *f2 * f3 * x);
}

double PionSpectrum::crossSection(double en) const {
	double Tp = en - mass_proton * c_squared;
	double Tpthr = 2.797e8 * eV;
	if (Tp >= Tpthr and Tp < 2e9 * eV) {
		return onePionCrossSection(en) + twoPionCrossSection(en);
	} else {
		double x = Tp / Tpthr;
	    double piCS = (30.7 - 0.96 * log(x) + 0.18 * log(x) * log(x)) * pow(1. - pow(1 / x, 1.9), 3);
	    return 1e-31 * piCS * neutralPionMultiplicity(en);
	}
}

double PionSpectrum::onePionCrossSection(double en) const {
	double s = centreOfMassEnergySquared(en);
	double sigma0 = 7.66e-3 * 1e-31;
	double mpi = 134.9766 * 1e6 * eV;
	double mp = mass_proton * c_squared;
	double eta = sqrt(pow(s - mpi * mpi * 4 * mp * mp, 2) - 16 * mpi * mpi * mp * mp) / (2 * mpi * sqrt(s));
	return sigma0 * pow(eta, 1.95) * (1 + eta + pow(eta, 5)) * pow(breitWigner(s), 1.86);
}

double PionSpectrum::twoPionCrossSection(double en) const {
	double Tp = en - mass_proton * c_squared;
	return 5.7e-31 / (1 + exp(-9.3 * (Tp - 1.4e9 * eV)));
}

double PionSpectrum::breitWigner(double s) const {
	double mp = mass_proton * c_squared;
	double mres = 1.1883e9 * eV;
	double gammares = 0.2264e9 * eV;
	double gamma = sqrt(mres * mres * (mres * mres + gammares * gammares));
	double k = sqrt(8) * mres * gammares * gamma / (M_PI * sqrt(mres * mres + gamma));
	double f1 = pow(sqrt(s) - mp, 2) - mres * mres;
	double f2 = mres * gammares;
	return mp * k / (f1 * f1 + f2 * f2);
}

double PionSpectrum::centreOfMassEnergySquared(double en) const {
	double Tp = en - mass_proton * c_squared;
	return 2 * mass_proton * c_squared * (Tp + 2 * mass_proton * c_squared);
}

double PionSpectrum::pionMaximumEnergyLab(double en) const {
	double Tp = en - mass_proton * c_squared;
	double mp = mass_proton * c_squared;
	double mpi = 134.9766 * 1e6 * eV;
	double s = centreOfMassEnergySquared(en);
	double gammaCM = lorentzFactorCMF(en);
	double betaCM = sqrt(1 - 1 / (gammaCM * gammaCM));
	double EpiCM = pionEnergyCMF(en);
	double PpiCM = sqrt(EpiCM * EpiCM - mpi * mpi);
	return gammaCM * (EpiCM + PpiCM * betaCM);
}

double PionSpectrum::pionEnergyCMF(double en) const {
	double s = centreOfMassEnergySquared(en);
	double mp = mass_proton * c_squared;
	double mpi = 134.9766 * 1e6 * eV;
	return (s - 4 * mp * mp + mpi * mpi) / (2 * sqrt(s));
}

double PionSpectrum::lorentzFactorCMF(double en) const {
	double Tp = en - mass_proton * c_squared;
	double mp = mass_proton * c_squared;
	double mpi = 134.9766 * 1e6 * eV;
	double s = centreOfMassEnergySquared(en);
	return (Tp + 2 * mp) / sqrt(s);
}

double PionSpectrum::lorentzFactorLab(double en) const {
	double mpi = 134.9766 * 1e6 * eV;
	return pionMaximumEnergyLab(en) / mpi;
}

GammaSpectrum::GammaSpectrum() {
	pionSpec = new PionSpectrum();
	initSpectrum();
	// std::cout << "g frac: " << frac.size() << "  " << frac[0] << " " << frac[frac.size()-1] << std::endl;
	// std::cout << "g prob: " << prob.size() << "  " << prob[0] << " " << prob[prob.size()-1] << std::endl;
	// std::cout << "g ener: " << energy.size() << "  " << energy[0] / eV << " " << energy[energy.size()-1] / eV << std::endl;

	// for (int i = 0; i < frac.size(); i++)
	// 	std::cout << frac[i] << std::endl;
}

void GammaSpectrum::initSpectrum() {
	// Initiates the table computed following Eq. 8 of Kafexhiu et al. 2014.
	// clear previously loaded tables
	frac.clear();
	prob.clear();
	energy.clear();

	const int nx = 600;
	const int ne = 110; 
	const int np = 600;
	double frac_tmp[nx][ne];
	double prob_tmp[nx];
	double energy_tmp[ne];
	const double emin = 1e10 * eV;
	const double emax = 1e21 * eV;
	const double fmin = 1e-6;
	const double fmax = 1;
	std::vector<double> f_tmp;

	for (int i = 0; i < np; i++) {
		double dp = log10(fmax / fmin) / np;
		double p = log10(fmin) + (double) (i + 1) * dp;
		p = pow(10, p);
		prob.push_back(p);
	}

	for (int j = 0; j < ne; j++ ) {	
		double de = log10(emax / emin) / ne;
		double Ep = log10(emin) + (double) (j + 1) * de;
		Ep = pow(10, Ep);
		energy.push_back(Ep);
		double p = 0;
		double pmax = 0;
		std::vector<double> p_tmp;
		std::vector<double> x_tmp;
		for (int i = 0; i < nx; i++) {
			double dx = log10(fmax / fmin) / nx;
			double x = log10(fmin) + (double) (i + 1) * dx;
			x = pow(10, x);
		    double f = differentialCrossSection(Ep, x * Ep) * x;
		    if (f < std::numeric_limits<double>::max() and f > std::numeric_limits<double>::min())
			    p += f;
		    // std::cout << x << " " << f << " " << p << std::endl;
		    x_tmp.push_back(x);
		    p_tmp.push_back(p);
		    if (i == nx - 1)
		    	pmax = p;
		}

		for (int k = 0; k < p_tmp.size(); k++) {
			p_tmp[k] /= pmax;
		}

		for (int k = 0; k < prob.size(); k++) {
			double y = interpolate(prob[k], p_tmp, x_tmp);	
			frac.push_back(y);
			// std::cout << Ep / eV << " " << y << " " << prob[k] << std::endl;
		}
	}
}

double GammaSpectrum::differentialCrossSection(double en, double Eg) const {
	// std::cout << en / eV << " " << Eg / eV << " " << Fdist(en, Eg) << std::endl;
	return Amax(en) * Fdist(en, Eg);
}

double GammaSpectrum::Amax(double en) const {
	// GEANT 4.10.0 below 100 TeV, SYBILL 2.1 above
	// See Eq. 12 Kafexhiu et al. 2014, and Table VII.

	double Tp = en - mass_proton * c_squared;
	double Tpthr = 2.797e8 * eV;

	double b0, b1, b2, b3;
	if (Tp >= 1e9 * eV and Tp < 5e9 * eV) {
		b1 = 9.53;
		b2 = 0.52;
		b3 = 0.054;
	} else if (Tp >= 5e9 * eV and Tp < 1e14 * eV) {
		b1 = 9.13;
		b2 = 0.35;
		b3 = 9.7e-3;
	} else {
		b1 = 10.77;
		b2 = 0.412;
		b3 = 0.01264;
	}
	double Epimax = pionSpec->pionMaximumEnergyLab(en);

	double Amax = 0;
	if (Tp >= Tpthr and Tp < 1e9 * eV) {
		double b0 = 5.9;
		Amax = b0 * pionSpec->crossSection(en) / Epimax;
	} else {
		double t = Tp / (mass_proton * c_squared);
		Amax = b1 * pow(t, -b2) * exp(b3 * log(t) * log(t)) * pionSpec->crossSection(en) / (mass_proton * c_squared);
	}	
	// cheating to fix differential cross sections, which are wrong by a log(10)^2 factor.
	// note that this won't affect the results, since dSigma/dE_gamma won't be directly used, only when normalised.
	return Amax / pow(log(10), 2);
}

double GammaSpectrum::Fdist(double en, double Eg) const {
	// GEANT 4.10.0 below 100 TeV, SYBILL 2.1 above
	// See Eq. 11 Kafexhiu et al. 2014, and Table V.
	double mp = mass_proton * c_squared;
	double mpi = 134.9766 * 1e6 * eV;
	double lfCMF = pionSpec->lorentzFactorCMF(en);
	double lfLab = pionSpec->lorentzFactorLab(en);
	double betaLab = sqrt(1 - 1 / (lfLab * lfLab));
	double Egmax = 0.5 * mpi * lfLab * (1 + betaLab);
	double Yg = Eg + mpi * mpi / (4 * Eg);
	double Ygmax = Egmax + mpi * mpi / (4 * Egmax);
	double Xg = (Yg - mpi) / (Ygmax - mpi);
	
	double Tp = en - mass_proton * c_squared;
	double Tpthr = 2.797e8 * eV;
	double t = Tp / (mass_proton * c_squared);
	double kappa = 3.29 - 0.2 * pow(t, -1.5);
	double q = (Tp - 1e9 * eV) / mp;
	double mu = 5 / 4 * pow(q, 5 / 4) * exp(- 5 / 4 * q);
	double alpha, beta, gamma, lambda;
	if (Tp >= Tpthr and Tp < 1e9 * eV) {
		lambda = 1;
		alpha = 1.;
		beta = kappa;
		gamma = 0;
	} else if (Tp >= 1e9 * eV and Tp < 4e9 * eV) {
		lambda = 3;
		alpha = 1.;
		beta = mu + 2.45;
		gamma = mu + 1.45;
	} else if (Tp >= 4e9 * eV and Tp < 2e10 * eV) {
		lambda = 3;
		alpha = 1.;
		beta = 1.5 * mu + 4.95;
		gamma = mu + 1.50;
	} else if (Tp >= 2e10 * eV and Tp < 1e14 * eV) {
		lambda = 3;
		alpha = 0.5;
		beta = 4.2;
		gamma = 1.;
	} else {
		lambda = 3.55;
		alpha = 0.5;
		beta = 3.6;
		gamma = 1;
	}
	double C = lambda * mpi / Ygmax;
	return pow(1 - pow(Xg, alpha), beta) / pow(1 + Xg / C, gamma);
}
