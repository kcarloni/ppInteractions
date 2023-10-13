#include "ppInteractions/DecayChargedPion.h"

using namespace crpropa;


DecayChargedPion::DecayChargedPion(bool muons, bool neutrinos, double thinning, double limit) {
	setHaveMuons(muons);
	setHaveNeutrinos(neutrinos);
	setThinning(thinning);
	setLimit(limit);
}

void DecayChargedPion::setHaveMuons(bool b) {
	haveMuons = b;
}

void DecayChargedPion::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void DecayChargedPion::setLimit(double l) {
	limit = l;
}

void DecayChargedPion::setThinning(double thinning) {
	thinning = thinning;
}

double DecayChargedPion::lossLength(const double& lf) const {
	// Returns the loss length in the lab frame.
	const static double lifetime = tauChargedPion;
	return c_light * lifetime * lf;
}

double DecayChargedPion::energyFractionMuon() const {

	// // approximation: beta = 1 and the rest frame opening angle always = 1. (not realistic)
	// return pow_integer<2>(mMuon / mChargedPion);

	// // approximation: beta = 1:
	// sample random rest frame angle:
	Random &random = Random::instance();
	double x = random.randUniform(0., 1.);
	return 0.5 * ( (1-x) + (1+x) * pow_integer<2>(mMuon/mChargedPion) );
}

void DecayChargedPion::performInteraction(Candidate* candidate) const {

	int id = candidate->current.getId();
	double sign = (id > 0) ? 1 : ((id < 0) ? -1 : 0);
	double E = candidate->current.getEnergy();  

	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// particle disappears after decay
	candidate->setActive(false);

	double f = energyFractionMuon();
	if (haveNeutrinos) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			candidate->addSecondary(sign * 13, E * (1 - f), pos, w);
		} 
	}
	if (haveMuons) {
		if (random.rand() < pow(f, thinning))  {
			double w = 1. / pow(f, thinning);
			candidate->addSecondary(-sign * 14, E * f, pos, w);
		} 
	}
}
