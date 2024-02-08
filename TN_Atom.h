/*******************************************
 *   Copyright (C) 2022 by Ben F McLean    *
 *   drbenmclean@gmail.com                 *
 *                                         *
 ******************************************/

#include <vector>
#include <random>
#include "TN_Math.h"

#ifndef TN_ATOM
#define TN_ATOM

/*
TN_Atom.h defines four classes:

1. A base class of TN_AtomicParticle. All atomic particles, whether atomic nucleii or electrons,
have a name and a location, mass, charge, velocity and acceleration, and potentially a force acting
upon it. Therefore it makes sense to have a base class defining all these properties, and derive
specific atomic particle classes from it. Further, knowing mass and charge allows for the calculation
of between-particle forces, hence the forceBetweenParticles function can be a member function here.

2. Deriving from TN_AtomicParticle, nucleii are defined by the addition of merely 3 member
variables - the number of inner-shell electrons, the number of protons, and the number of neutrons.

3. Electrons are also TN_AtomicParticles, that in addition have a spin, and a radius.

4. Finally, a TN_Atom class is declared. Just like physical atoms, a TN_Atom always has one TN_Nucleus,
and can have any number of TN_Electrons stored in a vector.

TN_Atom.h also contains the declaration of the all-important function forceBetweenParticles(p1, p2).
This function allows for the calculation of forces on p1 due to interaction with p2.

*/

extern double g_CoulombConst;

class TN_AtomicParticle {
	
	public:
	TN_Point m_location;
	TN_Vector m_velocity;
	TN_Vector m_force;
	TN_Vector m_acceleration;
	double m_mass;
	double m_charge;

	TN_AtomicParticle();
	virtual ~TN_AtomicParticle(); //include a virtual fn to enable dynamic casting
};


class TN_Nucleus : public TN_AtomicParticle {
	
	public:
	int m_nInnerShellElectrons = 0;
	int m_nProtons;
	int m_nNeutrons;

	TN_Nucleus();
	TN_Nucleus(int nNeutrons, int nProtons, double cx=0.0,  double cy=0.0, double cz=0.0);

	void addInnerShellElectrons(int n);

};


class TN_Electron : public TN_AtomicParticle {
	
	public:
	const TN_Nucleus* m_nucleus; //pointer to nucleus with which electron is associated
	double m_radius;
	double m_spin;

	TN_Electron(const TN_Nucleus& nucleus, double radius, double spin, double cx=0.0,  double cy=0.0, double cz=0.0);

};


class TN_Atom{
	
	public:
	TN_Nucleus m_nucleus; //note atom has no coords of its own, uses nucleus coords instead.
	std::vector<TN_Electron> m_electrons;
	
	TN_Atom(int nNeutrons, int nProtons, double cx=0.0,  double cy=0.0, double cz=0.0);
	
	void addElectron(double radius, double spin, double cx=0.0,  double cy=0.0, double cz=0.0);
	const TN_Point& location() const;
	const double& nucleusCharge() const;
	void addInnerShellElectrons(int n);
	double charge();
	void collapseElectronToRadius(int e);
	
};


//declaration of non-member function that calculates the force between any pair of atomic particles.
void forceBetweenParticles(TN_AtomicParticle &a, TN_AtomicParticle &b);



#endif