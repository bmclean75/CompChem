/*******************************************
 *   Copyright (C) 2022 by Ben F McLean    *
 *   drbenmclean@gmail.com                 *
 *                                         *
 ******************************************/

#include "TN_Atom.h"

double g_CoulombConst = 8.9875e9;

TN_AtomicParticle::TN_AtomicParticle(){};
TN_AtomicParticle::~TN_AtomicParticle(){};


TN_Nucleus::TN_Nucleus(){};

TN_Nucleus::TN_Nucleus(int nNeutrons, int nProtons, double cx,  double cy, double cz){
		
	m_nProtons = nProtons;
	m_nNeutrons = nNeutrons;
	
	//mass and initial charge
	m_charge = 1.602176634e-19 * m_nProtons; //coulombs
	m_mass = m_nNeutrons * 939.565 + m_nProtons * 938.272; //units:MeV/c^2. Ignores mass loss due to binding energy.
	m_mass = m_mass * 1.6606e-27; //kg
	
	//initial location
	m_location.m_cx = cx;
	m_location.m_cy = cy;
	m_location.m_cz = cz;
};

void TN_Nucleus::addInnerShellElectrons(int n){
	m_nInnerShellElectrons += n;
	m_charge -= n * 1.602176634e-19;
	m_mass += n * 9.1093837015e-31;
};


TN_Electron::TN_Electron(const TN_Nucleus& nucleus, double radius, double spin, double cx,  double cy, double cz){
		
	m_nucleus = &nucleus;
	m_radius = radius;
	m_spin = spin;
	
	//mass and initial charge
	m_charge = -1.602176634e-19; //coulombs
	m_mass = 9.1093837015e-31; //kg
	
	//initial location
	m_location.m_cx = nucleus.m_location.m_cx + cx;
	m_location.m_cy = nucleus.m_location.m_cy + cy;
	m_location.m_cz = nucleus.m_location.m_cz + cz;
	
};


TN_Atom::TN_Atom(int nNeutrons, int nProtons, double cx,  double cy, double cz) : m_nucleus(nNeutrons,nProtons,cx,cy,cz) {};

void TN_Atom::addElectron(double radius, double spin, double cx,  double cy, double cz){
	TN_Electron electron(m_nucleus,radius,spin,cx,cy,cz);
	m_electrons.push_back(electron);
	collapseElectronToRadius(m_electrons.size()-1);
};

const TN_Point& TN_Atom::location() const{
	return m_nucleus.m_location;
};

const double& TN_Atom::nucleusCharge() const{
	return m_nucleus.m_charge;
};

void TN_Atom::addInnerShellElectrons(int n){
	m_nucleus.addInnerShellElectrons(n);
};

double TN_Atom::charge(){
	double c = m_nucleus.m_charge;
	for(unsigned int i=0;i<m_electrons.size();i++){
		c += m_electrons[i].m_charge;
	}
	return c;
};

void TN_Atom::collapseElectronToRadius(int e){
	//https://stackoverflow.com/questions/9604132/how-to-project-a-point-on-to-a-sphere#:~:text=In%20greater%20detail%3A,point%20you%20are%20looking%20for.
	//find projection of point on sphere of radius m_radius
	//x = s + r*(p-s)/(norm(p-s))
	
	double dist = distance(m_electrons[e].m_location, m_nucleus.m_location);
	if(dist != 0.0)
		m_electrons[e].m_location = m_nucleus.m_location
			+ (1.0 / dist) * m_electrons[e].m_radius*(m_electrons[e].m_location-m_nucleus.m_location);

};


void forceBetweenParticles(TN_AtomicParticle &a, TN_AtomicParticle &b){ //force of b on a
	
	if(std::addressof(a) != std::addressof(b)){ //ie not same particle
		
		if(a.m_location == b.m_location){
			
			std::mt19937 gen(clock());
			std::uniform_real_distribution<double> epsilon(-1.0e-13,1.0e-13);
			a.m_location.m_cx = a.m_location.m_cx + epsilon(gen);
			a.m_location.m_cy = a.m_location.m_cy + epsilon(gen);
			a.m_location.m_cz = a.m_location.m_cz + epsilon(gen);
			b.m_location.m_cx = b.m_location.m_cx + epsilon(gen);
			b.m_location.m_cy = b.m_location.m_cy + epsilon(gen);
			b.m_location.m_cz = b.m_location.m_cz + epsilon(gen);
		}
		
		TN_Vector unitvec = normalize(makevector(a.m_location, b.m_location)); //direction of one particle from another
		double separation = distance(a.m_location, b.m_location);

		if(separation==0.0){
			std::cout << "particles have separation zero!" << std::endl;
		}else{
			
			// The electrostatic force bween atoms is defined by Coulomb's law, and governed by their charges
			// and their separation distance. This may be repulsive (like charges) or attractive (differing charges).
			// The minus sign ensures repulsion of like (+/+ or -/-) charges.
			double electrostaticForce = -g_CoulombConst*((a.m_charge * b.m_charge)/pow(separation,2.0));
			
			// Repulsion is caused at very short distances by electron cloud overlap and the Pauli exclusion principle.
			// Repulsion requires definition of the equilibrium bond length (the distance between the nuclei of the bonded atoms)
			// which is the distance at which the attractive and repulsive forces balance out.
			double equBondLength = 4.0e-11;
			// Repulsion is then calculated from an empirically-chosen function that mirrors and balances attraction,
			// with slight differences that allow for attraction to dominate above the equilibrium bond length, and
			// repulsion to dominate at distances less than the equilibrium bond length.
			// Note a better empirically-chosen function for repulsion may allow for more accurate bond lengths to be used.
			double PauliRepulsion = -std::abs(g_CoulombConst*a.m_charge*b.m_charge)*std::pow(equBondLength,2.0)/std::pow(separation,4.0); //always -ve (repulsive)
			
			//find total force on atom
			a.m_force += electrostaticForce*unitvec + PauliRepulsion*unitvec; //allow forces to accumulate

			#ifdef TN_DEBUG
				std::cout << "****************************************" << std::endl;
				std::cout << "Coulomb = " << g_CoulombConst << std::endl;
				std::cout << "mass1 = " << a.m_mass << std::endl;
				std::cout << "mass2 = " << b.m_mass << std::endl;
				std::cout << "charge1 = " << a.m_charge << std::endl;
				std::cout << "charge2 = " << b.m_charge << std::endl;
				std::cout << "force magnitude = " << length(attraction*unitvec + repulsion*unitvec) << std::endl;
				std::cout << "radius = " << separation << std::endl;
				std::cout << "Angstrom radius = " << separation/1e-10 << std::endl;
				std::cout << "attraction = " << attraction << std::endl;
				std::cout << "repulsion = " << repulsion << std::endl;
				std::cout << "****************************************" << std::endl;
			#endif

			// Attempt to downcast to TN_Electron using dynamic_cast, fails if particles are not electrons
			TN_Electron* e_a = dynamic_cast<TN_Electron*>(&a);
			TN_Electron* e_b = dynamic_cast<TN_Electron*>(&b);
			
			if (e_a && e_b) { // Check if both are electrons, if so, add spin interaction to forces
				//electromagnetic force due to electron spin
				double z = 2.0e-10; //the distance at which repulsion equals attraction - shorter range than for coulomb forces
				double numerator = -g_CoulombConst * a.m_charge * b.m_charge * pow(z,2.0) * e_a->m_spin * e_b->m_spin * 4.0;
				//numerator uses same factor as repelling_force to ensure similar scale, even tho not a coulomb effect!
				double electromagnetic_force = numerator / pow(separation,4.0);
				a.m_force = a.m_force + (electromagnetic_force * unitvec);
			}
		}
	}
};
