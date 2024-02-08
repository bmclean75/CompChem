/*******************************************
 *   Copyright (C) 2022 by Ben F McLean    *
 *   drbenmclean@gmail.com                 *
 *                                         *
 ******************************************/

#include <vector>
#include <iostream>
#include <fstream>

#include "TN_Atom.h"

#ifndef TN_EXPERIMENT
#define TN_EXPERIMENT

/*
Class TN_Experiment defines an experiment by the physical space it contains, and the atoms
present within the experiment. It then provides functions for containing all atoms within
the experiment, adding atoms, providing atoms with a thermally-defined initial velocity,
and for stepping the experiment forward in time.
*/

extern double g_Viscosity;
extern double g_ElectronViscosity;

class TN_Experiment {

	public:
	
	double m_cx_min, m_cx_max, m_cy_min, m_cy_max, m_cz_min, m_cz_max;
	std::vector<TN_Atom> m_atoms;

	// Default constructor
	TN_Experiment();

	// Sets the spatial boundaries of the experiment
	void setBoundaries(double cx_min, double cx_max, double cy_min, double cy_max, double cz_min, double cz_max);

	// Ensures all atoms are within the defined spatial boundaries
	void constrainToBoundaries();

	// Bounces any atoms moving across a boundary back into the experiment
	void bounceFromBoundaries();

	// Adds a TN_Atom into the experiment
	void addAtom(const TN_Atom &atom);

	// Provides all atoms in the experiment with a random velocity
	// Default represents a typical ambient temperature of 20 deg Celcius
	void giveThermalVelocities(double tempKelvin = 293.15);

	// Moves the atoms in the experiment forward by one timestep
	void timestep(double dt);
	
	// Saves the experiment data to a specified text file, in csv format
	void saveExperimentData(const std::string& filename);

};

#endif