/*******************************************
 *   Copyright (C) 2022 by Ben F McLean    *
 *   drbenmclean@gmail.com                 *
 *                                         *
 ******************************************/

#include <iostream>
#include <random>
#include <cassert>
#include <filesystem>

#include "TN_Experiment.h"

using namespace std;

/*
main() currently generates an experiment modelling the burning of carbon in excess oxygen.
This is a covalent bonding process, but rather than explicitly modelling the sharing of
electrons, the modelling uses a heuristic - initially modelling carbon and oxygen atoms as
varying in electronegativity. Electronegativity was simulated by giving the reaction
constitutents a partial charge approximating the charge they would assume once covalently
bonded. This charge allows for the attraction and bonding of atoms.
*/


int main(){

	cout << "define output directory..." << endl;
	string output_path = "./output/";
	
	cout << "create experiment..." << endl;
	TN_Experiment experiment01;
	
	cout << "establish experiment boundaries..." << endl;
	experiment01.setBoundaries(-300.0e-12,300.0e-12,-300.0e-12,300.0e-12,-300.0e-12,300.0e-12);

	//atom rand coords
	std::mt19937 gen(std::random_device{}());
	uniform_real_distribution<double> rand_coord(-300.0e-12,300.0e-12);
	
	cout << "create one carbon and two hydrogen atoms each loop..." << endl;
	for(int i=0;i<12;i++){
		
		//create carbon atom
		TN_Atom carbon_01(6,6,rand_coord(gen),rand_coord(gen),rand_coord(gen));
		carbon_01.addInnerShellElectrons(2); //add less electrons to simulate lower electronegativity

		//create oxygen atoms
		TN_Atom oxygen_01(8,8,rand_coord(gen),rand_coord(gen),rand_coord(gen));
		TN_Atom oxygen_02(8,8,rand_coord(gen),rand_coord(gen),rand_coord(gen));
		oxygen_01.addInnerShellElectrons(10); //add more electrons to simulate higher electronegativity
		oxygen_02.addInnerShellElectrons(10);

		//add atoms to the experiment
		experiment01.addAtom(carbon_01);
		experiment01.addAtom(oxygen_01);
		experiment01.addAtom(oxygen_02);
	}
	
	cout << "give thermal velocities..." << endl;
	experiment01.giveThermalVelocities(2000.0); //A red Bunsen burner flame burns at 2000 deg Kelvin
	
	cout << "timestep..." << endl;
	int nt = 1000;
	int snapshotInterval = 2;
	double dt = 5.0e-15;
	/*Note: timesteps must allow for incremental changes in position leading to incremental changes in
	interatomic forces. If dt is slightly too large, attractive atoms may interact by oscillating back
	and forth, as they repeatedly attract, move too close, and repel, before attracting again. If dt is
	much too large, atoms may instead interact by 'slingshotting' off one another, as attraction leads
	them into extremely close proximity in one step, resulting in extreme repulsion and dispersion. If
	dt is too small, neither oscillation nor slingshotting occurs, rather the atoms will move very slowly.
	*/

	int snap = 0;
	for(int i=0;i<nt;i++){
		experiment01.timestep(dt);
		if(i%snapshotInterval==0){
			string file = output_path + to_string(snap) + ".csv";
			experiment01.saveExperimentData(file);
			snap++;
		}
	}
	
  	return EXIT_SUCCESS;
};

