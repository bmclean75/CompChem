/*******************************************
 *   Copyright (C) 2022 by Ben F McLean    *
 *   drbenmclean@gmail.com                 *
 *                                         *
 ******************************************/

#include "TN_Experiment.h"

double g_Viscosity = 0.6;
double g_ElectronViscosity = 0.6;

TN_Experiment::TN_Experiment(){};

void TN_Experiment::setBoundaries(double cx_min, double cx_max, double cy_min, double cy_max, double cz_min, double cz_max){
	
	m_cx_min = cx_min;
	m_cx_max = cx_max;
	m_cy_min = cy_min;
	m_cy_max = cy_max;
	m_cz_min = cz_min;
	m_cz_max = cz_max;
	
};


void TN_Experiment::constrainToBoundaries(){
	for(unsigned int i=0;i<m_atoms.size();i++){
		if(m_atoms[i].m_nucleus.m_location.m_cx < m_cx_min){
			m_atoms[i].m_nucleus.m_location.m_cx = m_cx_min;
		}
		if(m_atoms[i].m_nucleus.m_location.m_cx > m_cx_max){
			m_atoms[i].m_nucleus.m_location.m_cx = m_cx_max;
		}
		
		if(m_atoms[i].m_nucleus.m_location.m_cy < m_cy_min){
			m_atoms[i].m_nucleus.m_location.m_cy = m_cy_min;
		}
		if(m_atoms[i].m_nucleus.m_location.m_cy > m_cy_max){
			m_atoms[i].m_nucleus.m_location.m_cy = m_cy_max;
		}
		
		if(m_atoms[i].m_nucleus.m_location.m_cz < m_cz_min){
			m_atoms[i].m_nucleus.m_location.m_cz = m_cz_min;
		}
		if(m_atoms[i].m_nucleus.m_location.m_cz > m_cz_max){
			m_atoms[i].m_nucleus.m_location.m_cz = m_cz_max;
		}
	}
};

void TN_Experiment::bounceFromBoundaries(){
	
	std::mt19937 gen(std::random_device{}());
	std::uniform_real_distribution<double> epsilon(0,std::min(std::min(m_cx_max-m_cx_min, m_cy_max-m_cy_min), m_cz_max-m_cz_min)/100.0);
	//(make epsilon a uniform real distribution over 1/100th of the minimum dimensions of the experiment)
	//(needed so that returning an errant atom to inside the experiment does not result in multiple atoms
	//being returned to the one corner point, and being co-located)
	
	for(unsigned int i=0;i<m_atoms.size();i++){
		
		bool bouncing = false;
		int multibounce = 0;
		TN_Vector nvec; //the normal vector to the wall off which atom is bouncing
		if(m_atoms[i].m_nucleus.m_location.m_cx < m_cx_min){
			
			//bounce off plane x = m_cx_min, likewise below...
			m_atoms[i].m_nucleus.m_location.m_cx = m_cx_min + epsilon(gen);
			nvec.set(1.0,0.0,0.0);
			bouncing = true;
			multibounce++;
		}
		if(m_atoms[i].m_nucleus.m_location.m_cx > m_cx_max){
			m_atoms[i].m_nucleus.m_location.m_cx = m_cx_max - epsilon(gen);
			nvec.set(-1.0,0.0,0.0);
			bouncing = true;
			multibounce++;
		}
		
		if(m_atoms[i].m_nucleus.m_location.m_cy < m_cy_min){
			m_atoms[i].m_nucleus.m_location.m_cy = m_cy_min + epsilon(gen);
			nvec.set(0.0,1.0,0.0);
			bouncing = true;
			multibounce++;
		}
		if(m_atoms[i].m_nucleus.m_location.m_cy > m_cy_max){
			m_atoms[i].m_nucleus.m_location.m_cy = m_cy_max - epsilon(gen);
			nvec.set(0.0,-1.0,0.0);
			bouncing = true;
			multibounce++;
		}
		
		if(m_atoms[i].m_nucleus.m_location.m_cz < m_cz_min){
			m_atoms[i].m_nucleus.m_location.m_cz = m_cz_min + epsilon(gen);
			nvec.set(0.0,0.0,1.0);
			bouncing = true;
			multibounce++;
		}
		if(m_atoms[i].m_nucleus.m_location.m_cz > m_cz_max){
			m_atoms[i].m_nucleus.m_location.m_cz = m_cz_max - epsilon(gen);
			nvec.set(0.0,0.0,-1.0);
			bouncing = true;
			multibounce++;
		}
		
		if(multibounce>1){ 	//atom is out of bounds in 2 or 3 of the x, y, and z directions
							//so send vector straight back from where it came...
			nvec.set(-m_atoms[i].m_nucleus.m_velocity.m_cx,-m_atoms[i].m_nucleus.m_velocity.m_cy,-m_atoms[i].m_nucleus.m_velocity.m_cz);
			nvec = normalize(nvec);
		}
		
		if(bouncing == true){
			
			// Reflect the atom by giving it a reoriented velocity, reflected in plane defined by the normal vector.
			// r = d -2.0*dotprod(d.n)*n;
			// https://math.stackexchange.m_com/questions/13261/how-to-get-a-reflection-vector
			m_atoms[i].m_nucleus.m_velocity = m_atoms[i].m_nucleus.m_velocity
				- 2.0*dotproduct(m_atoms[i].m_nucleus.m_velocity,nvec)*nvec;
			
			for(unsigned int j=0;j<m_atoms[i].m_electrons.size();j++){
				// electrons must reflect off same wall also:
				m_atoms[i].m_electrons[j].m_velocity = m_atoms[i].m_electrons[j].m_velocity
					- 2.0*dotproduct(m_atoms[i].m_electrons[j].m_velocity,nvec)*nvec;
				m_atoms[i].collapseElectronToRadius(j);
			}
		}
	}
};

void TN_Experiment::addAtom(const TN_Atom &atom){
	m_atoms.push_back(atom);
};

void TN_Experiment::giveThermalVelocities(double tempKelvin){
	
	// Produce random numbers via the Mersenne twister algorithm
	std::mt19937 gen(std::random_device{}());
	
	// Create a uniform distribution for generating unit vectors of direction
	std::uniform_real_distribution<double> dirDistribution(-1.0,1.0);
	
	for(unsigned int i=0;i<m_atoms.size();i++){
		
		// Calculate the mean magnitude velocity of such a particle, dependent on mass and temperature
		// https://en.wikipedia.org/wiki/Thermal_velocity
		double meanVelocity = std::sqrt((8.0*1.380649e-23*tempKelvin)/(M_PI * m_atoms[i].m_nucleus.m_mass));

		// Create a normal distribution of velocities around meanVelocity, and find
		// a random velociuty for this particle
		std::normal_distribution<double> velocityDistribution(0.0, 2.0 * meanVelocity);
		double thisVelocity = velocityDistribution(gen);

		// Generate a random-direction unit vector
		TN_Vector randvec(dirDistribution(gen),dirDistribution(gen),dirDistribution(gen));
		randvec = normalize(randvec);

		// Give particle the random velocity in the random direction
		m_atoms[i].m_nucleus.m_velocity = thisVelocity*randvec;
	}
	
};


void TN_Experiment::timestep(double dt){
	
	//zero all forces so they dont accumulate timestep-to-timestep
	//(but forces still allowed to accumulate particle-to-particle, within each timestep)
	for(unsigned int a1=0;a1<m_atoms.size();a1++){
		m_atoms[a1].m_nucleus.m_force.set(0.0);
		
		for(unsigned int e1=0;e1<m_atoms[a1].m_electrons.size();e1++){
			m_atoms[a1].m_electrons[e1].m_force.set(0.0);
		}
	}

	//force on nuclei from nuclei
	for(unsigned int a1=0;a1<m_atoms.size();a1++){
		for(unsigned int a2=0;a2<m_atoms.size();a2++){
			forceBetweenParticles(m_atoms[a1].m_nucleus,m_atoms[a2].m_nucleus);
		}
	}
	
	//force of nuclii on electrons and vice versa
	for(unsigned int a1=0;a1<m_atoms.size();a1++){
		for(unsigned int e1=0;e1<m_atoms[a1].m_electrons.size();e1++){
			for(unsigned int a2=0;a2<m_atoms.size();a2++){
				if(a1!=a2){
					//force on nuclei from electrons
					forceBetweenParticles(m_atoms[a2].m_nucleus,m_atoms[a1].m_electrons[e1]);
					//force on electrons from nuclei
					forceBetweenParticles(m_atoms[a1].m_electrons[e1],m_atoms[a2].m_nucleus);
				}
			}
		}
	}
	
	//force on electrons from electrons
	for(unsigned int a1=0;a1<m_atoms.size();a1++){
		for(unsigned int e1=0;e1<m_atoms[a1].m_electrons.size();e1++){
			for(unsigned int a2=0;a2<m_atoms.size();a2++){
				for(unsigned int e2=0;e2<m_atoms[a2].m_electrons.size();e2++){
					forceBetweenParticles(m_atoms[a1].m_electrons[e1],m_atoms[a2].m_electrons[e2]);
				}
			}
		}
	}
	
	//force on electrons transmitted to own nuclei, and vice versa
	for(unsigned int a1=0;a1<m_atoms.size();a1++){
		for(unsigned int e1=0;e1<m_atoms[a1].m_electrons.size();e1++){
			TN_Vector force_on_electrons = m_atoms[a1].m_nucleus.m_force;
			TN_Vector force_on_nucleus = m_atoms[a1].m_electrons[e1].m_force;
			
			m_atoms[a1].m_electrons[e1].m_force += force_on_electrons;
			m_atoms[a1].m_nucleus.m_force += force_on_nucleus;
		}
	}
	
	for(unsigned int i=0;i<m_atoms.size();i++){
		
		// Calculate new acceleration, position, and velocity, using equations of motion
		
		//acceleration of nuclii, a = F/m
		m_atoms[i].m_nucleus.m_acceleration = m_atoms[i].m_nucleus.m_force / m_atoms[i].m_nucleus.m_mass;
		
		//new position of nuclii, s1 = s0 + V0.t + 1/2.a.t^2
		m_atoms[i].m_nucleus.m_location = m_atoms[i].m_nucleus.m_location + m_atoms[i].m_nucleus.m_velocity*dt
			+ (m_atoms[i].m_nucleus.m_acceleration*std::pow(dt,2.0))/2.0;
			
		//new velocity of nuclii, v1 = v0 + a.t
		m_atoms[i].m_nucleus.m_velocity = m_atoms[i].m_nucleus.m_velocity + m_atoms[i].m_nucleus.m_acceleration*dt;
		
		//viscosity term to dampen velocities
		m_atoms[i].m_nucleus.m_velocity = m_atoms[i].m_nucleus.m_velocity*g_Viscosity;
	}

	for(unsigned int i=0;i<m_atoms.size();i++){
		for(unsigned int j=0;j<m_atoms[i].m_electrons.size();j++){
			
			//acceleration of electrons, a = F/m
			m_atoms[i].m_electrons[j].m_acceleration = m_atoms[i].m_electrons[j].m_force / m_atoms[i].m_electrons[j].m_mass;
			
			//new position of electrons, s1 = s0 + V0.t + 1/2.a.t^2
			m_atoms[i].m_electrons[j].m_location = m_atoms[i].m_electrons[j].m_location + m_atoms[i].m_electrons[j].m_velocity*dt
				+ (m_atoms[i].m_electrons[j].m_acceleration*pow(dt,2.0))/2.0;
				
			//collapse electron position to radius of shell
			//_atoms[i].m_electrons[j].m_collapse_electron_to_radius();
			m_atoms[i].collapseElectronToRadius(j);
			
			//new velocity of electrons, v1 = v0 + a.t
			m_atoms[i].m_electrons[j].m_velocity = m_atoms[i].m_electrons[j].m_velocity + m_atoms[i].m_electrons[j].m_acceleration*dt;
			
			//viscosity term to dampen velocities
			m_atoms[i].m_electrons[j].m_velocity = m_atoms[i].m_electrons[j].m_velocity*g_ElectronViscosity;
		}
	}

	//ensure experimental particles are retained within experimental container by
	//returning to container wall and reflecting velocity back into container
	bounceFromBoundaries();

};

//write experiment to file
void TN_Experiment::saveExperimentData(const std::string& filename){
	std::ofstream out;
	out.open(filename,std::ios::out);
	if(out.bad()){
		std::cerr << "ERROR: cannot open outfile " << std::endl;
		exit(8);
	}
	out << "x,y,z,e,mass,relsize" << std::endl;
	
	//relative but not-to-scale sizing of atoms
	std::vector<double> rel_size;
	double max_size = 10.0;
	double min_size = 8.0;
	double max = m_atoms[0].m_nucleus.m_mass;
	double min = m_atoms[0].m_nucleus.m_mass;
	
	for(unsigned int i=1;i<m_atoms.size();i++){
		if(m_atoms[i].m_nucleus.m_mass > max){max = m_atoms[i].m_nucleus.m_mass;}
		if(m_atoms[i].m_nucleus.m_mass < min){min = m_atoms[i].m_nucleus.m_mass;}
	}
	for(unsigned int i=0;i<m_atoms.size();i++){
		if(max==min){ //all nuclii the same size
			rel_size.push_back(10.0);
		}else{
			rel_size.push_back(min_size + (max_size-min_size) * (m_atoms[i].m_nucleus.m_mass-min) / (max-min));
		}
	}
	
	// Print out atoms
	for(unsigned int i=0;i<m_atoms.size();i++){
		
		out << m_atoms[i].m_nucleus.m_location.m_cx << "," 
			<< m_atoms[i].m_nucleus.m_location.m_cy << "," 
			<< m_atoms[i].m_nucleus.m_location.m_cz 
			<< "," << m_atoms[i].m_nucleus.m_charge*1.0e19 
			<< "," << m_atoms[i].m_nucleus.m_mass*1.0e23 
			<<  "," << rel_size[i] << std::endl;
		
		
		//for drawing electrons
		for(unsigned int j=0;j<m_atoms[i].m_electrons.size();j++){
				
			double atomRadius = rel_size[i];	//this is the radius given to the atom for display purposes...
												//electrons should be displayed protruding from this surface

			//find a scaled vector from the center of the atom to the displayed surface of the atom
			TN_Vector vectorToElectron = makevector(m_atoms[i].m_nucleus.m_location, m_atoms[i].m_electrons[j].m_location);
			TN_Vector scaledVectorToElectron = atomRadius * normalize(vectorToElectron);
			TN_Point electronLocation = m_atoms[i].m_nucleus.m_location + scaledVectorToElectron;

			out << electronLocation.m_cx << "," 
				<< electronLocation.m_cy << "," 
				<< electronLocation.m_cz << "," 
				<< -1.602176634e-19 << "," 
				<< 9.1093837015e-31 << ","
				<< 1.0 << std::endl; //electrons have rel size = 1.0

		}
	}
	
	out << std::endl;
	out << std::flush;
	out.close();

};

