/*******************************************
 *   Copyright (C) 2022 by Ben F McLean    *
 *   drbenmclean@gmail.com                 *
 *                                         *
 ******************************************/
 
#ifndef TN_MATH
#define TN_MATH

#include <iostream>
#include <cmath>


/*
TN_MATH provides a very simple TN_Point and TN_Vector class, along with non-class inlined functions
for all necessary mathematical functions between them.

Points and vectors are both defined by a triplet of signed floating-point values. Hence a base class 
TN_Triplet is defined, from which both TN_Point and TN_Vector are derived.

Note, the TN_Point and TN_Vector class are identical, as they both derive from the TN_Triplet class,
however defining them as different classes allows strict control over what functions can be called for
what quantities, e.g., you can't find the dot product of two points, only two vectors.
*/

// Base struct representing a triplet of values.
struct TN_Triplet {
	
	double m_cx, m_cy, m_cz;

	// Constructor to initialize the triplet.
	TN_Triplet() : m_cx(0.0), m_cy(0.0), m_cz(0.0) {}
	TN_Triplet(double cx, double cy, double cz) : m_cx(cx), m_cy(cy), m_cz(cz) {}

	// Function to display the triplet values.
	inline void display() const {
		std::cout << "(" << m_cx << ", " << m_cy << ", " << m_cz << ")" << std::endl;
	}

	// Set an already-existing triplet's coordinates
	inline void set(double cx, double cy, double cz){
		m_cx = cx;
		m_cy = cy;
		m_cz = cz;
	}

	// Set an already-existing triplet's coordinates
	inline void set(double d){
		m_cx = d;
		m_cy = d;
		m_cz = d;
	}

	template<class tn_class>
	inline TN_Triplet& operator+=(const tn_class& t2) {
		m_cx += t2.m_cx;
		m_cy += t2.m_cy;
		m_cz += t2.m_cz;
		return *this;
	}

	// Function to check equality of two triplets
	inline bool operator==(const TN_Triplet &t2) const {
		if(this->m_cx == t2.m_cx && this->m_cy == t2.m_cy && this->m_cz == t2.m_cz)
			return true;
		return false;
	}

};

// Derived class representing a point, inherits from TN_Triplet.
class TN_Point : public TN_Triplet {
	
	public:

	// Constructor to initialize a point, uses base class constructor.
	TN_Point() : TN_Triplet() {}
	TN_Point(double cx, double cy, double cz) : TN_Triplet(cx, cy, cz) {}

	// Additional methods specific to points can be added here.
	inline void displayPoint() const {
		std::cout << "Point: ";
		display(); // Call base class display method.
	}
};

// Derived class representing a vector, inherits from TN_Triplet.
class TN_Vector : public TN_Triplet {
	
	public:
	
	// Constructor to initialize a vector, uses base class constructor.
	TN_Vector() : TN_Triplet() {}
	TN_Vector(double cx, double cy, double cz) : TN_Triplet(cx, cy, cz) {}

	// Additional methods specific to vectors can be added here.
	inline void displayVector() const {
		std::cout << "Vector: ";
		display(); // Call base class display method.
	}

};


/***************************************
*		   Math Operations for         *
*		TN_Points and TN_Vectors       *
***************************************/

// Basic arithmetic operators, templated for either point or vector
//*****************************************************************

template<class tn_class>
inline tn_class operator+(const tn_class& t1, const tn_class& t2) {
	return tn_class(t1.m_cx + t2.m_cx, t1.m_cy + t2.m_cy, t1.m_cz + t2.m_cz);
}

template<class tn_class>
inline tn_class operator-(const tn_class& t1, const tn_class& t2) {
	return tn_class(t1.m_cx - t2.m_cx, t1.m_cy - t2.m_cy, t1.m_cz - t2.m_cz);
}

template<class tn_class>
inline tn_class operator*(double scalar, const tn_class& t) {
	return tn_class(scalar * t.m_cx, scalar * t.m_cy, scalar * t.m_cz);
}

template<class tn_class>
inline tn_class operator*(const tn_class& t, double scalar) {
	return tn_class(scalar * t.m_cx, scalar * t.m_cy, scalar * t.m_cz);
}

template<class tn_class>
inline tn_class operator/(const tn_class& t, double scalar) {
	return tn_class(t.m_cx / scalar, t.m_cy / scalar, t.m_cz / scalar);
}

//calculate the distance between two points
//*****************************************

inline double distance(const TN_Point &p1, const TN_Point &p2){
	return sqrt(pow(p2.m_cx - p1.m_cx, 2.0) + pow(p2.m_cy - p1.m_cy, 2.0) + pow(p2.m_cz - p1.m_cz, 2.0));
};

//generate a vector from two points
//*********************************

inline TN_Vector makevector(const TN_Point &p1, const TN_Point &p2){
	
	TN_Vector vec;
	vec.m_cx = p2.m_cx - p1.m_cx;
	vec.m_cy = p2.m_cy - p1.m_cy;
	vec.m_cz = p2.m_cz - p1.m_cz;
	
	return vec;
};

// Point addition of a vector
//***************************

inline TN_Point operator+( const TN_Point& p, const TN_Vector& v) {
	return TN_Point(p.m_cx + v.m_cx, p.m_cy + v.m_cy, p.m_cz + v.m_cz);
}

//calculate the length of a vector
//********************************

inline double length(const TN_Vector &v){
	return sqrt(pow(v.m_cx,2.0) + pow(v.m_cy,2.0) + pow(v.m_cz,2.0));
};

//calculate the normalized version of a vector
//********************************************

inline TN_Vector normalize(const TN_Vector &v){
	
	TN_Vector normvec;
	double length = sqrt(pow(v.m_cx,2.0) + pow(v.m_cy,2.0) + pow(v.m_cz,2.0));
	if(length != 0.0)
		normvec = (1.0/length) * v;
	else
		normvec = v;
	return normvec;
};

//calculate cross product of two vectors
//**************************************
//(note: the cross product of two vectors is the normal of the plane defined by the two vectors)

inline TN_Vector crossproduct(const TN_Vector &v1, const TN_Vector &v2){
	
	TN_Vector cpvec;
	cpvec.m_cx = v1.m_cy * v2.m_cz - v1.m_cz * v2.m_cy;
	cpvec.m_cy = v1.m_cz * v2.m_cx - v1.m_cx * v2.m_cz;
	cpvec.m_cz = v1.m_cx * v2.m_cy - v1.m_cy * v2.m_cx;
	
	return cpvec;
};

//calculate dot product of two vectors
//************************************
//a<dot>b = a*b*cos(theta)

inline double dotproduct(const TN_Vector &v1, const TN_Vector &v2){
	
	double dpvec;
	dpvec = v1.m_cx * v2.m_cx + v1.m_cy * v2.m_cy + v1.m_cz * v2.m_cz;
	return dpvec;
};


#endif
