#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include<iostream>
#include<cmath>
#include<string>
#include<iomanip>

using namespace std;


class Ellipsoid {
public:
	Ellipsoid(); // default parameters for GRS80 ellipsoid
	Ellipsoid(double a, double b);
	Ellipsoid(Ellipsoid & theOther);
	~Ellipsoid() {}

	void setMajorAxis(double a);
	void setMinorAxis(double b);

	double get_a();
	double get_b();

	double calc_f(void);
	double calc_e2(void);
	double calc_eu2(void);

	void printParameters(ostream& );

	double calc_N(const double lat);
	double calc_M(const double lat);

	static double deg2rad(const double deg);
	static double rad2deg(const double rad);

private:
	double a, b;                    // major and semi-major axis
	double f, e2, eu2, N, M;       // flattening, primary eccentcicity, second eccentricity (e')^2,
	static const double PI;       // radius of Curvature in Prime Vertical and Radius of Curvature in Meridian

	void calcEllipsoidParameters();
};

#endif

