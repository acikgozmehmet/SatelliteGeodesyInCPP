#include "Ellipsoid.h"

const double Ellipsoid::PI = 4 * atan(1.0);

Ellipsoid::Ellipsoid() :a(6378137), b(6356752.314140284) { calcEllipsoidParameters(); } // Default Ellipsoid - GRS80
Ellipsoid::Ellipsoid(double a, double b) : a(a), b(b) { calcEllipsoidParameters(); }
Ellipsoid::Ellipsoid(Ellipsoid & theOther) : a(theOther.a), b(theOther.b) { calcEllipsoidParameters(); }


void Ellipsoid::calcEllipsoidParameters(){
	 f = (a - b) / a; // flattening
	 e2 = ( pow(a, 2) - pow(b,2) )  / pow(a,2);
     eu2 = ( pow(a, 2) - pow(b, 2) ) / pow(b,2); // second eccentricity (e')^2
}

void Ellipsoid::setMajorAxis(double a) { this->a = a; }
void Ellipsoid::setMinorAxis(double b) { this->b = b; }

double Ellipsoid::get_a() { return this->a; }
double Ellipsoid::get_b() { return this->b; }

double Ellipsoid::calc_f(void) {  return this->f; }
double Ellipsoid::calc_e2(void) { return this->e2; }
double Ellipsoid::calc_eu2(void) { return this->eu2; }

double Ellipsoid::calc_N(const double lat) {
	double W = sqrt(1 - e2 * sin(lat) * sin(lat));
	N = a / W;
	return N; // radiusOfCurInPrimeVer
}

double  Ellipsoid::calc_M(const double lat) {
	double W = sqrt(1 - e2 * sin(lat) * sin(lat));
	M = N * (1 - e2) / pow(W, 2);
	return M; // radiusOfCurInMeridian
}


 double Ellipsoid::deg2rad(const double deg) {
	 double RHO = 180.0 / PI;
	 return deg / RHO;

}
 double Ellipsoid::rad2deg(const double rad) {
	 double RHO = 180.0 / PI;
	 return rad * RHO;
 }


 void Ellipsoid::printParameters(ostream& w) {
	 w << fixed << setprecision(15);
	 w << "a   : " << this->get_a() << endl;
	 w << "b   : " << this->get_b() << endl;
	 w << "f   : " << this->calc_f() << endl;
	 w << "e2  : " << this->calc_e2() << endl;
	 w << "eu2 : " << this->calc_eu2() << endl;
 }

 /*
 Ellipsoid::Ellipsoid(string ellName) {
 switch (ellName) {
 case "GRS80":
 a = 6378137.0;
 //b = 6356752.31427833;
 b = 6356752.314140284


 calcEllipsoidParameters();
 break;
 case "WGS84":
 a = 6378137.0;
 b = 6356752.314245;
 calcEllipsoidParameters();
 break;
 case "ED50":
 a = 6378388;
 b = 6356911.946;
 calcEllipsoidParameters();
 }
 }
 */
