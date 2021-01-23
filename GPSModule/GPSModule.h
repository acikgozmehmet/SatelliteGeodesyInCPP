#ifndef  GPSMODULE_H
#define GPSMODULE_H

#include<iostream>
#include<cmath>
#include<vector>
#include<iomanip>
#include "..\Constants\Constants.h"
#include "..\Ellipsoid2\Ellipsoid.h"
#include "..\Matrix\Matrix.h"
#include "..\BrdcEph\BrdcEph.h"
#include "..\RinexObsFile\RinexObsFile.h"


class GPSModule {
public:
	static vector<double> formP3(vector<double>& P1, vector<double>& P2);
	static double ArcTan(double x, double y);

	static Geo Car2Geo(Ellipsoid& elpsd, const Car& point);
	static Geo Car2Geo(const Car& point);

	static Car Geo2Car(Ellipsoid& elpsd, const Geo& point);
	static Car Geo2Car(const Geo& point);

	static Matrix xyz2neu(Car& point1, Car& p2);

	static double elevAngleSat(Car& point1, Car& point2);
	static double elevAngleSat(Matrix& NEU);

	static double azimuthOfSat(Car& point1, Car& point2);
	static double azimuthOfSat(Matrix& NEU);

	static double Saastamoinen(double height, double elev);

	static double CelciusToKelvin(double celcius);
	static double KelvinToCelcius(double celcius);
	
	static double interpolateBValue(double x);

	//static Matrix earthRotationCorrection(double t);

	static bool sppEpochbyEpoch(BrdcEph & brdc, RinexObsFile & rinexo, ostream& w = cout);
	
	static vector<double> lightTimeIteration(const vector<int>& svList, const int& GPSWeek,
		                                    const double& tReceiverTrue, BrdcEph& brdc, Car& appPos);

private:
	static const double PI;

};


#endif // ! GPSMODULE_H