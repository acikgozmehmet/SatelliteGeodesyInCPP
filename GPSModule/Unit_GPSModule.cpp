#include<iostream>
#include "..\Constants\Constants.h"
#include "GPSModule.h"

const double PI = atan(1.0) * 4;
using namespace std;

int main() {

	{
		cout << "Test-1" << endl;

		Car P1, P3;
		Geo P2;

		Ellipsoid ell(6378137.0000, 6356752.3142);
		// ANKR 20805M002    4121948.52404  2652187.90210  4069023.75864
		P1.x = 4121948.52404;
		P1.y = 2652187.90210;
		P1.z = 4069023.75864;
		cout << fixed << endl;
		P1.print(cout); cout << endl;

		P2 = GPSModule::Car2Geo(ell, P1);
		P2.print(cout); cout << endl;

		P3 = GPSModule::Geo2Car(ell, P2);
		P3.print(cout); cout << endl;
	}

	{
		cout << "Test-2" << endl;


		Car sta, sat;
		sta.x = -1283387.6627;
		sta.y = -4713015.5351;
		sta.z = 4090190.4704;
		cout << "Position 1 :"; 
		sta.print(cout);
		
		Geo P2;
		Ellipsoid ell;


		P2 = GPSModule::Car2Geo(ell, sta);
		P2.print(cout); cout << endl;
		cin.get();


		sat.x = -11918707.0438537;
		sat.y = -23755879.2115591;
		sat.z = -75141.8219498;
		
		cout << "Position 2 :";
		sat.print(cout); cout << endl;


		Matrix NEU = GPSModule::xyz2neu(sta, sat);
		NEU.printToScreen();
		double alpha = GPSModule::azimuthOfSat(NEU);
		double elev = GPSModule::elevAngleSat(NEU);
		cout << setprecision(8) << alpha * (180/PI)<< endl;
		cout << elev * (180 / PI) << endl;


		cout << "elevation Angle of the Satellite : " << GPSModule::elevAngleSat(sta, sat)*(180 / PI);
		cout << endl;
		cout << "azimuth of the Satellite : " << GPSModule::azimuthOfSat(sta, sat)*(180/PI);

		cout << "Saastamonien : " << GPSModule::Saastamoinen(P2.hgt, elev) << endl;


	}


	{ cout << endl << "Test - 3" << endl;


	cout << " Ineterpolation of B values from the table given an input of height ranging from 0 to 5.4 km" << endl;
	/*
	for (double x = 0.0; x < 5.4; x += 0.1) {
		cout << x << " " << GPSModule::interpolateBValue(x) << endl;
	}
	*/
	
	
	/*
	for (double T = 0; T < 100; T++) {
		double e = Rh * exp(-37.2465 + 0.213166 * T - 0.000256908 * pow(T, 2));

	}
	*/

	
	cin.ignore();
	cin.get();


	}

	
	
	
	return 0;
}