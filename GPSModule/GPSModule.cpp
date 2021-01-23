#include "GPSModule.h"
#include "..\AdjustmentModule2\AdjustmentModule.h"


const double GPSModule::PI = atan(1.0) * 4;

vector<double> GPSModule::formP3(vector<double>& P1, vector<double>& P2) {
	vector<double> P3;
	
	double FREQ1 = constants::FREQ1;
	double FREQ2 = constants::FREQ2;

	for (unsigned int i = 0; i < P1.size(); i++) {
		double p = (pow(FREQ1, 2)* P1[i] - pow(FREQ2, 2)*P2[i]) / (pow(FREQ1, 2) - pow(FREQ2, 2));
		P3.push_back(p);
	}

	return P3;
}



/**
Converts cartesian coordinates (x,y,z) into
ellipsoidal coordinates (lat,lon,alt) on the given ellipsoid
according to iterative method.
param@ Ellipsoid
param@ Cartesian coordinates (in meters)
return Ellipsoidal coordinates (lat, lon in radians, hgt in meters)
*/
Geo GPSModule::Car2Geo(Ellipsoid& elpsd, const Car& point) {

	Geo T;

	const double RHO = 180.0 / PI;
	const double EPS = 1e-10;

	double a = elpsd.get_a();
	double b = elpsd.get_b();
	double e2 = elpsd.calc_e2();
	double eu2 = elpsd.calc_eu2();

	double p = sqrt(pow(point.x, 2) + pow(point.y, 2));

	double tanu = (point.z / p) * (a / b);

	double cosu, sinu, tanuPrev, lat;

	cout << fixed << setprecision(15);

	do {
		tanuPrev = tanu;
		cosu = sqrt(1 / (1 + pow(tanu, 2)));
		sinu = sqrt(1 - pow(cosu, 2));
		lat = ArcTan((point.z + eu2 * b*pow(sinu, 3)), (p - e2 * a * pow(cosu, 3)));
		tanu = (b / a) * tan(lat);
	} while (abs(tanu - tanuPrev) > EPS);


	double N = elpsd.calc_N(lat);
	T.lat = lat;
	T.lon = ArcTan(point.y, point.x);

	if (abs(lat - PI / 2) * RHO  > EPS) { // as long as lat != 90 
		T.hgt = p / cos(lat) - N;
	}
	else if (abs(lat) * RHO > EPS) { // as long as lat != 0
		T.hgt = point.z / sin(lat) - N + e2 * N;
	}

	return T;
}

/**
Converts cartesian coordinates (x,y,z) into
ellipsoidal coordinates (lat,lon,alt) on the default ellipsoid- GRS80 
according to iterative method.
param@ Cartesian coordinates (in meters)
return Ellipsoidal coordinates (lat, lon in radians, hgt in meters)
*/
Geo GPSModule::Car2Geo(const Car& point) {
	Geo T;
	Ellipsoid elpsd; // Default Elipsoid - GRS80
	T = Car2Geo(elpsd, point);
	return T;
}


/**
Converts ellipsoidal coordinates (lat,lon,alt) into
cartesian coordinates (x,y,z) on the given ellipsoid
param@ Ellipsoid
param@ Ellipsoidal coordinates (lat, lon in radians, hgt in meters)
return Cartesian coordinates (in meters)
*/
Car GPSModule::Geo2Car(Ellipsoid& elpsd, const Geo& point) {
	Car T;

	double a = elpsd.get_a();
	double b = elpsd.get_b();
	double N = elpsd.calc_N(point.lat);

	T.x = ( N + point.hgt) * cos(point.lat) * cos(point.lon);
	T.y = ( N + point.hgt) * cos(point.lat) * sin(point.lon);
	T.z = (N * pow(b, 2) / pow(a, 2) + point.hgt) * sin(point.lat);

	return T;
}


/**
Converts ellipsoidal coordinates (lat,lon,alt) into
cartesian coordinates (x,y,z) on the default ellipsoid - GRS80
param@ Ellipsoidal coordinates (lat, lon in radians, hgt in meters)
return Cartesian coordinates (in meters)
*/
Car GPSModule::Geo2Car(const Geo& point) {
	Car T;
	Ellipsoid elpsd;
	T = Geo2Car(elpsd, point);
	return T;
}




/*
returns the Arctan function without considering the region problem....
*/
double GPSModule::ArcTan(double x, double y) {
	const double RHO = 180.0 / PI;
	double teta;
	teta = atan2(x, y);
	(teta < 0) ? (teta = teta + 2 * PI) : (teta = teta);
	return teta;
}


/**
Converts the Cartesian Coordinate System (x, y, z) to Local Geodetic System (n,e,u)
Leick, Satellite Surveying, Page 64
*/

Matrix GPSModule::xyz2neu(Car& point1, Car& point2) {

	// NEU = G * DX

	Matrix NEU;
	Matrix DX(3, 1);
	DX.setValueAtElement(0, 0, point2.x-point1.x);
	DX.setValueAtElement(1, 0, point2.y-point1.y);
	DX.setValueAtElement(2, 0, point2.z-point1.z);


	Geo p1 = Car2Geo(point1);
	Matrix G(3, 3);
	double sf = sin(p1.lat);
	double cf = cos(p1.lat);
	double sl = sin(p1.lon);
	double cl = cos(p1.lon);

	G.setValueAtElement(0, 0, -sf * cl);
	G.setValueAtElement(0, 1, -sf * sl);
	G.setValueAtElement(0, 2, cf );

	G.setValueAtElement(1, 0, -sl);
	G.setValueAtElement(1, 1,  cl);

	G.setValueAtElement(2, 0, cf * cl);
	G.setValueAtElement(2, 1, cf * sl);
	G.setValueAtElement(2, 2, sf);

	NEU = G * DX;

	return NEU;
}

/**
Returns the elevation angle of a satellite
param@ point1, the cartesian coordinates of the station - ECEF (in meters)
param@ point2, the cartesian coordinates of the station - ECEF (in meters)
return elevation angle of the satellite (in radians)
Leick, Alfred. GPS Satellite Surveying, page 46, (Eq. 2.99)
*/
double GPSModule::elevAngleSat(Car& point1, Car& point2) {
	double dx = point2.x - point1.x;
	double dy = point2.y - point1.y;
	double dz = point2.z - point1.z;
	
	double s = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));

	Geo p1 = Car2Geo(point1);
	double term1 = cos(p1.lat) * cos(p1.lon) * dx + cos(p1.lat) * sin(p1.lon) * dy + sin(p1.lat) * dz;
	double elevAngle = asin(term1/s);
	
	return elevAngle;	
}

/**
Returns the elevation angle of a satellite
param@ Coordinates in LGS (NEU) at the station
return elevation angle of the satellite (in radians)
Leick, Alfred. GPS Satellite Surveying, page 45, (Eq. 2.93)
*/
double GPSModule::elevAngleSat(Matrix& NEU) {
	double n = NEU.getValueAtElement(0, 0);
	double e = NEU.getValueAtElement(1, 0);
	double u = NEU.getValueAtElement(2, 0);

	double elev = ArcTan(u, sqrt(pow(e, 2) + pow(n, 2)));
	return elev;
}


/**
Returns the azimuth of a satellite
param@ point1, the cartesian coordinates of the station - ECEF (in meters)
param@ point2, the cartesian coordinates of the station - ECEF (in meters)
return@ azimuth of the satellite (in radians)
Leick, Alfred. GPS Satellite Surveying, page 46, (Eq. 2.100)
*/
double GPSModule::azimuthOfSat(Car& point1, Car& point2) {
	double dx = point2.x - point1.x;
	double dy = point2.y - point1.y;
	double dz = point2.z - point1.z;

	double s = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));

	Geo p1 = Car2Geo(point1);

	double term1 = -sin(p1.lon) * dx + cos(p1.lon) * dy;
	double term2 = -sin(p1.lat) *cos(p1.lon) * dx - sin(p1.lat) * sin(p1.lon) * dy + cos(p1.lat) * dz;
	
	double azimuth = ArcTan(term1, term2);

	return azimuth;
}


/**
Returns the Azimuth of a satellite
param@ Coordinates in LGS (NEU) at the station
return elevation angle of the satellite (in radians)
Leick, Alfred. GPS Satellite Surveying, page 45, (Eq. 2.92)
*/
double GPSModule::azimuthOfSat(Matrix& NEU) {
	double n = NEU.getValueAtElement(0, 0);
	double e = NEU.getValueAtElement(1, 0);
	double u = NEU.getValueAtElement(2, 0);

	double azimuth = ArcTan(e, n);
	return azimuth;
}

/**
returns the troposheric path delay based onthe Modified Saastamoinen Troposhere Model (Saastamoinen 1972, 1973) 
param@ height, height of statin in meters
param@ elev, elevation angle of satellite in radians
return@ troposheric path delay in meters

Source : Xu Guochang, GPS Theory, Algorithms and Applications Page 81. 

// DOCU52 **
// Elementsof PPP
// ESA Vol I
*/

double GPSModule::Saastamoinen(double height, double elev) {

//	vector<vector<double>> A;
//	A[0][0] = 1;
	//const double PI = atan(1.0) * 4;
	
	double zenith = PI / 2 - elev;

	const double Rh0 = 0.50; // humidity 
	const double P0 = 1013.25; // standart pressure at height H0 in milibars (mbar)
	const double H0 = 0; // reference height (m)
	const double T0 = CelciusToKelvin(18); // 18 degrees in Celcius is converted to K.
	//cout << "T0 K : " << T0 << endl;

	// T will be transformed to Kelvin here to use in the further steps ....
	// http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/atmos/atmos.html
	//L -- Lapse rate(0.0065 K / m).
	
	double T = T0 - 0.0065 * (height - H0);
	double P = P0 * pow((1 - 0.0000266 * (height - H0)), -5.225);	
	double Rh = Rh0 * exp(-0.0006396 * (height - H0));


	// e: partial pressure of water vapour (in mb).	

	double e = Rh * exp(-37.2465 + 0.213166 * T - 0.000256908 * pow(T, 2)); 
	
	
	double B =  interpolateBValue( height / 1000);
	double dR = 0; // It is negligible...

	double tpd = 0.002277 / cos(zenith) * (P + (1255 / T + 0.05) * e - B * pow(tan(zenith), 2) ) + dR;

	return tpd;
}


double GPSModule::KelvinToCelcius(double kelvin) {
	return (kelvin - 273.15);
}

double GPSModule::CelciusToKelvin(double celcius) {
	return (celcius + 273.15);
}


/**
interpolates and returns the B coefficient value from the table to estimate the tropospheric zebit delay
@param,  height of the station in km.
@return the coeefficient, B, value

Model Source : Xu Guochang, GPS Theory, Algorithms and Applications Page 82.
*/

double GPSModule::interpolateBValue(double x) {
// H(km)		0.0    0.5		1.0		1.5 2.0		2.5	3.0    4.0   5.0
// B(mbar)		1.156 1.079		1.006 0.938 0.874 0.813 0.757 0.654 0.563
	const int numRow = 9;
	const int lastIndex = numRow - 1;
	double table[numRow][2] = {{0.0, 1.156}, {0.5, 1.079}, {1.0, 1.006}, {1.5, 0.938}, 
	                            {2.0, 0.874}, {2.5, 0.813}, {3.0, 0.757}, {4.0, 0.654}, {5.0, 0.563} };

	int loc = 0;
	if (x < table[0][0])
		return table[0][1];
	else if (x >= table[lastIndex][0])
		return table[lastIndex][1];
	else {
		for (int i = 0; i <= lastIndex-1; i++) {
			if ( ( x >= table[i][0] ) && ( x < table[i + 1][0] ) )
			{
				loc = i;
				break;
			}
		 }

		double dy = table[loc + 1][1] - table[loc][1];
		double dx = table[loc + 1][0] - table[loc][0];
		double slope = dy / dx;

		return ( table[loc][1] + slope * ( x - table[loc][0]) );
	}
}


bool GPSModule::sppEpochbyEpoch(BrdcEph & brdc, RinexObsFile & rinexo, ostream& w ) {

	const double C = constants::C;
	Epoch t;
	vector<int> svList;
	vector<double> obsList;
	vector<double> obsListP1;
	vector<double> obsListP2;



	//cout << brdc.getFilename() << endl;
	//cout << rinexo.getFilename() << endl;
	//cout << rinexo.getNumEpoch() << endl;

		
	Car staPos0 = rinexo.getAppCoord(); // Initial coordinate for the site.
	Geo staPos0Geo = Car2Geo(staPos0);
	cout << "Station Position: "; staPos0.print(); cout << endl;
	cin.ignore();
	cin.get();

	for (unsigned int cnt = 0; cnt < rinexo.getNumEpoch(); cnt++) { // for each epoch
		bool isIterated = false;
		t = rinexo.getEpoch(cnt);
		svList = rinexo.getSatSysXAtEpoch(cnt, 'G');
		obsListP1 = rinexo.getAnEpochOfObs(cnt, "P1", 'G');
		obsListP2 = rinexo.getAnEpochOfObs(cnt, "P2", 'G');
		obsList = GPSModule::formP3(obsListP1, obsListP2);
		obsListP1.clear();
		obsListP2.clear();

		double secWeek;
		int GPSWeek, dayOfWeek;
		TimeSystems::DateToGPSTime(t.date.year, t.date.month, t.date.day, t.time.hour, t.time.min, t.time.sec, secWeek, GPSWeek, dayOfWeek);
		
		// Test needed for checking the quality of observations. in case of blunder, it shoud be removed...
		// Do this part later...





		// Calculate the pos. and clk of sat
		vector<Car> svPos(svList.size(), { 0, 0, 0 });
		vector<double> svClk(svList.size(), 0.0);
		vector<double> tropDelay(svList.size(), 0.0);  // setting all the values 0.0 for the initial troposheric zenith delays.
		vector<double> tSatTrue(svList.size());

		// do iterate 2 times from here... later

		for (unsigned int i = 0; i < svList.size(); i++) {
			if (brdc.calcSatPos(svList[i], t, svPos[i])) {
				svClk[i] = (brdc.calcSatClkBias(svList[i], t));
				double elev = GPSModule::elevAngleSat(staPos0, svPos[i]);
				tropDelay[i] = GPSModule::Saastamoinen(staPos0Geo.hgt, elev);
				//cout << i << "  "; svPos[i].print(); 
				//cout << "   SvClk:"; cout << fixed << setprecision(12) << svClk[i] << endl;
				//cout << i << "tropDelay : " << tropDelay[i] << endl;
			}
			else {
				// deleting the data which has no ephemeris ...
				cout << "There is no ephemeris data for Satellite " << i <<"\n This will erase Sv" << i << " " << svList[i] << " from the observation"; // endl;
				svList.erase(svList.begin() + i);
				obsList.erase(obsList.begin() + i);
				cout << svList.size() << " " << obsList.size() << endl;
			}			    
		}

		int iterCounter = 0;
		while (!isIterated){
			iterCounter++;



			// setting the a priori values for receiver clock bias 
			double range = AdjustmentModule::calcRange(staPos0, svPos[0]);
			//const double C = 299792458;
			double dtReceiver = (obsList[0] - range + svClk[0] * C) / C;  // dtReceiver is the receiver clock bias

			vector<double> aprioriUnknown;
			aprioriUnknown.push_back(staPos0.x);
			aprioriUnknown.push_back(staPos0.y);
			aprioriUnknown.push_back(staPos0.z);
			aprioriUnknown.push_back(dtReceiver);
			cout << "a priori size" << aprioriUnknown.size() << endl << endl;



			Matrix X0 = AdjustmentModule::setAprioriMatrix(aprioriUnknown);
			cout << "A priori values for unknowns: \n"; X0.printToScreen();

			Matrix A = AdjustmentModule::setDesignMatrix(staPos0, svPos);
			cout << " Design Matrix : \n"; A.printToScreen();
			//		cin.get();

			Matrix W = AdjustmentModule::setClosureMatrix(aprioriUnknown, svPos, obsList, svClk, tropDelay);
			cout << "Clossure Matrix : \n"; W.printToScreen();
			cin.get();


			// NO adjustment if there are fewer than 4 observations . Go on with the next observation.
			if (obsList.size() < 4) {
				cout << "Not sufficient observation available. Skipping to the next epoch...";
				continue;
			}


			Matrix d = AdjustmentModule::estimateCorrections(A, W);
			cout << "\nA priori unknown corrections" << endl;
			d.printToScreen();
			cout << endl;
			cin.get();

			Matrix xk = X0 + d;
			cout << "\nAdjusted values : " << endl;
			xk.printToScreen();

			if (iterCounter == 2) {
				isIterated = true; // bunlardan biri gereksiz
				if (isIterated) {
					cout << endl << "Final Solution for epoch : " << cnt << "  ";  t.print();
					cout << endl << "Final Coordinate: " << endl; xk.printToScreen();
					Car p;
					p.x = xk.getValueAtElement(0, 0);
					p.y = xk.getValueAtElement(1, 0);
					p.z = xk.getValueAtElement(2, 0);
					Matrix dneu = GPSModule::xyz2neu(staPos0, p);
					cout << endl << "NEU: " << endl; dneu.printToScreen();
					cin.get();
				}
				break;
			}
				



			double tReceiverTrue = secWeek - xk.getValueAtElement(3, 0); // In GPS time...

			// lightTimeIteration 

			tSatTrue = lightTimeIteration(svList, GPSWeek, tReceiverTrue, brdc, staPos0);
			cout << "End of Test-lightTimeIteration" << endl;
			cin.get();

			// Recalculcating the positions of the Satellites with tSatTrue
			for (unsigned int i = 0; i < svList.size(); i++) {
				Matrix R3(3, 3);
				//double tFly = tReceiverTrue - tSatTrue[i];
				double WT = constants::OMEGAE_DOT * (tReceiverTrue - tSatTrue[i]);

				R3.setValueAtElement(0, 0, cos(WT));
				R3.setValueAtElement(0, 1, sin(WT));
				R3.setValueAtElement(1, 0, -sin(WT));
				R3.setValueAtElement(1, 1, cos(WT));
				R3.setValueAtElement(2, 2, 1);
				//cout << "R3 Matrix" << endl;
				R3.printToScreen();
				if (brdc.calcSatPos(svList[i], GPSWeek, tSatTrue[i], svPos[i])) {

					Matrix S(3, 1);
					S.setValueAtElement(0, 0, svPos[i].x);
					S.setValueAtElement(1, 0, svPos[i].y);
					S.setValueAtElement(2, 0, svPos[i].z);
					Matrix R3S = R3 * S;
					Car S_F;
					S_F.x = R3S.getValueAtElement(0, 0);
					S_F.y = R3S.getValueAtElement(1, 0);
					S_F.z = R3S.getValueAtElement(2, 0);
					svPos[i] = S_F;


					svClk[i] = (brdc.calcSatClkBias(svList[i], GPSWeek, tSatTrue[i]));

					/// Tropospehere
					double elev = GPSModule::elevAngleSat(staPos0, svPos[i]);
					double tropD = GPSModule::Saastamoinen(staPos0Geo.hgt, elev);
					tropDelay[i] = tropD;

				}
			} // Recalculcating the positions of the Satellites with tSatTrue

		} // iteration





	} // for each epoch



	return true;
}

	vector<double> GPSModule::lightTimeIteration(const vector<int>& svList, const int& GPSWeek, 
		                      const double& tReceiverTrue, BrdcEph& brdc, Car& appPos) {		
		
		const double C = constants::C;
		const int NUM_ITER = 10;
		vector<double> tSatTrue;// in GPStime
		vector<double> tTemp;
		

		for (int i = 0; i < svList.size(); i++) {
			tTemp.push_back(tReceiverTrue);
			Car satPos0 = { 0,0,0 };
			//cout << "Satellite # : " << svList[i] << endl;
			for (int j = 1; j < NUM_ITER; j++) {
				brdc.calcSatPos(svList[i], GPSWeek, tTemp[j - 1], satPos0);
				cout << satPos0.x << " " << satPos0.y << " " << satPos0.z << endl;
				double range = AdjustmentModule::calcRange(appPos, satPos0);
				double tt = tReceiverTrue - range / C;
				cout << "True Sat Time" << tt << endl << endl;
				tTemp.push_back(tt);
			}
			cout << endl << endl << endl;
			//double zz = tTemp[NUM_ITER-1];// .pop_back();
			tSatTrue.push_back(tTemp[NUM_ITER - 1]);
			//tSatTrue.push_back(tTemp.pop_back());
			tTemp.clear();

		}

		return tSatTrue;
	}
