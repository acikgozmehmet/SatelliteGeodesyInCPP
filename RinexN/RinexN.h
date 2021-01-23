#include<iostream>
#include<string>
#include<fstream>
#include<stdexcept>
#include<iomanip>
#include<cmath>

#include "..\GPSModule\GPSModule.h"
#include "..\TimeSystems\TimeSystems.h"

using namespace std;


struct Ephemeris {
	int svprn;
	double af0, af1, af2;
	double crs, deltan, M0;
	double cuc, ecc, cus, roota;
	double toe, cic, Omega0, cis;
	double i0, crc, omega, Omegadot;
	double idot, tom;

	bool hasData;

	Ephemeris() {
		hasData = false;
	}
};


// Maximum number of Satellite Vehicles (SV) + 1
const int MAX_SV = 32 + 1;

// Maximum number of epochs for each SV in the file. 
const int MAX_EPOCH = 13; // max number of epochs 

class RinexN {
public:
	RinexN();
	RinexN(string filename);
	~RinexN();


	void setFilename(string filename);
	string getFilename(void) const;
	void loadRinexN();
  
	bool calcSatPos(int svprn, double tObs, CarCoord& satPos);

	double calcSatClkBias(int svprn, double tObs);  // This has to be implemented...
	//friend CarCoord calcSatPos(int svprn, double tObs);

private:

	double m_toeStart;
	string m_filename;
	Ephemeris m_data[MAX_SV][MAX_EPOCH];
	//ofstream outFile; // for output


// Bernsedeki gibi out dosyatisndakinuzantilari artiracak sekide bir static int tanimlanabili ... L01. etc
	 
	void readFromFile(const string& filename);
	int hashFunction(double toe);	
	int getEphIndex(double timeOfObs);
	Ephemeris getEphData(int svprn, int ephIndex);
	void printGraphics(ostream& w);

	/**
	Changes the Fortran Scientific notation for double numbers (D-, D+) to 
	C++ notation (E-, E+)
	@param text, a line of string which contains only numbers 
	*/
	void changeSciNotation(string& text);	
	/*
	To check the time 
	*/
	void checkTk(double& tk);
	/*
	returns the ArcTan value 
	*/
	double ArcTan(double x, double y); 

	/* Returns the 2 digit year to 4 digit year 
	@param, 2 digit year yy
	@return 4 digit year yy
	*/
	void convertYYToYYYY(int & yy);
	void printLoadingBar(int counter);


};

void RinexN::loadRinexN() {
	cout << m_filename << " Loading ..." << endl;
	readFromFile(m_filename);
}

RinexN:: RinexN(){
	setFilename("");
}

RinexN::RinexN(string filename){
	setFilename(filename);
}
RinexN::~RinexN() {
}

void RinexN::setFilename(string filename) {
	m_filename = filename;
}

string RinexN::getFilename(void) const {
	return m_filename;
}

void RinexN::readFromFile(const string& filename) {
	ifstream inFile;
	string line;

	try {
		inFile.open(filename);
		if (!inFile)
			throw runtime_error(filename + " does not exist !");	
	}
	catch (runtime_error& ex) {
		cout << endl << "Error ! : " << ex.what() << endl;
		exit(1);
	}

	// For the default output file...
	string outFilename = filename.substr(0, filename.find(".")) + ".out";
	ofstream outFile;
	outFile.open(outFilename);



	// Getting the name of the brdc file and writing it to outFile
	string stationName = filename.substr(0, 4);
	for (unsigned int i = 0; i < stationName.length(); i++)
		stationName.at(i) = toupper(stationName.at(i));
	outFile << "\t\t\t" <<stationName << endl;

	// getting the date (year and doy) from the filename
	int year = atoi(filename.substr(9, 2).c_str());
	// to convert the 2 digit year to 4 digit year (the following logic is valid until 2080)
	convertYYToYYYY(year);
	int doy = atoi(filename.substr(4, 3).c_str());

	int month, day, hour = 0, minute = 0;
	double second = 0;
	TimeSystems::DoyToDate(year, doy, month, day, hour, minute, second); 
    

	
	int gpsweek, dow;
	TimeSystems::DateToGPSTime(year, month, day, hour, minute, second, m_toeStart, gpsweek, dow);
	outFile << endl << "Rinex Navigation File Information  : "<< endl;
	outFile << "Date        : " << setw(2) << month << "/" << setw(2) << day << "/" << year << endl;
	outFile << "GPS week    : " << setw(10) << gpsweek << endl;
	outFile << "Day of Week : " << setw(10) << dow << endl;
	outFile << "GPS time    : " << setw(10) << m_toeStart << endl << endl;

	// this part of code is to skip the header part of the file
	do {
		getline(inFile, line);
	} while (line.substr(60, 13) != "END OF HEADER");

	int counter = 0;
	while ( getline(inFile, line) && line.size() > 0) {

		changeSciNotation(line); // to get rid of the fortran style double sign D+2 -> E+2

		Epoch toc;
		// first line of the nav. message : prn, epoch, svn clk
		int prn = atoi(line.substr(0, 2).c_str());
		toc.yy = atoi(line.substr(3, 2).c_str());
		toc.mm = atoi(line.substr(6, 2).c_str());
		toc.dd = atoi(line.substr(9, 2).c_str());
		toc.hh = atoi(line.substr(12, 2).c_str());
		toc.min = atoi(line.substr(15, 2).c_str());
		toc.sec = atof(line.substr(17, 5).c_str());
		convertYYToYYYY(toc.yy);

		Ephemeris eph;
		eph.svprn = prn;

		eph.af0 = atof(line.substr(22, 19).c_str());
		eph.af1 = atof(line.substr(41, 19).c_str());
		eph.af2 = atof(line.substr(60, 19).c_str());

		// broadcast orbit 1
		getline(inFile, line);
		changeSciNotation(line);

		//    double IODE, crs, deltan, M0;
		//IODE = atof(line.substr(3,22).c_str());
		eph.crs = atof(line.substr(22, 19).c_str());
		eph.deltan = atof(line.substr(41, 19).c_str());
		eph.M0 = atof(line.substr(60, 19).c_str());

		// broadcast orbit 2
		getline(inFile, line);
		changeSciNotation(line);

		//    cuc, ecc, cus, roota;
		eph.cuc = atof(line.substr(3, 19).c_str());
		eph.ecc = atof(line.substr(22, 19).c_str());
		eph.cus = atof(line.substr(41, 19).c_str());
		eph.roota = atof(line.substr(60, 19).c_str());


		// broadcast orbit 3
		getline(inFile, line);
		changeSciNotation(line);

		//    toe, cic, Omega0, cis;
		eph.toe = atof(line.substr(3, 19).c_str());
		eph.cic = atof(line.substr(22, 19).c_str());
		eph.Omega0 = atof(line.substr(41, 19).c_str());
		eph.cis = atof(line.substr(60, 19).c_str());


		// broadcast orbit 4
		getline(inFile, line);
		changeSciNotation(line);

		//    i0, crc, omega, Omegadot;
		eph.i0 = atof(line.substr(3, 19).c_str());
		eph.crc = atof(line.substr(22, 19).c_str());
		eph.omega = atof(line.substr(41, 19).c_str());
		eph.Omegadot = atof(line.substr(60, 19).c_str());


		// broadcast orbit 5
		getline(inFile, line);
		changeSciNotation(line);

		//   idot, codes, weekno, L2flag;
		eph.idot = atof(line.substr(3, 19).c_str());
		//codes = atof(line.substr(22,41).c_str());
		//weekno = atof(line.substr(41,60).c_str());
		//L2flag = atof(line.substr(60,79).c_str());


		// broadcast orbit 6
		getline(inFile, line);
		changeSciNotation(line);

		//   svaccur, svhealth, tgd, iodc;
		/*   svaccur = atof(line.substr(3,22).c_str());
		svhealth = atof(line.substr(22,41).c_str());
		tgd = atof(line.substr(41,60).c_str());
		iodc = atof(line.substr(60,79).c_str());
		*/

		// broadcast orbit 7
		getline(inFile, line);
		changeSciNotation(line);
		//cout << line << endl;


		//   tom, spare, spare, spare;
		eph.tom = atof(line.substr(3, 19).c_str());
		//cout << eph.tom << endl;


		eph.hasData = true;
		int index = hashFunction(eph.toe);
		outFile << "Indexing SV # " <<setw(2) << eph.svprn <<" at the epoch " << eph.toe <<"....index : " << setw(2) << index << endl;


		if (index >= 13) {
			outFile << "FATAL ERROR while indexing !!!!: SV# " <<eph.svprn << "  " << index << endl;
		}

		
		if ( ! m_data[prn][index].hasData ) {
			m_data[prn][index] = eph;
		}
		else
		{
			outFile << "\t*Collision for the SV # " << setw(2) << m_data[prn][index].svprn << " and epoch : " << m_data[prn][index].toe ;
			//outFile << "\t*Epoch of previous data and new data respectivley : " << m_data[prn][index].toe << "  " <<eph.toe << endl;
			if (eph.toe > m_data[prn][index].toe) {
				m_data[prn][index] = eph;				
			    outFile << "  Indexing updated with the most recent one" <<endl;
			}
		}
		// Reporting the status
		counter++;
		printLoadingBar(counter);
	}
	cout << " 100 %" << endl;  // for printLoadingBar to display that all the data is loaded...
	printGraphics(outFile);
	inFile.close();		
	outFile.close();
}


int RinexN::hashFunction(double toe) {
	const double INTERVAL = 2 * 60; // 2 minutes interval
	const double UPDATE_INTERVAL = 2 * 60 * 60;// 2 hours 
	double toeRef;

	for (int i = 0; i < MAX_EPOCH; i++) {	

		toeRef = m_toeStart + i * UPDATE_INTERVAL; // increasing the frame each time for 2 hours 	
		
		if ( (toe >= (toeRef - INTERVAL)) && (toe <= (toeRef + INTERVAL)) ) {
			return i;
		}
	}
	cout << "No index was found for this toe." << endl;
	return -1;
}

/**
Returns the index for the Ephemeris that will be used in the calculations
@param, time in GPStime (secOfWeek)
@return index number of found or -1 if not found
*/
int RinexN::getEphIndex(double timeOfObs) {	
	const double INTERVAL = 1 * 60; // 1 minutes interval

	return ( (timeOfObs + INTERVAL) >= m_toeStart) ? (timeOfObs + INTERVAL - m_toeStart) / 7200 : -1;
}

void RinexN::printGraphics(ostream& w) {
	w << endl << "\t\t\t\tRinex Navigation File Data Availibilty Graph" << endl<< endl;
	w << "       ";

	for (int i = 1; i < MAX_SV; i++)
		w << setw(3) << i;// % 10;
	w << endl;

	for (int i = 0; i < MAX_EPOCH; i++) {
		w << setw(2) << (2 * i) % 24 << " - " << setw(2) << (2 * i + 2)  % 24 ; // << "  ";
		for (int j = 1; j < MAX_SV; j++) {
			w << (m_data[j][i].hasData ? "  X" : "   ");
		}
		w << endl;
	}
}


void RinexN::changeSciNotation(string& text) {
	for (unsigned int i = 0; i < text.length(); i++) {
		if (toupper(text.at(i) == 'D') ) 
			text.at(i) = 'E';
	}
}


Ephemeris RinexN::getEphData(int svprn, int ephIndex) {
	Ephemeris temp;

	if (m_data[svprn][ephIndex].hasData) {
		temp = m_data[svprn][ephIndex];
	}
	else if (m_data[svprn][ephIndex + 1].hasData) {
		temp = m_data[svprn][ephIndex + 1];
	}
	
	return temp;	
}




void RinexN::checkTk(double& tk) {
	if (tk > 302400)
		tk -= 604800;
	else if (tk < -302400)
		tk += 604800;
}


double RinexN::ArcTan(double x, double y) {
	const double PI = std::atan(1.0) * 4;
	const double RHO = 180.0 / PI;
	double teta;
	teta = atan2(x, y);
	teta < 0 ? teta = teta + 2 * PI : teta = teta;
	return teta;
}


/*
GPS data processing: code and phase
Algorithms, Techniques and Recipes
Research group of Astronomy and GEomatics (gAGE/UPC)
*/

/**
Returns true if the Satellite Coordinate is calculated at any given time, or else false
@param svprn, Satellite #
@param tObs, time (second of GPSWeek) at which Sat. Coordinates are to be calculated. 
@param, Satellite Position to calculate
*/

bool RinexN::calcSatPos(int svprn, double tObs, CarCoord& satPos) {

	satPos = { 0,0,0 };
	const double PI = std::atan(1.0) * 4;

    //const double GM = 3.986005e14;  // earth's universal gravitational parameter m^3/s^2
	const double  GM = 3.986004418e14;
	const double  Omegae_dot = 7.2921151467e-5; // % earth rotation rate, rad/s

	/*
	% Earth's universal gravitational parameter, m^3/s^2
	GM = 3.986004418e14;

	% earth rotation rate, rad/s
	Omegae_dot = 7.2921151467e-5;

	% speed of light
	c = 0.299792458e9;


	BERNESE CONST
	C      =  299792458.D0    VELOCITY OF LIGHT                M/SEC
	GM     = 398.6004415D12   GRAVITY CONSTANT*EARTH MASS      M**3/SEC**2
	OMEGA  = 7292115.1467D-11 ANGULAR VELOCITY OF EARTH        RAD/SEC
	FROM INTERNET
	IERS 3.986004418×10-14 7.2921150×10-5
	*/


	// First step to get the ephemerides ...
	int ephIndex = getEphIndex(tObs);
	if (ephIndex == -1) {
		cout << "Error! :No Solution available since lack of ephemeris data for the epoch  " << tObs << endl;
		cin.ignore();
		cin.get();
		exit(1);
	}
	Ephemeris eph = getEphData( svprn, ephIndex);
	if ( (svprn != eph.svprn) && (eph.hasData)) {
		cout << svprn << " ? " << eph.svprn << endl;
		cout << "Error! :Probable error in indexing the ephemeris "<< endl;
		cin.ignore();
		cin.get();
		exit(1);
	}

	if (!eph.hasData) { // if there is no ephemeris for 4 hour period , no use calculating....
		return false;
	}

		// tk
		double tk = tObs - eph.toe;
		checkTk(tk);

		// mean anomaly for tk
		double A = eph.roota * eph.roota;
		double M = eph.M0 + (sqrt(GM / pow(A, 3)) + eph.deltan) * tk;
		//cout <<" A= " << A << "  " << M<< endl;

		// Solving (iteratively) Kepler equation for the eccentricity anomaly Ek:
		double E, Eold, dE, v;
		E = M;
		do {
			Eold = E;
			E = M + eph.ecc * sin(E);
			dE = abs(Eold - E);
		} while (dE > 1e-12);

		//Calculation of real anomaly vk:
		double uu1 = sqrt(1 - pow(eph.ecc, 2)) * sin(E);
		double uu2 = (cos(E) - eph.ecc);
		v = ArcTan(uu1, uu2);



		// Calculation of the argument of latitude uk from the argument of
		// perigee ω, real anomaly vk and corrections cuc and cus:
		double phi = v + eph.omega;
		phi = fmod(phi, (2 * PI));
		double u = phi + eph.cuc*cos(2 * phi) + eph.cus*sin(2 * phi);

		//Calculation of the radial distance rk, considering corrections crc and crs:
		double r = A * (1 - eph.ecc*cos(E)) + eph.crc*cos(2 * phi) + eph.crs*sin(2 * phi);

		// Calculation of inclination ik of the orbital plane from the inclination
		// io at the reference time toe, and corrections cic and cis:
		double i = eph.i0 + eph.idot*tk + eph.cic*cos(2 * phi) + eph.cis*sin(2 * phi);

		/*Calculation of the longitude of the ascending node λk (referring
		to Greenwich), using its right ascension
		o at the beginning
		of the current week, corrected from the apparent sidereal time
		variation in Greenwich between the beginning of the week and
		reference time tk = t − toe, and the change in longitude of the
		ascending node from the reference time toe.
		*/

		double Omega = eph.Omega0 + (eph.Omegadot - Omegae_dot)*tk - Omegae_dot * eph.toe;
		Omega = fmod((Omega + 2 * PI), (2 * PI));


		// Position in orbital plane
		double x1, y1;
		x1 = cos(u) * r;
		y1 = sin(u) * r;

		/*
		Calculation of coordinates in CTS frame, by applying three rotations
		(about uk, ik, λk):
		*/
		satPos.x = x1 * cos(Omega) - y1 * cos(i)*sin(Omega);
		satPos.y = x1 * sin(Omega) + y1 * cos(i)*cos(Omega);
		satPos.z = y1 * sin(i);

		return true;
}


void RinexN::convertYYToYYYY(int & yy) {
	if (yy < 100) {
		if (yy < 80) {
			yy += 2000;
		}
		else
			yy += 1900;
	}
}

/**
Returns the Satellite clock bias
param svprn, Satellite number
param tObs, observation time in GPStime.
*/

double RinexN::calcSatClkBias(int svprn, double tObs) {
	//af0+ af1*t + af2* t^2      t -> tk
	int index = getEphIndex(tObs);
	double tk = tObs - m_data[svprn][index].toe; //% time elapsed since toe
	return m_data[svprn][index].af0 + (m_data[svprn][index].af1 * tk) + (m_data[svprn][index].af2 *pow(tk, 2));
}


void RinexN::printLoadingBar(int counter) {

	int refVal = round(double(MAX_EPOCH * (MAX_SV-1) ) / 100) * 2;

	if (counter == 0)
		cout << "Loading Navigation File...  " << m_filename << endl;
	else if (counter % refVal == 0)
		cout << "*";

}
