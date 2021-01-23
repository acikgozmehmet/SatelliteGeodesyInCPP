#include<iostream>
#include "RinexN.h"

//#include "GPSGeodesy.h"

using namespace std;

int main() {

	RinexN A;
 	A.setFilename("epgg0100.02n");
	A.loadRinexN();
	cout << "All done .." << endl;


	int prn = 31;
	double time = 346500;
	CarCoord satPos;
	cout << "Position of the Satellite Vehicle #" << setw(2) << prn << "  ";


	A.calcSatPos(prn, time,satPos);
	//x = calcSatPos(prn, time);
	//x = calcSatPos(A, prn, time);
	cout << std::fixed << setprecision(4);
	cout << satPos.x << "  " << satPos.y << "   " << satPos.z << endl;
	cout << "Sat clk bias : " << setprecision(10) << A.calcSatClkBias(31, time);




	//calcSatPos(1, 418800);
	//cout << "Printing ...." << endl;
	//cout << "Saving to a file...." << endl;
	//A.saveToFile();

	
	//RinexO obs("ankr0200.17o");
	//obs.skipHeader();
	//string line = obs.getData();
	//cout << line << endl;
	






	cin.ignore();
	cin.get();




	return 0;
}