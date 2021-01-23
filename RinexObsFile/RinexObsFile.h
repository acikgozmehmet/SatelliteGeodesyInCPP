#ifndef  RINEXOBSFILE_H
#define RINEXOBSFILE_H

#include<iostream>
#include<string>
#include<fstream>
#include<stdexcept>
#include<iomanip>
#include<cmath>
#include<vector>

#include "..\Constants\Constants.h"


struct ObsData {
	Epoch toc;
	int flag;
	int numSatInCurEpoch;
	string * svIDInCurEpoch = nullptr;

	double ** data = nullptr;
	int ** lli = nullptr; // Loss of lock indicator (LLI)
	int ** ss = nullptr;  //Signal strength

	ObsData() {
		numSatInCurEpoch = 0;
		flag = 0;
	}
};

class RinexObsFile {
public:
	RinexObsFile();
	RinexObsFile(string filename);
	~RinexObsFile();

	void setFilename(string filename);
	string getFilename();
	void loadFromFile();

	int getNumEpoch();
/**
   Returns the epoch of the observation
   @param atIndex is the index of the observation
   @return the epoch
*/
	Epoch getEpoch(int atIndex);

	/*
	Returns the List of the SV in the SatType at current epoch
	@param atIndex is the index of the epoch
	@param satType is the Satellite System indicator ( G = NAVSTAR-GPS (default), R GLONASS etc)
	@return the List of the SV in the SatType at current epoch
	*/
	vector<int> getSatSysXAtEpoch(int atIndex, const char& satType = 'G');

/**
	Returns the desired observation at the current epoch
	@param atIndex, index of the obervation epoch
	@param satType, Satellite Type desired. ( G = NAVSTAR-GPS,R GLONASS etc)
	@param observationType, desired observation type. (Ex: L1, L2, C1, C2, P1, P2 etc)
	#return the desired observation at the current epoch
*/

	vector<double> getAnEpochOfObs(int atIndex, string observableType, const char& satType = 'G' );

	Car getAppCoord();


private:
	string m_filename; 
	int m_numEpoch;                                      // number of epochs in the file.
	int m_countEpoch;                                    // counter for epoch

// data members in the header part
	double m_version;
	char m_fileType, m_satSys;

	string m_markerName, m_markerNumber;
	string m_recNum, m_recType, m_recVers;
	string m_antNum, m_antType;
	Car m_appPos;                                       // approcimate position of the site.
	double m_antHgtH, m_antHgtE, m_antHgtN;             
	int m_waveLengthFacL1, m_waveLengthFacL2;
	int m_numObservables;                               // number of observables in the file.
	string *m_observables = nullptr;                    // the observables in the file.
	int m_interval;                                     // interval of the observation
	Epoch m_timeFirstObs, m_timeLastObs;
	string m_timeFirstObsTimeSys, m_timeLastObsTimeSys;
	int m_leapSecond;

	// data member for the obs data
	ObsData * obsList = nullptr;



	// helper methods
	bool readHeaderPart(ifstream & inFile);
	void readObsPart(ifstream & inFile);

/** Converts  the 2 digit year to 4 digit year format.
	@param yy is the 2 digit year.
	@return None.
*/
	void convertYYToYYYY(int & yy);

	void printLoadingBar(int counter);

/**
	Returns the indices of the desired Satellite System (ie. GPS (G), GLONASS (R) etc) at the current epoch.
	@param atIndex is the index for the each epoch for obsList.
	@param satType is Satellite System indicator, ( like G: NAVSTAR-GPS (default), R: GLONASS)
	@return the indices of the desired Satellite System Satellites.
*/
	vector<int> getSatSysXIndicesAtEpoch(int atIndex, const char& satType = 'G');


/** 
    Returns the index of the observable in the file
    @param observable is the type of the desired observable from the user (ie. L1, P1, L2 etc).
	@return returns the index if available, otherwise returns -1.
*/
	int getIndexOfObservable(const string& observable);





	

};




#endif // ! RINEXOBSFILE_H


