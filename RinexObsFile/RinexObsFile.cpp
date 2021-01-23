#include "RinexObsFile.h"

RinexObsFile::RinexObsFile(){}
RinexObsFile::RinexObsFile(string filename):m_filename(filename){ }
RinexObsFile::~RinexObsFile(){
	for (unsigned int i = 0; i < m_numEpoch; i++) {
		for (unsigned int j = 0; j < obsList[i].numSatInCurEpoch; j++) {
			delete[] obsList[i].data[j];
			delete[] obsList[i].lli[j];
			delete[] obsList[i].ss[j];
			//delete[] obsList[i].svIDInCurEpoch[j];
		}
		delete[] obsList[i].data;
		delete[] obsList[i].lli;
		delete[] obsList[i].ss;
		delete[] obsList[i].svIDInCurEpoch;

		obsList[i].data = nullptr;
		obsList[i].lli = nullptr;
		obsList[i].ss = nullptr;
		obsList[i].svIDInCurEpoch = nullptr;
	}
	// Still not sure if obs has been dealloacted ?? Do I have to use sth like "delete obs" ?? When I do it, I have runtime error
	delete[] m_observables;
	m_observables = nullptr;
	obsList = nullptr;
}

void RinexObsFile::setFilename(string filename) { m_filename = filename; }
string RinexObsFile::getFilename() { return m_filename; }

int RinexObsFile::getNumEpoch() { return m_countEpoch; }

Epoch RinexObsFile::getEpoch(int atIndex) { return obsList[atIndex].toc; }


void RinexObsFile::loadFromFile() {
	ifstream inFile;
	inFile.open(m_filename);
	
	if (!inFile) {
		cout << "Error: " << m_filename<< "does NOT exist !" << endl;
		cin.get();
		exit(1);
	}
	else {

		readHeaderPart(inFile); // first the header 
		// do some math and calculate the number for the m_numEpoch
		m_numEpoch = (24 * 60 * 60) / m_interval;

		obsList = new ObsData[m_numEpoch];

		readObsPart(inFile);

		//cin.get();

	}



}


bool RinexObsFile::readHeaderPart(ifstream & inFile){
	string line;
	do {
		getline(inFile, line);

		if (line.substr(60, 20) == "RINEX VERSION / TYPE") {
			m_version = atof(line.substr(0, 9).c_str());
			m_fileType = line.at(20);
			m_satSys = line.at(40);
		}
		else if ( line.substr(60,11) == "MARKER NAME" ) {
			m_markerName = line.substr(0, 4).c_str();
		}
		else if ( line.substr(60,13) == "MARKER NUMBER" ) {
			m_markerNumber = line.substr(0, 20).c_str();
		}
		else if ( line.substr(60,19) == "REC # / TYPE / VERS" )  {
			m_recNum = line.substr(0, 19).c_str();
			m_recType = line.substr(20, 19).c_str();
			m_recVers = line.substr(40, 19).c_str();
		}
		else if (line.substr(60,12) == "ANT # / TYPE") {
			m_antNum = line.substr(0, 19).c_str();
			m_antType = line.substr(20, 19).c_str();			
		}
		else if ( line.substr(60,19) == "APPROX POSITION XYZ") {
			m_appPos.x = atof(line.substr(0, 14).c_str());
			m_appPos.y = atof(line.substr(14, 14).c_str());
			m_appPos.z = atof(line.substr(28, 14).c_str());
		}
		else if (line.substr(60, 20) == "ANTENNA: DELTA H/E/N") {
			 m_antHgtH = atof(line.substr(0, 14).c_str());
			 m_antHgtE = atof(line.substr(14, 14).c_str());
			 m_antHgtN = atof(line.substr(28, 14).c_str());
		}
		else if ( line.substr(60,19) == "WAVELENGTH FACT L1/2" ) {
			m_waveLengthFacL1 = atoi(line.substr(5,1).c_str());
            m_waveLengthFacL2 = atoi(line.substr(11,1).c_str());
		}
		else if ( line.substr(60,19) == "# / TYPES OF OBSERV" ) {
			m_numObservables = atoi(line.substr(0, 6).c_str());
			m_observables = new string[m_numObservables];
			int numObsPerLine = 9;

			int numLineToRead = 1;
			int numObs = m_numObservables;

			int i = 0;
			while (numObs > numObsPerLine) {
				for (unsigned int j = 0; j < numObsPerLine; j++)
					m_observables[i * numObsPerLine + j] = line.substr(10 + 6 * j, 2);

				getline(inFile, line);
				i++;
				numObs -= numObsPerLine;
			}
			for (unsigned int j = 0; j < numObs; j++) // for the remaining part which is fewer than 9
				m_observables[i * numObsPerLine + j] = line.substr(10 + 6 * j, 2);					
			
			// for test
			//for (int k = 0; k < m_numObservables; k++)  
				//cout << k << " " << m_observables[k] << endl;
		} //  else if "# / TYPES OF OBSERV" 
		else if ( line.substr(60, 8) == "INTERVAL" ) {
			m_interval = atof(line.substr(0, 10).c_str());
		}
		else if (line.substr(60, 12) == "LEAP SECONDS") {
			m_leapSecond = atoi(line.substr(0, 6).c_str());
		}
		else if ( line.substr(60, 17) == "TIME OF FIRST OBS" || line.substr(60, 16) == "TIME OF LAST OBS") {
			Epoch temp;
			string tempSys;

			temp.date.year  = atoi(line.substr(0, 6).c_str());
			temp.date.month = atoi(line.substr(6, 6).c_str());
			temp.date.day   = atoi(line.substr(12, 6).c_str());
			temp.time.hour  = atoi(line.substr(18, 6).c_str());
			temp.time.min   = atoi(line.substr(24, 6).c_str());
			temp.time.sec   = atof(line.substr(30, 13).c_str());
			tempSys = line.substr(48, 3).c_str();

			if (line.substr(60, 17) == "TIME OF FIRST OBS") {
				m_timeFirstObs = temp;
				m_timeFirstObsTimeSys = tempSys;
			}
			else {
				m_timeLastObs = temp;
				m_timeLastObsTimeSys = tempSys;
			}
		}// end of if cases
	}while (line.substr(60, 13) != "END OF HEADER");

	return true;  // if successfully read
}


void RinexObsFile::readObsPart(ifstream & inFile){
	string line;
	m_countEpoch = 0;
	getline(inFile, line);
	while ( !inFile.eof() ) {
		obsList[m_countEpoch].toc.date.year = atoi(line.substr(1, 2).c_str()); 
		convertYYToYYYY(obsList[m_countEpoch].toc.date.year);
		obsList[m_countEpoch].toc.date.month = atoi(line.substr(3, 3).c_str());
		obsList[m_countEpoch].toc.date.day = atoi(line.substr(6, 3).c_str());
		obsList[m_countEpoch].toc.time.hour = atoi(line.substr(9, 3).c_str());
		obsList[m_countEpoch].toc.time.min = atoi(line.substr(12, 3).c_str());
		obsList[m_countEpoch].toc.time.sec = atof(line.substr(15, 11).c_str());
		obsList[m_countEpoch].flag = atoi(line.substr(26, 3).c_str());
		obsList[m_countEpoch].numSatInCurEpoch = atoi(line.substr(29, 3).c_str());

		if (obsList[m_countEpoch].flag != 0) { // If this epoch of obs is not good, just go on with the next one.
			for (int i = 0; i < obsList[m_countEpoch].numSatInCurEpoch + 1; i++)
				getline(inFile, line);
			continue;
		}

		// Storing the ids of the SV...
		obsList[m_countEpoch].svIDInCurEpoch = new string[obsList[m_countEpoch].numSatInCurEpoch];

		const int MAX_NUM_SAT_IN_ROW = 12;
		int k = 0;
		int numSatCount = obsList[m_countEpoch].numSatInCurEpoch;
		//for (int i = 0; i < obsList[m_countEpoch].numSatInCurEpoch; i++)
		//	cout << obsList[m_countEpoch].svIDInCurEpoch[i] << endl;


		while (numSatCount > MAX_NUM_SAT_IN_ROW) {			
			for (unsigned int i = 0; i < MAX_NUM_SAT_IN_ROW; i++) {
				obsList[m_countEpoch].svIDInCurEpoch[k * MAX_NUM_SAT_IN_ROW + i] = line.substr(32 + i * 3, 3);
				//cout << i << " " << obsList[m_countEpoch].svIDInCurEpoch[k * MAX_NUM_SAT_IN_ROW + i] << endl;
			}
			getline(inFile, line);
			k++;
			numSatCount -= MAX_NUM_SAT_IN_ROW;
		}
		for (unsigned int i = 0; i < numSatCount; i++) {
			obsList[m_countEpoch].svIDInCurEpoch[k * MAX_NUM_SAT_IN_ROW + i] = line.substr(32 + i * 3, 3);
			//cout << i << " " << obsList[m_countEpoch].svIDInCurEpoch[k * MAX_NUM_SAT_IN_ROW + i] << endl;
		}
			
		// creating the matrix for the obs data
		obsList[m_countEpoch].data = new double *[obsList[m_countEpoch].numSatInCurEpoch];
		for (int i = 0; i < obsList[m_countEpoch].numSatInCurEpoch; i++) {
			obsList[m_countEpoch].data[i] = new double[m_numObservables];
		}

		//  Loss of lock indicator (LLI)
		obsList[m_countEpoch].lli = new int *[obsList[m_countEpoch].numSatInCurEpoch];
		for (int i = 0; i < obsList[m_countEpoch].numSatInCurEpoch; i++) {
			obsList[m_countEpoch].lli[i] = new int[2];
		}
		//  Signal strength
		obsList[m_countEpoch].ss = new int *[obsList[m_countEpoch].numSatInCurEpoch];
		for (int i = 0; i < obsList[m_countEpoch].numSatInCurEpoch; i++) {
			obsList[m_countEpoch].ss[i] = new int[2];
		}

		// reading one epoch each time
		const int MAX_REC_PERLINE = 5;
		for (int i = 0; i < obsList[m_countEpoch].numSatInCurEpoch; i++) {
			getline(inFile, line);
			for (int j = 0; j < m_numObservables; j++) {
				if ( (j > 0) && ((j % MAX_REC_PERLINE) == 0) )
					getline(inFile, line);

				if (line.length() < 80) { // to fill the line with whitespace when its length is fewer than 80 characters
					for (int i = line.length(); i <= 80; i++)
						line.insert(i, " ");
				}

				obsList[m_countEpoch].data[i][j] = atof(line.substr(16 * (j % MAX_REC_PERLINE), 14).c_str());

				// LLI and SS part
				if (m_observables[j] == "L1" || m_observables[j] == "L2") {
					int val_lli = atoi(line.substr(14 + 16 * (j % MAX_REC_PERLINE), 1).c_str());
					int val_ss = atoi(line.substr(15 + 16 * (j % MAX_REC_PERLINE), 1).c_str());
					if (m_observables[j] == "L1") {
						obsList[m_countEpoch].lli[i][0] = val_lli;
						obsList[m_countEpoch].ss[i][0] = val_ss;
					}
					else
					{
						obsList[m_countEpoch].lli[i][1] = val_lli;
						obsList[m_countEpoch].ss[i][1] = val_ss;
					}
				} // end of LLI and SS part
			} // for (int j = 0; j < m_numObservables; j++) 
		} // for (int i = 0; i < obsList[m_countEpoch].numSatInCurEpoch; i++)
	printLoadingBar(m_countEpoch);
	m_countEpoch++;
	getline(inFile, line);
	}
	
	cout << " 100 %" << endl;  // for printLoadingBar to display that all the data is loaded...
	inFile.close();

}


void RinexObsFile::convertYYToYYYY(int & yy) {
	if (yy < 100) {
		if (yy < 80) {
			yy += 2000;
		}
		else
			yy += 1900;
	}
}


void RinexObsFile::printLoadingBar(int counter) {
	int refVal = round(double(m_numEpoch) / 100) * 2;

	if (counter == 0)
		cout << "Loading Rinex Observation ...  " << m_filename << endl;
	else if (counter % refVal == 0)
		cout << "*";

}


vector<int> RinexObsFile::getSatSysXIndicesAtEpoch(int atIndex, const char& satType) {
	vector<int> Temp;
	Temp.clear();
	for (unsigned int i = 0; i < obsList[atIndex].numSatInCurEpoch; i++) {
		
		if (obsList[atIndex].svIDInCurEpoch[i].at(0) == ' ')  // for the compatabilty for the rinex version prior to 2.x. 
			obsList[atIndex].svIDInCurEpoch[i].at(0) == 'G';  // Since ' ' is used for G, it is chaneged to 'G' 
		
		if (obsList[atIndex].svIDInCurEpoch[i].at(0) == satType)
			Temp.push_back(i);
	}
	return Temp;
}

vector<int> RinexObsFile::getSatSysXAtEpoch(int atIndex, const char& satType ) {
	vector<int> svIDList;
	vector<int> svIDIndices = getSatSysXIndicesAtEpoch(atIndex, satType);

	for (int i = 0; i < svIDIndices.size(); i++) {
		int id = atoi(obsList[atIndex].svIDInCurEpoch[svIDIndices[i]].substr(1, 2).c_str());
		svIDList.push_back(id);
	}
	return svIDList; // The List of the SV in the SatType at current epoch
}

int RinexObsFile::getIndexOfObservable(const string& observable) {

	for (int i = 0; i < m_numObservables; i++) {
		if (observable == m_observables[i])
			return i;
	}
	return -1; // if not available. 
}


/**
Returns the desired observation at the current epoch
@param atIndex, index of the obervation epoch
@param satType, Satellite Type desired. ( G = NAVSTAR-GPS,R GLONASS etc)
@param observationType, desired observation type. (Ex: L1, L2, C1, C2, P1, P2 etc)
#return the desired observation at the current epoch
*/
vector<double> RinexObsFile::getAnEpochOfObs(int atIndex, string observableType, const char& satType) {
	vector<double> obsListAtEpoch;

	int indexOfObservable = getIndexOfObservable(observableType);

	vector<int> svIndices = getSatSysXIndicesAtEpoch(atIndex, satType);

	for (int i = 0; i < svIndices.size(); i++) {
		double obsTemp = obsList[atIndex].data[svIndices[i]][indexOfObservable];
		obsListAtEpoch.push_back(obsTemp);
	}

	return obsListAtEpoch;
}

Car RinexObsFile::getAppCoord() { return m_appPos; }