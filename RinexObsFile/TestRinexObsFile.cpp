#include<iostream>
#include "RinexObsFile.h"
using namespace std;

int main() {

	RinexObsFile Test;
	Test.setFilename("wtzr0200.17o");
	Test.loadFromFile();


	cout << "Number of epochs in the file is " << Test.getNumEpoch() << endl;

	vector<int> svList;
	vector<double> obs;



	for (unsigned int j = 0; j < Test.getNumEpoch(); j++) {
		//svList.clear();
		//obs.clear();

		

		
		svList = Test.getSatSysXAtEpoch(j);
        obs = Test.getAnEpochOfObs(j, "C1", 'G');
		cout << "Epoch : " << j << "  "; Test.getEpoch(j).print(); cout << endl;

		for (unsigned int i = 0; i < svList.size(); i++)
			cout << svList[i] << "    " << fixed << setprecision(10) << obs[i] << endl;
	    
		cout << endl;

	}


	cin.ignore();
	cin.get();


	return 0;
}