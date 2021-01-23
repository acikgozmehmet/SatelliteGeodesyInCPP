#include<iostream>
#include "..\BrdcEph\BrdcEph.h"
#include "..\RinexObsFile\RinexObsFile.h"
#include "..\GPSModule2\GPSModule.h"

using namespace std;

int menu();
void program(int choice);

string filename;
BrdcEph brdc;
RinexObsFile rinexo;

int main() {
	cout << " Single Point Positioning with GPS Data" << endl;

	int choice = 0;
	while (!(choice >= 1 && choice <= 5)) {
	   choice = menu();
	   program(choice);
	}

	cin.get();
	return 0;
}

int menu(){
    int choice;
    bool valid = false;
    while (!valid){
        cout << "************************************************************" << endl;
        cout << "* 1. Load Broadcast Ephemeris                              *" << endl;
        cout << "* 2. Load Precise Ephemeris                                *" << endl;
        cout << "* 3. Load Rinex Observation File                           *" << endl;
        cout << "* 4. Single Point Positioning (Epochwise)                  *" << endl;
		cout << "* 5. Exit                                                  *" << endl;
        cout << "************************************************************" << endl;
        cout << endl;
        cout << ">> ";
        cin >> choice;

        if ( (choice >= 1) && (choice <= 5) )
            return choice;
        else
            cout <<endl<< "Invalid choice. Please enter a number between 1 and 4" << endl << endl;
    }
}

void program(int choice) {

	switch(choice) {
	case 1:{
		   brdc.~BrdcEph(); // to clear the memory for the previous attempts
		   cout << "Enter the broadcast ephemeris filename: ";
		   cin >> filename;
		   brdc.setFilename(filename);
		   brdc.load();

		   cout << brdc.getFilename() << endl;

		   program(menu());
		   break;
	       }
	case 2:{
		   program(menu());
		   break;
	       }
	case 3:{
		    rinexo.~RinexObsFile(); // to clear the memory for the previous attempts
		    cout << "Enter the rinex observation filename: ";
		    cin >> filename;
		    rinexo.setFilename(filename);
		    rinexo.loadFromFile();

			cout << rinexo.getFilename() << endl;
			cout << rinexo.getNumEpoch() << endl;

		    program(menu());
		    break;
	        }
	case 4:{
		    // here is the estimation part.... update GPS module
			cout << brdc.getFilename() << endl;
			cout << rinexo.getFilename() << endl;
			cout << rinexo.getNumEpoch() << endl;

		    GPSModule::sppEpochbyEpoch(brdc, rinexo);

		    program(menu());
		    break;
	        }
	case 5:{
		    cout << "Bye..." << endl;
		    exit(0);
	        }
	default:
		cout << "Invalid input. Please enter the correct key" << endl;
	}
}