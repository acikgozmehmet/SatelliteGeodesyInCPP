#include<iostream>
#include<iomanip>
#include "Ellipsoid.h"

using namespace std;

int main() {
	{
		Ellipsoid elpsd;
		elpsd.printParameters(cout);
	}
	cout << endl;
	{
		Ellipsoid ell(6378137.0000, 6356752.3142);
		ell.printParameters(cout);
	}
	cout << endl;
	{
		Ellipsoid ell(6378137.0000, 6356752.3142);
		Ellipsoid ell2(ell);
		ell2.printParameters(cout);
	}

	cin.ignore();
	cin.get();

	return 0;
}
