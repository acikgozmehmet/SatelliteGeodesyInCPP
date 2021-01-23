#include<iostream>
#include<iomanip>
#include "TimeSystems.h"

using namespace std;


int main() {

	
		cout << "Test-1" << endl;
		int yy = 2017;
		int mm = 3;
		int dd = 13;
		int hh = 0;
		int minute = 0;
		double sec = 0;
		cout << fixed;
		cout << "Date to convert is : " << yy << " " << mm << " " << dd << " " << hh << " " << minute << " " << sec << endl;
		double jd = TimeSystems::DateToJd(yy, mm, dd, hh, minute, sec);
		cout <<"Julian Day : " << jd << endl;
		double mjd = TimeSystems::DateToMjd(yy, mm, dd, hh, minute, sec);
		cout << "Modified Julian Day :" << mjd << endl;
		cout << "Julian Day to confirm : " << TimeSystems::MjdToJd(mjd) << endl;

		int y2, m2, d2, h2, mn2;
		double sc2;
		TimeSystems::JdToDate(jd, y2, m2, d2, h2, mn2, sc2);
		cout << "Date to confirm !: " << y2 << " " << m2 << " " << d2 << " " << h2 << " " << mn2 << " " << sc2 << endl;


	
	//int mjd = TimeSystems::DateToMjd(yy, mm, dd, hh, minute, sec);
	//double jd = TimeSystems::MjdToJd(mjd);
	//cout << fixed << jd << endl;

	cout << endl;
	int yy2, mm2, dd2, hh2, min2;
	double sec2;
	TimeSystems::JdToDate(jd, yy2, mm2, dd2, hh2, min2, sec2);
	cout << yy2 << " " << mm2  << "  " << dd2 << endl;
	cout << hh2 << " " << min2 << "  " << sec2 << endl;
	cout << endl;

	int doy = TimeSystems::DateToDoy(yy2, mm2, dd2, hh2, min2, sec2);
	cout << doy << endl;
	int yy3, mm3, dd3, hh3, min3;
	double sec3;
	TimeSystems::DoyToDate(yy2, doy, mm3, dd3, hh2, min2, sec2);

	cout << yy2 << " " << mm3 << "  " << dd3 << endl;
	cout << hh2 << " " << min2 << "  " << sec2 << endl;

	cout << endl;

	double sow;
	int gweek, dow;
	// 17  1 19 23 59 44.0
	TimeSystems::DateToGPSTime(yy, mm, dd,hh,minute,sec, sow, gweek, dow);
	//TimeSystems::DateToGPSTime(yy2, mm2, dd2, hh2, min2, sec2, sow, gweek, dow);
	cout << "GPS Week : " << gweek <<" Second of Week : "<< sow <<" Day of Week : "<< dow << endl;
	cout << TimeSystems::GPSWeekToGPSCycle(gweek) << endl;



	cout << "Remaining :" << (int(sow) % 86400) / 7200 << endl;

	int yy4, mm4, dd4, hh4, min4;
	double sec4;


	TimeSystems::GPSTimeToDate(gweek, sow, yy4, mm4, dd4, hh4, min4, sec4);
	cout << yy4 << " " << mm4 << "  " << dd4 << endl;
	cout << hh4 << " " << min4 << "  " << sec4 << endl;






	cin.ignore();
	cin.get();

	return 0;
}