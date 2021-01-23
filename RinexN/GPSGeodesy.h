#include<iostream>

using namespace std;

struct GregorianTime {
	int year, month, day, hour, minute;
	double sec;
};

struct GpsTime {
	int gpsweek, dow;  // dow = day of week
	double sow; // second of week
};


class GPSGeodesy {
public: 
	static double ArcTan(double x, double y);
	
	static double julday(GregorianTime t);
	static GpsTime  CalculateGpsTime(GregorianTime t);
	static void MjdToDate(long Mjd, long *Year, long *Month, long *Day);
	static long DateToMjd(long Year, long Month, long Day);
	long GpsToMjd(long GpsCycle, long GpsWeek, long GpsSeconds);
	static long MjdToJd(long Mjd);
	static long JdToMjd(long jd);


};





double GPSGeodesy::ArcTan(double x, double y)
{
	const double PI = std::atan(1.0) * 4;
	const double RHO = 180.0 / PI;
	double teta;
	teta = atan2(x, y);
	teta < 0 ? teta = teta + 2 * PI : teta = teta;
	return teta;
}


double GPSGeodesy::julday(GregorianTime t)
{
	double result;

	if (t.month <= 2)
	{
		t.year = t.year - 1;
		t.month = t.month + 12;
	}

	result = int(365.25 * t.year) + int(30.6001 * (t.month + 1)) + t.day +
		(double(t.hour) / 24) + (double(t.minute) / 60 / 24) + (t.sec / 3600 / 24) + 1720981.5;

	//cout << "******julday: "<< result << endl;

	return  result;
}


GpsTime  GPSGeodesy::CalculateGpsTime(GregorianTime t)
{
	GpsTime ans;

	const double SEC_PER_DAY = 86400.0;
	const GregorianTime JAN061980 = { 1980,1,6, 0, 0, 0.0 };
	double jan61980 = julday(JAN061980); // start of gpstime in terms of jd
	double mjd = julday(t) - 2400000.5;

	//    cout << " JAN061980 :" << jan61980 << endl;
	//    cout << " MJD       : " << mjd << endl;

	double jd = julday(t);
	int gps_week = (jd - jan61980) / 7;
	int day_of_week = int(jd - jan61980) % 7;
	double sec_of_week = ((jd - jan61980) - gps_week * 7) * SEC_PER_DAY;

	//    cout << " JD :" << jd << endl;
	//    cout << " GPS_WEEK " << gps_week << endl;
	//    cout << " DAY_OF_ WEEK " << day_of_week << endl;
	//    cout << sec_of_week << endl;

	ans.gpsweek = gps_week;
	ans.sow = sec_of_week;
	ans.dow = day_of_week;

	return ans;
}


/*
* Convert Modified Julian Day to calendar date.
* - Assumes Gregorian calendar.
* - Adapted from Fliegel/van Flandern ACM 11/#10 p 657 Oct 1968.
*/
void GPSGeodesy::MjdToDate(long Mjd, long *Year, long *Month, long *Day) {
	long J, C, Y, M;

	J = Mjd + 2400001 + 68569;
	C = 4 * J / 146097;
	J = J - (146097 * C + 3) / 4;
	Y = 4000 * (J + 1) / 1461001;
	J = J - 1461 * Y / 4 + 31;
	M = 80 * J / 2447;
	*Day = J - 2447 * M / 80;
	J = M / 11;
	*Month = M + 2 - (12 * J);
	*Year = 100 * (C - 49) + Y + J;
}

/*
* Return Modified Julian Day given calendar year,
* month (1-12), and day (1-31).
* - Valid for Gregorian dates from 17-Nov-1858.
* - Adapted from sci.astro FAQ.
*/

long GPSGeodesy::DateToMjd(long Year, long Month, long Day) {
	return
		367 * Year
		- 7 * (Year + (Month + 9) / 12) / 4
		- 3 * ((Year + (Month - 9) / 7) / 100 + 1) / 4
		+ 275 * Month / 9
		+ Day
		+ 1721028
		- 2400000;
}

/*
* Convert GPS Week and Seconds to Modified Julian Day.
* - Ignores UTC leap seconds.
*/

long GPSGeodesy::GpsToMjd(long GpsCycle, long GpsWeek, long GpsSeconds) {
	long GpsDays;

	GpsDays = ((GpsCycle * 1024) + GpsWeek) * 7 + (GpsSeconds / 86400);
	return DateToMjd(1980, 1, 6) + GpsDays;
}


long GPSGeodesy::MjdToJd(long Mjd) {
	return Mjd + 2400000.5;

}

long GPSGeodesy::JdToMjd(long jd) {
	return jd - 2400000.5;

}

