#ifndef ADJUSTMENT_MODULE_H
#define ADJUSTMENT_MODULE_H
#include<vector>
#include "..\GPSModule2\GPSModule.h"
#include "..\Matrix\Matrix.h"


class AdjustmentModule {
public:
	static Matrix setAprioriMatrix(Car & appPos, double dt);
	static Matrix setAprioriMatrix(vector<double> apriori);
	static Matrix setDesignMatrix(Car & appPos, vector<Car>& svPos);
	static Matrix setClosureMatrix(const double dt, const Car& appPos, const vector<Car>& svPos, const vector<double>& obs, vector<double>& satClk, vector<double>& tropDelay);
	static Matrix setClosureMatrix(const vector<double>& apriori, const vector<Car>& svPos, const vector<double>& obs, vector<double>& satClk, vector<double>& tropDelay);

	//static Matrix setClosureMatrix(const Matrix& X0, const vector<CarCoord>& svPos, const vector<double>& obs, vector<double>& satClk);

	static double calcRange(Car p1, Car p2);
	static Matrix estimateCorrections(const Matrix &A, const Matrix& W);
	static Matrix estimateUnknowns(const Matrix& X0, const Matrix &A, const Matrix& W);

private:
	static const double C;
};

#endif