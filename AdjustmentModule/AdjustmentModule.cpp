#include "AdjustmentModule.h"

const double AdjustmentModule::C = constants::C;

/** 
    Returns the a priori values for the unknowns
	@param appPos is the approximate cartesian coordinate of station in m.
	@param dt is the receiver clock bias in sec.
	@return Matrix form of the a priori values for the unknown matrix.
*/

Matrix AdjustmentModule::setAprioriMatrix(Car& appPos, double dt) {
	Matrix T(4, 1);
	T.setValueAtElement(0, 0, appPos.x);
	T.setValueAtElement(1, 0, appPos.y);
	T.setValueAtElement(2, 0, appPos.z);
	T.setValueAtElement(3, 0, dt);

	return T;
}

/**
	Returns the a priori values for the unknowns
	@param apriori is a vector of 4 elements consisting of position and receiver clock bias.
	@return Matrix form of the a priori values for the unknown matrix.
*/
Matrix AdjustmentModule::setAprioriMatrix(vector<double> apriori) {
	Matrix T(4, 1);
	for (unsigned int i = 0; i < apriori.size(); i++) {
		T.setValueAtElement(i, 0, apriori.at(i));
	}
	   
	return T;
 }



/**
Returns the design matrix A
@param appPos, approximate position in Cartezian Coordinates
@param svPos, Satellite Positions in Cartezian Coordinates
*/
 Matrix AdjustmentModule::setDesignMatrix(Car & appPos, vector<Car>& svPos) {
	 const int NUM_OF_UNKNOWN = 4;
	 //const double C = 299792458;
	 Matrix Temp(svPos.size(), NUM_OF_UNKNOWN);
	 //vector<double> range;

	 for (int i = 0; i < svPos.size(); i++) {
		 double range = calcRange(appPos, svPos[i]) ;
		 double a1 = (appPos.x - svPos[i].x) / range;
		 double a2 = (appPos.y - svPos[i].y) / range;
		 double a3 = (appPos.z - svPos[i].z) / range;

		 Temp.setValueAtElement(i, 0, a1);
		 Temp.setValueAtElement(i, 1, a2);
		 Temp.setValueAtElement(i, 2, a3);
		 Temp.setValueAtElement(i, 3, C);
	 }
	 return Temp;
 }


 Matrix AdjustmentModule::setClosureMatrix(const double dt, const Car& appPos,
	                      const vector<Car>& svPos, const vector<double>& obs, vector<double>& satClk, vector<double>& tropDelay){
	 const double C = 299792458;
	 Matrix T(svPos.size(), 1);


	 for (int i = 0; i < svPos.size(); i++) {
		 double range = calcRange(appPos, svPos[i]);
		 double w1 = ( range + C * (dt-satClk[i]) + tropDelay[i] )- obs[i];
		 T.setValueAtElement(i, 0, w1);
	 }

	 return T;
 }

 Matrix AdjustmentModule::setClosureMatrix(const vector<double>& apriori, const vector<Car>& svPos, 
	                      const vector<double>& obs, vector<double>& satClk, vector<double>& tropDelay) {
	 
	 Matrix T(obs.size(), 1);

	 const double C = constants::C;

	 Car appPos;
	 appPos.x = apriori.at(0);
	 appPos.y = apriori.at(1);
	 appPos.z = apriori.at(2);
	 double dt = apriori.at(3);

	 for (int i = 0; i < obs.size(); i++) {
		 double range = calcRange(appPos, svPos[i]);
		 double w1 = (range + C * (dt - satClk[i]) + tropDelay[i]) - obs[i];
		 T.setValueAtElement(i, 0, w1);
	 }

	 return T;



 }

 /*
 static Matrix setClosureMatrix(Matrix& X0, const vector<CarCoord>& svPos, const vector<double>& obs, vector<double>& satClk){
	 const double C = 299792458;
	 Matrix T(svPos.size(), 1);
	 CarCoord appPos;
	 appPos.x = X0.getValueAtElement(0, 0);
	 appPos.y = X0.getValueAtElement(1, 0);
	 appPos.z = X0.getValueAtElement(2, 0);
	 double dt = X0.getValueAtElement(3, 0);

	 for (int i = 0; i < svPos.size(); i++) {
		 double range = calcRange(appPos, svPos[i]);
		 double w1 = (range + C * (dt - satClk[i])) - obs[i];
		 T.setValueAtElement(i, 0, w1);
	 }

	 return T;
 }
 */

 /**
 Returns the range between two points
 @param p1, Point 1
 @param p2, Point 2
 */
 double AdjustmentModule::calcRange(Car p1, Car p2) {
	 //double pr = sqrt(pow((appPos.x - satPos.x), 2) + pow((appPos.y - satPos.y), 2) + pow((appPos.z - satPos.z), 2));
	 return sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2) + pow((p1.z - p2.z), 2));
 }

 Matrix AdjustmentModule::estimateCorrections(const Matrix &A, const Matrix& W) {
	 Matrix N = (~A* A);
	 Matrix d = -1 * inv(N)*transpose(A)*W;//-1*inv(N)*~A*W;

	 return d;
 }

 Matrix AdjustmentModule::estimateUnknowns(const Matrix& X0, const Matrix &A, const Matrix& W) {
	 Matrix d = estimateCorrections(A, W);
	 return (X0 + d);
}
