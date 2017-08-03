#ifdef _DEBUG
//Debugモードの場合
#pragma comment(lib,"C:\\opencv\\build\\x86\\vc12\\lib\\opencv_world300d.lib")            // opencv_core
#else
//Releaseモードの場合
#pragma comment(lib,"C:\\opencv\\build\\x86\\vc12\\lib\\opencv_world300.lib") 
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
/*数値計算ライブラリ*/
#define _USE_MATH_DEFINES
#include <math.h>
/*可変長リストライブラリ*/
#include <vector>
/*OpenCVライブラリ*/
#include "opencv2/opencv.hpp"
#include "opencv2/highgui/highgui_c.h"

using namespace std;
using namespace cv;

class TripleClothoidCurve{
	/////////////////////////////////////////////////
	//
	// Private Variable
	//
	Point2d point_start;
	Point2d vector_start;
	double radius_start;
	double curvature_start;

	Point2d point_end;
	Point2d vector_end;
	double radius_end;
	double curvature_end;

	Point2d temp_point;

	vector<Point2d> curve;

	double h;
	double r;

	double max_line_length;
	double max_a12;

	double S1, S2;
	double step_draw;
	double step_integral_sin;
	double step_integral_cos;
	double step_search_a;
	double step_search_h;
	double accuracy_width;
	double accuracy_length;

	double t11, t12, t13, t21, t22, t23;
	double a10, a11, a12, a20, a21, a22, a30, a31, a32;

	double theta0;

	double phi0;
	double phi1;

	double k0;
	double k1;

	double S;
	double phi;
	//
	// Private Variable
	//
	/////////////////////////////////////////////////

	/////////////////////////////////////////////////
	//
	// Private Function
	//
	void CalcParamA(){
		a22 = (a12*(t11*t23 - t13*t21) + h*(k0*(t23 - t13) + k1*t13) - t23*(phi1 - phi0)) / (t13*t22 - t12*t23);
		a32 = (a12*(t11*t22 - t12*t21) + h*(k0*(t22 - t12) + k1*t12) - t22*(phi1 - phi0)) / (t12*t23 - t13*t22);

		a10 = phi0 - theta0;
		a11 = k0 * h;

		a20 = a10 + (a11 * S1) + (a12 * S1 * S1);
		a21 = a11 + (2 * a12 * S1);

		a30 = a20 + a21 * (S2 - S1) + a22 * (S2 - S1) * (S2 - S1);
		a31 = a21 + 2 * a22 * (S2 - S1);
	}
	double SinPhi1(double x){
		return sin(a10 + a11 * x + a12 * x * x);
	}
	double SinPhi2(double x){
		return sin(a20 + a21 * (x - S1) + a22 * (x - S1) * (x - S1));
	}
	double SinPhi3(double x){
		return sin(a30 + a31 * (x - S2) + a32 * (x - S2) * (x - S2));
	}
	double CosPhi1(double x){
		return cos(a10 + a11 * x + a12 * x * x);
	}
	double CosPhi2(double x){
		return cos(a20 + a21 * (x - S1) + a22 * (x - S1) * (x - S1));
	}
	double CosPhi3(double x){
		return cos(a30 + a31 * (x - S2) + a32 * (x - S2) * (x - S2));
	}
	double IntegralSin(){
		double ans = 0;
		for (double S = 0; S < S1; S += step_integral_sin){
			ans += (SinPhi1(S) + SinPhi1(S + step_integral_sin)) * step_integral_sin / 2;
		}
		for (double S = S1; S < S2; S += step_integral_sin){
			ans += (SinPhi2(S) + SinPhi2(S + step_integral_sin)) * step_integral_sin / 2;
		}
		for (double S = S2; S < 1; S += step_integral_sin){
			ans += (SinPhi3(S) + SinPhi3(S + step_integral_sin)) * step_integral_sin / 2;
		}
		return ans;
	}
	double IntegralCos(){
		double ans = 0;
		for (double S = 0; S < S1; S += step_integral_cos){
			ans += (CosPhi1(S) + CosPhi1(S + step_integral_cos)) * step_integral_cos / 2;
		}
		for (double S = S1; S < S2; S += step_integral_cos){
			ans += (CosPhi2(S) + CosPhi2(S + step_integral_cos)) * step_integral_cos / 2;
		}
		for (double S = S2; S < 1; S += step_integral_cos){
			ans += (CosPhi3(S) + CosPhi3(S + step_integral_cos)) * step_integral_cos / 2;
		}
		return ans;
	}
	double GetSolustionSin(){ return IntegralSin(); }
	double GetSolustionCos(){ return h * IntegralCos() - r; }
	double ErrorSumOfSquares(){ return GetSolustionSin() * GetSolustionSin() + GetSolustionCos() * GetSolustionCos(); }
	double DeltaESOSa12(double da12, double dh, double dx){
		a12 = da12 - dx;
		h = dh;
		CalcParamA();
		double fh1 = ErrorSumOfSquares();

		a12 = da12 + dx;
		h = dh;
		CalcParamA();
		double fh2 = ErrorSumOfSquares();

		return fh2 - fh1;
	}
	double DeltaESOSh(double da12, double dh, double dx){
		a12 = da12;
		h = dh - dx;
		CalcParamA();
		double fh1 = ErrorSumOfSquares();

		a12 = da12;
		h = dh + dx;
		CalcParamA();
		double fh2 = ErrorSumOfSquares();

		return fh2 - fh1;
	}
	bool ZeroCross(double in1, double in2){
		if (in1 * in2 <= 0){
			return true;
		}
		return false;
	}
	int SearchNestOnece(double nest_initial_a, double nest_initial_h, double dx, double num, double lr, vector<double> &list_a, vector<double> &list_h){
		double temp_a12;
		double temp_h;
		temp_a12 += -lr * DeltaESOSa12(nest_initial_a, nest_initial_h, dx);
		temp_h += -lr * DeltaESOSh(nest_initial_a, nest_initial_h, dx);

		for (int i = 0; i < num; ++i){
			a12 += -lr * DeltaESOSa12(temp_a12, h, dx);
			h += -lr * DeltaESOSh(temp_a12, h, dx);
		}
		list_a.push_back(a12);
		list_h.push_back(h);

		return 0;
	}
	int SearchParamAHNext(){
		vector<double> list_a;
		vector<double> list_h;

		////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////

		SearchNestOnece(0, r,
			0.0001, 10000, 0.5,
			list_a, list_h);

		if (list_a.empty()){ return -1; }
		a12 = list_a[0];
		h = list_h[0];
		CalcParamA();
		return 0;
	}
	//
	// Private Function
	//
	/////////////////////////////////////////////////
public:

	/////////////////////////////////////////////////
	//
	// Set
	//
	void SetMaxLineLength(double length){
		max_line_length = length;
	}
	void SetMaxa12(double a12){
		max_a12 = a12;
	}
	void SetPoint(Point2d start, Point2d end){
		point_start = start;
		point_end = end;

		r = sqrt((point_end.x - point_start.x) * (point_end.x - point_start.x) + (point_end.y - point_start.y) * (point_end.y - point_start.y));
		theta0 = atan2(point_end.y - point_start.y, point_end.x - point_start.x);
	}
	void SetAngle(double start, double end){
		radius_start = start;
		radius_end = end;

		vector_start = Point2d(cos(radius_start), sin(radius_start));
		vector_end = Point2d(cos(radius_end), sin(radius_end));

		phi0 = radius_start;
		phi1 = radius_end;

	}
	void SetCurvature(double start, double end){
		curvature_start = start;
		curvature_end = end;

		k0 = curvature_start;
		k1 = curvature_end;
	}
	void SetAccuracy(double width, double length, double integral_sin, double integral_cos, double search_a, double search_h){
		accuracy_width = width;
		accuracy_length = length;
		step_integral_sin = integral_sin;
		step_integral_cos = integral_cos;
		step_search_a = search_a;
		step_search_h = search_h;
	}
	void SetAccuracy(double width, double length, double integral_sin, double integral_cos){
		accuracy_width = width;
		accuracy_length = length;
		step_integral_sin = integral_sin;
		step_integral_cos = integral_cos;
		step_search_a = max_a12 * 2 / 100;
		step_search_h = ((max_line_length * r) - r) / 100;
	}
	void SetAccuracy(double width, double length){
		accuracy_width = width;
		accuracy_length = length;
		step_integral_sin = 0.001;
		step_integral_cos = 0.001;
		step_search_a = max_a12 * 2 / 100;
		step_search_h = ((max_line_length * r) - r) / 100;
	}
	void SetDrawStep(double step){
		step_draw = step;
	}
	void SetS1S2(double in1, double in2){
		S1 = in1;
		S2 = in2;

		t11 = -S1 * S1 + 2 * S1;
		t12 = -S2 * S2 + 2 * S2 + S1 * S1 - 2 * S1;
		t13 = S2 * S2 - 2 * S2 + 1;
		t21 = 2 * S1;
		t22 = 2 * S2 - 2 * S1;
		t23 = 2 - 2 * S2;

		a12 = 0;
		h = 0;
	}
	//
	// Set
	//
	/////////////////////////////////////////////////////

	/////////////////////////////////////////////////////
	//
	// Create
	//
	int CreateCurve(vector<Point2d> &list_curve = vector<Point2d>(0)){
		S = 0;
		phi = 0;
		double step = step_draw;
		temp_point = point_start;

		curve.clear();

		if (SearchParamAHNext()){ return -1; }
		curve.push_back(temp_point);
		for (double S = step; S < S1; S += step){
			phi = theta0 + a10 + a11 * S + a12 * S * S;

			temp_point += h * Point2d(step * cos(phi), step * sin(phi));
			curve.push_back(temp_point);
		}
		for (double S = S1; S < S2; S += step){
			phi = theta0 + a20 + a21 * (S - S1) + a22 * (S - S1) * (S - S1);

			temp_point += h * Point2d(step * cos(phi), step * sin(phi));
			curve.push_back(temp_point);
		}
		for (double S = S2; S < 1; S += step){
			phi = theta0 + a30 + a31 * (S - S2) + a32 * (S - S2) * (S - S2);

			temp_point += h * Point2d(step * cos(phi), step * sin(phi));
			curve.push_back(temp_point);
		}
		phi = theta0 + a30 + a31 * (1 - S2) + a32 * (1 - S2) * (1 - S2);
		temp_point += h * Point2d(step * cos(phi), step * sin(phi));
		curve.push_back(temp_point);

		list_curve = curve;

		return 0;
	}
	void CorrectEndPoint(){
		Point2d difference = point_end - curve.back();

		for (size_t i = 0; i < curve.size(); i++){
			curve[i] = curve[i] + (difference  * (i / double(curve.size() - 1)));
		}
	}
	//
	// Create
	//
	/////////////////////////////////////////////////////

	/////////////////////////////////////////////////////
	//
	// OutPut
	//
	void PrintParam(){
		cout << "S1  " << S1 << ", " << endl;
		cout << "S2  " << S2 << ", " << endl;
		cout << "a10 " << a10 << ", " << phi0 - theta0 << endl;
		cout << "a11 " << a11 << ", " << k0 * h << endl;
		cout << "a12 " << a12 << ", " << endl;
		cout << "a20 " << a20 << ", " << a10 + a11 * S1 + a12 * S1 * S1 << endl;
		cout << "a21 " << a21 << ", " << a11 + 2 * a12 * S1 << endl;
		cout << "a22 " << a22 << ", " << (a12 * (t11 * t23 - t13 * t21) + h * (k0 * (t23 - t13) + k1 * t13) - (phi1 - phi0) * t23) / (t13 * t22 - t12 * t23) << endl;
		cout << "a30 " << a30 << ", " << a20 + a21 * (S2 - S1) + a22 * (S2 - S1) * (S2 - S1) << endl;
		cout << "a31 " << a31 << ", " << a21 + 2 * a22 * (S2 - S1) << endl;
		cout << "a32 " << a32 << ", " << (a12 * (t11 * t22 - t12 * t21) + h * (k0 * (t22 - t12) + k1 * t12) - (phi1 - phi0) * t22) / (t12 * t23 - t13 * t22) << endl;
		cout << "phi0 " << phi0 << ", " << a10 << endl;
		cout << "phi1 " << phi1 << ", " << a30 + a31 * (1 - S2) + a32 * (1 - S2) * (1 - S2) << ", " << phi << endl;
		cout << "theta0 " << theta0 << endl;
		cout << "h   " << h << ", " << endl;
		cout << "r   " << r << ", " << endl;
		cout << "k0  " << k0 << ", " << a11 / h << endl;
		cout << "k1  " << k1 << ", " << (a31 + 2 * a32 * (1 - S2)) / h << endl;
		cout << "k0.1 " << (a11 + 2 * a12 * S1) / h << endl;
		cout << "k0.2 " << (a21 + 2 * a22 * (S2 - S1)) / h << endl;

		cout << "intesin " << IntegralSin() << endl;
		cout << "intecos " << IntegralCos() << endl;
		cout << "absh " << abs(r - h * IntegralCos()) << endl;
		cout << "absa " << abs(IntegralSin()) << endl;
		cout << "arival Point" << curve.back() << endl;
	}
	void GetParam(
		double &S1_out, double &S2_out,
		double &a10_out, double &a11_out, double &a12_out,
		double &a20_out, double &a21_out, double &a22_out,
		double &a30_out, double &a31_out, double &a32_out,
		double &phi0_out, double &phi1_out,
		double &theta0_out,
		double &h_out,
		double &r_out,
		double &k0_out, double &k1_out){

		S2_out = S2;
		S1_out = S1;
		a10_out = phi0 - theta0;
		a11_out = k0 * h;
		a12_out = a12;
		a20_out = a10 + a11 * S1 + a12 * S1 * S1;
		a21_out = a11 + 2 * a12 * S1;
		a22_out = (a12 * (t11 * t23 - t13 * t21) + h * (k0 * (t23 - t13) + k1 * t13) - (phi1 - phi0) * t23) / (t13 * t22 - t12 * t23);
		a30_out = a20 + a21 * (S2 - S1) + a22 * (S2 - S1) * (S2 - S1);
		a31_out = a21 + 2 * a22 * (S2 - S1);
		a32_out = (a12 * (t11 * t22 - t12 * t21) + h * (k0 * (t22 - t12) + k1 * t12) - (phi1 - phi0) * t22) / (t12 * t23 - t13 * t22);
		phi0_out = a10;
		phi1_out = a30 + a31 * (1 - S2) + a32 * (1 - S2) * (1 - S2);
		theta0_out = theta0;
		h_out = h;
		r_out = r;
		k0 = a11 / h;
		k1 = (a31 + 2 * a32 * (1 - S2)) / h;

		return;
	}
	void PrintCurve(){
		cout << curve << endl;
	}
	void GetCurve(vector<Point2d> &out_curve){
		out_curve = curve;
	}
	//
	// OutPut
	//
	/////////////////////////////////////////////////////


};
int main(){
	TripleClothoidCurve clothoid;
	clothoid.SetMaxLineLength(3);
	clothoid.SetMaxa12(100);
	clothoid.SetDrawStep(0.00001);
	clothoid.SetS1S2(0.25, 0.75);
	clothoid.SetPoint(Point2d(13.8179, 61.248), Point2d(120, 20));
	clothoid.SetAccuracy(0.1, 1, 0.01, 0.01);
	clothoid.SetAngle(-M_PI / 2, 0);
	clothoid.SetCurvature(0, 0);

	vector<Point2d> curve;
	clothoid.CreateCurve(curve);
	clothoid.GetCurve(curve);


	//画像出力
	Mat image(Size(1024, 1024), CV_8U, Scalar::all(0));
	line(image, Point(0, 511), Point(1023, 511), Scalar(127), 1, CV_AA);
	line(image, Point(511, 0), Point(511, 1023), Scalar(127), 1, CV_AA);
	for (vector<Point2d>::iterator itr = curve.begin(); itr < curve.end(); itr++){
		image.at<unsigned char>(Point2d(itr->x, -itr->y) + Point2d(511, 511)) = 255;
	}
	imshow("ClothoidCurve", image);


	//ファイル出力
	ofstream ofs("C:\\Users\\0133752\\Desktop\\Clothoid.csv");
	for (int i = 0; i < curve.size(); ++i){
		ofs << setprecision(16) << curve[i].x << "," << setprecision(16) << curve[i].y << endl;
	}

	//	curve3.CorrectEndPoint();
	//	clothoid.PrintCurve();
	clothoid.PrintParam();

	waitKey(0);
	return 0;
}