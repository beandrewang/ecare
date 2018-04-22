// ec_ellipse.hpp
// Andrew Wang, beandrewang@gmail.com

#ifndef ECARE_ELLIPSE_Hh_
#define ECARE_ELLIPSE_Hh_

#include <dlib/matrix.h>

namespace ecare
{
	using namespace std;
	using namespace dlib;

	class ellipse
	{
		struct feature
		{
			double					cx;
			double					cy;
			double					a;
			double 					b;
			double					theta;
		};
		typedef matrix<double> mat;
		typedef matrix<double> vec;
	public:
		ellipse();
		ellipse(const mat &points);
		feature read_features() { return f; }
	private:
		bool ellipse_fitting();
		bool generate_points();
		bool design_matrix(const vec &x, const vec &y, mat &S);
		bool solve_equation(const mat &S, vec &A);
		bool normalize(vec &x, vec &y);
		bool unnormalize(const vec &A, vec &par);
		bool computer_geometry(const vec &par);
	private:
		feature f;
		mat scatter;
		double mx;
		double my;
		double sx;
		double sy;
	};
}

#endif