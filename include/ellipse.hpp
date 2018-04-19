// ec_ellipse.hpp
// Andrew Wang, beandrewang@gmail.com

#ifndef ECARE_ELLIPSE_Hh_
#define ECARE_ELLIPSE_Hh_

#include <armadillo>

using namespace std;
using namespace arma;

namespace ecare
{
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

	public:
		ellipse();
		ellipse(const mat &points);
		feature read_features() { return f; }
	private:
		template <typename T> int sgn(T val) {
		    return (T(0) < val) - (val < T(0));
		}
		bool ellipse_fitting(const mat &points);
		bool generate_points(mat &points);
	private:
		feature f;
	};
}

#endif