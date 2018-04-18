// ec_ellipse.hpp
// Andrew Wang, beandrewang@gmail.com

#ifndef ECARE_ELLIPSE_Hh_
#define ECARE_ELLIPSE_Hh_

#include <dlib/matrix.h>

namespace ecare
{
	class ellipse
	{
		struct feature
		{
			dlib::matrix<double, 2, 1> center;
			double					a;
			double 					b;
			double					theta;
		};

		template <typename T> int sgn(T val) {
		    return (T(0) < val) - (val < T(0));
		}
	public:
		ellipse();
		ellipse(const dlib::matrix<double> points);
		feature read_features() { return f; }
	private:
		bool ellipse_fitting(const dlib::matrix<double> points);
		bool generate_points(dlib::matrix<double> &points);
	private:
		feature f;
	};
}

#endif