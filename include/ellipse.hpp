// ec_ellipse.hpp
// Andrew Wang, beandrewang@gmail.com

#ifndef ECARE_ELLIPSE_Hh_
#define ECARE_ELLIPSE_Hh_

#include <dlib/geometry.h>

namespace ecare
{
	class ellipse
	{
		struct feature
		{
			dlib::vector<double, 2> center;
			double					axis_ratio;
			double					rotation;
		};
	public:
		ellipse(dlib::matrix<double> points);
		feature read_features() { return f }
	private:
		bool ellipse_fitting(dlib::matrix<double> points);
	private:
		feature f;
	};
}

#endif