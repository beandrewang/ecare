#include "ellipse.hpp"

ellipse::ellipse(dlib::matrix<double> points)
{
	f = ellipse_fitting(points);
}

bool ellipse::ellipse_fitting(dlib::matrix<double> points)
{
	
}