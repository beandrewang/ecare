#include "ellipse.hpp"

using namespace ecare;

ellipse::ellipse(const dlib::matrix<double> points)
{
	ellipse_fitting(points);
}

bool ellipse::ellipse_fitting(const dlib::matrix<double> points)
{
	dlib::matrix<double> X = dlib::colm(points, 0);
	dlib::matrix<double> Y = dlib::colm(points, 1);
}