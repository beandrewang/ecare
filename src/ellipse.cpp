#include <iostream>
#include <cmath>
#include "ellipse.hpp"

using namespace ecare;

ellipse::ellipse()
{
	dlib::matrix<double> points;
	generate_points(points);
	ellipse_fitting(points);
}

ellipse::ellipse(const dlib::matrix<double> points)
{
	ellipse_fitting(points);
}

bool ellipse::generate_points(dlib::matrix<double> &points)
{
	// Create an ellipse
	dlib::matrix<double> t = dlib::linspace(0, 2*M_PI, 100);
	double Rx = 300.0;
	double Ry = 200.0;
	double Cx = 250.0;
	double Cy = 150.0;
	double Rotation = M_PI/2 + 0.1; //Radians

	double NoiseLevel = 1.0; // Will add Gaussian noise of this std.dev. to points

	dlib::matrix<double> x = Rx * dlib::cos(t);
	dlib::matrix<double> y = Ry * dlib::cos(t);

	dlib::matrix<double> nx = x*cos(Rotation)-y*sin(Rotation) + Cx + dlib::randm(t.size(), 1)*NoiseLevel;
	dlib::matrix<double> ny = x*sin(Rotation)+y*cos(Rotation) + Cy + dlib::randm(t.size(), 1)*NoiseLevel;

	points = dlib::join_cols(nx, ny);
	points = dlib::trans(points);

	return true;
}

bool ellipse::ellipse_fitting(const dlib::matrix<double> points)
{
	//std::cout << "points = \n" << points << std::endl;
	dlib::matrix<double> X = dlib::colm(points, 0);
	dlib::matrix<double> Y = dlib::colm(points, 1);
	//std::cout << "X = \n" << X << "\n" << "Y = \n" << Y << std::endl;

	// normalize data
	double mx = dlib::mean(X);
	double my = dlib::mean(Y);
	double sx = (dlib::max(X) - dlib::min(X)) / 2;
	double sy = (dlib::max(Y) - dlib::min(Y)) / 2;
	//std::cout << "mx = \n" << mx << "\n" << "my = \n" << my << std::endl;

	dlib::matrix<double> x = (X - mx) / sx;
	dlib::matrix<double> y = (Y - my) / sy;
	//std::cout << "x = \n" << x << "\n" << "y = \n" << y << std::endl;

	// build design matrix
	dlib::matrix<double> D;
	D = dlib::join_rows(dlib::pointwise_multiply(x, x), dlib::pointwise_multiply(x, y));
	D = dlib::join_rows(D, dlib::pointwise_multiply(y, y));
	D = dlib::join_rows(D, x);
	D = dlib::join_rows(D, y);
	D = dlib::join_rows(D, dlib::ones_matrix<double>(x.size(), 1));
	//std::cout << "D = \n" << D << std::endl;

	// Build scatter matrix
	dlib::matrix<double> S = dlib::trans(D) * D;
	std::cout << "S = \n" << S << std::endl;
	std::cout << "pinv(S) = \n" << dlib::pinv(S) << std::endl;
	std::cout << "I = \n" << S * dlib::pinv(S) << std::endl;

	// Build 6x6 constraint matrix
	dlib::matrix<double> C = dlib::zeros_matrix<double>(6, 6);

	C(0, 2) = 2; C(1, 1) = -1; C(2, 0) = 2;
	std::cout << "C = \n" << C << std::endl;
	#if 0
	// genralized eig
	dlib::matrix<double, 3, 3> tmpA = dlib::subm(S, dlib::range(0, 2), dlib::range(0, 2));
	dlib::matrix<double, 3, 3> tmpB = dlib::subm(S, dlib::range(0, 2), dlib::range(3, 5));
	dlib::matrix<double, 3, 3> tmpC = dlib::subm(S, dlib::range(3, 5), dlib::range(3, 5));
	dlib::matrix<double, 3, 3> tmpD = dlib::subm(C, dlib::range(0, 2), dlib::range(0, 2));
	dlib::matrix<double, 3, 3> tmpE = dlib::pinv(tmpC) * tmpB;

	dlib::eigenvalue_decomposition<dlib::matrix<double>> ed(dlib::inv(tmpD) * (tmpA - tmpB * tmpE)); 
	
	// Find the positive (as det(tmpD) < 0) eigenvalue
	dlib::matrix<double> d = ed.get_pseudo_d();
	dlib::matrix<double> v = ed.get_pseudo_v();
	//std::cout << d << std::endl;
	//std::cout << v << std::endl;
	d = dlib::real(dlib::diag(d));
	std::cout << d << std::endl;

	int i;
	for(i = 0; i < d.nr(); i++)
	{
		if (d(i, 0) < 1e-8 && std::isinf(d(i, 0)))
		{
			break;
		}
	}

	// Extract eigenvector corresponding to negative eigenvalue
	dlib::matrix<double> A = dlib::real(v(:, i))

	// Recover the bottom half...
	dlib::matrix<double> y = -tmpE * A;
	A = 
	#else
	dlib::eigenvalue_decomposition<dlib::matrix<double>> ed(dlib::pinv(S) * C);
	dlib::matrix<double> d = ed.get_pseudo_d();
	dlib::matrix<double> v = ed.get_pseudo_v();
	int i;
	for(i = 0; i < d.nr(); i++)
	{
		if (d(i, 0) > 1e-8 && std::isinf(d(i, 0)))
		{
			break;
		}
	}

	dlib::matrix<double> A = dlib::colm(v, i);
	#endif
	std::cout << "d = \n" << d << std::endl;
	std::cout << "v = \n" << v << std::endl;
	std::cout << "A = \n" << A << std::endl;

	#if 1
	// unnormalize
	dlib::matrix<double, 6, 1> par;
	par(0, 0) = A(0, 0)*sy*sy;
	par(1, 0) = A(1, 0)*sx*sy;
	par(2, 0) = A(2, 0)*sx*sx;
	par(3, 0) = -2*A(0, 0)*sy*sy*mx - A(1, 0)*sx*sy*my + A(3, 0)*sx*sy*sy;
	par(4, 0) = -A(1, 0)*sx*sy*mx - 2*A(2, 0)*sx*sx*my + A(4, 0)*sx*sx*sy;
	par(5, 0) = A(0, 0)*sy*sy*mx*mx + A(1, 0)*sx*sy*mx*my + A(2, 0)*sx*sx*my*my - A(3, 0)*sx*sy*sy*mx - A(4, 0)*sx*sx*sy*my + A(5, 0)*sx*sx*sy*sy;

	//std::cout << "par = \n" << par << std::endl;
	#endif

	// Convert to geometric radii, and centers
	double thetarad = 0.5 * atan2(par(2), par(1) - par(3));
	double cost = cos(thetarad);
	double sint = sin(thetarad);
	double sin_squared = sint * sint;
	double cos_squared = cost *cost;
	double cos_sin = sint * cost;

	double Ao = par(6);
	double Au =   par(4) * cost + par(5) * sint;
	double Av = - par(4) * sint + par(5) * cost;
	double Auu = par(1) * cos_squared + par(3) * sin_squared + par(2) * cos_sin;
	double Avv = par(1) * sin_squared + par(3) * cos_squared - par(2) * cos_sin;

	// % ROTATED = [Ao Au Av Auu Avv]
	double tuCentre = - Au/(2*Auu);
	double tvCentre = - Av/(2*Avv);
	double wCentre = Ao - Auu*tuCentre*tuCentre - Avv*tvCentre*tvCentre;

	double uCentre = tuCentre * cost - tvCentre * sint;
	double vCentre = tuCentre * sint + tvCentre * cost;

	double Ru = -wCentre/Auu;
	double Rv = -wCentre/Avv;

	Ru = sqrt(fabs(Ru))*sgn(Ru);
	Rv = sqrt(fabs(Rv))*sgn(Rv);

	f.center = uCentre, 
			   vCentre;
	f.a      = Ru;
	f.b 	 = Rv;
	f.theta  = thetarad;

	std::cout << f.center << "\n" << f.a << "\n" << f.b << "\n" << f.theta << std::endl;

	return true;
}