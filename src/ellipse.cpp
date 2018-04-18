#include <iostream>
#include <cmath>
#include "ellipse.hpp"

using namespace ecare;

ellipse::ellipse()
{
	mat points;
	generate_points(points);
	ellipse_fitting(points);
}

ellipse::ellipse(const mat &points)
{
	ellipse_fitting(points);
}

bool ellipse::generate_points(mat &points)
{
	// Create an ellipse
	vec t = linspace<vec>(0, 2*M_PI, 100);
	double Rx = 300.0;
	double Ry = 200.0;
	double Cx = 250.0;
	double Cy = 150.0;
	double Rotation = M_PI/2 + 0.1; //Radians

	double NoiseLevel = 1.0; // Will add Gaussian noise of this std.dev. to points

	vec x = Rx * cos(t);
	vec y = Ry * cos(t);

	vec nx = x*cos(Rotation)-y*sin(Rotation) + Cx + randn<vec>(t.size())*NoiseLevel;
	vec ny = x*sin(Rotation)+y*cos(Rotation) + Cy + randn<vec>(t.size())*NoiseLevel;

	points = join_cols(nx, ny);
	points = trans(points);

	return true;
}

bool ellipse::ellipse_fitting(const mat &points)
{
	//std::cout << "points = \n" << points << std::endl;
	mat X = points.col(0);
	mat Y = points.col(1);
	//std::cout << "X = \n" << X << "\n" << "Y = \n" << Y << std::endl;

	// normalize data
	double mx = mean(X);
	double my = mean(Y);
	double sx = (max(X) - min(X)) / 2;
	double sy = (max(Y) - min(Y)) / 2;
	//std::cout << "mx = \n" << mx << "\n" << "my = \n" << my << std::endl;

	mat x = (X - mx) / sx;
	mat y = (Y - my) / sy;
	//std::cout << "x = \n" << x << "\n" << "y = \n" << y << std::endl;

	// build design matrix
	mat D;
	D = join_rows(x % x, x % y);
	D = join_rows(D, y % y);
	D = join_rows(D, x);
	D = join_rows(D, y);
	D = join_rows(D, ones<vec>(size(x));
	//std::cout << "D = \n" << D << std::endl;

	// Build scatter matrix
	mat S = trans(D) * D;

	std::cout << "S = \n" << S << std::endl;
	std::cout << "pinv(S) = \n" << pinv(S) << std::endl;
	std::cout << "I = \n" << S * pinv(S) << std::endl;

	// Build 6x6 constraint matrix
	mat C = zeros<mat>(6, 6);

	#if 0
	C(0, 2) = -2; C(1, 1) = 1; C(2, 0) = -2;
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
	d = dlib::diag(d);
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
	dlib::matrix<double> A = dlib::colm(v, i);

	// Recover the bottom half...
	dlib::matrix<double> yvec = -tmpE * A;
	A = dlib::join_cols(A, yvec);
	#else
	C(0, 2) = -2; C(1, 1) = 1; C(2, 0) = -2;

	cx_vec geval;
	cx_mat gevec;
	eig_pair(geval, gevec, S, C);

	// Todo, here should be modified later
	I = find(real(diag(geval)) < 1e-8));
	
	// Extract eigenvector corresponding to negative eigenvalue
	vec A = v.colm(I);
	#endif
	std::cout << "C = \n" << C << std::endl;
	std::cout << "d = \n" << d << std::endl;
	std::cout << "v = \n" << v << std::endl;
	std::cout << "A = \n" << A << std::endl;

	#if 1
	// unnormalize
	vec par;
	par(0) = A(0)*sy*sy;
	par(1) = A(1)*sx*sy;
	par(2) = A(2)*sx*sx;
	par(3) = -2*A(0)*sy*sy*mx - A(1)*sx*sy*my + A(3)*sx*sy*sy;
	par(4) = -A(1)*sx*sy*mx - 2*A(2)*sx*sx*my + A(4)*sx*sx*sy;
	par(5) = A(0)*sy*sy*mx*mx + A(1)*sx*sy*mx*my + A(2)*sx*sx*my*my - A(3)*sx*sy*sy*mx - A(4)*sx*sx*sy*my + A(5)*sx*sx*sy*sy;

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

	std::cout << f.cx << f.cy << "\n" << f.a << "\n" << f.b << "\n" << f.theta << std::endl;

	return true;
}