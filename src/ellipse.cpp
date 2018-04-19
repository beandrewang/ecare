#include <iostream>
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

	points = join_rows(nx, ny);
	cout << points << endl;
	return true;
}

bool ellipse::ellipse_fitting(const mat &points)
{
	//std::cout << "points = \n" << points << std::endl;
	vec X = points.col(0);
	vec Y = points.col(1);
	std::cout << "X = \n" << X << "\n" << "Y = \n" << Y << std::endl;

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
	D = join_rows(D, ones<vec>(size(x)));
	std::cout << "D = \n" << D << std::endl;

	// Build 6x6 constraint matrix
	mat C = zeros<mat>(6, 6);
	// Build scatter matrix
	mat S;
	S = D.t() * D;
	#if 0
	mat	S; 
   	S << 3.7243e+01 << -3.7361e+01 <<  3.7482e+01 <<  4.8177e-01 << -5.0223e-01 <<  4.9913e+01 << endr
  	  <<-3.7361e+01 <<  3.7482e+01 << -3.7604e+01 << -5.0223e-01 <<  5.2224e-01 << -5.0074e+01 << endr
   	  << 3.7482e+01 << -3.7604e+01 <<  3.7728e+01 <<  5.2224e-01 << -5.4180e-01 <<  5.0238e+01 << endr
   	  << 4.8177e-01 << -5.0223e-01 <<  5.2224e-01 <<  4.9913e+01 << -5.0074e+01 << -3.2085e-14 << endr
  	  <<-5.0223e-01 <<  5.2224e-01 << -5.4180e-01 << -5.0074e+01 <<  5.0238e+01 <<  3.2196e-15 << endr
      << 4.9913e+01 << -5.0074e+01 <<  5.0238e+01 << -3.2085e-14 <<  3.2196e-15 <<  1.0000e+02 << endr;
    #endif
	std::cout << "S = \n" << S << std::endl;

	#if 1
	C(0, 2) = -2; C(1, 1) = 1; C(2, 0) = -2;
	// genralized eig
	mat tmpA = S(span(0, 2), span(0, 2));
	cout << "tmpA = \n" << tmpA << endl;
	mat tmpB = S(span(0, 2), span(3, 5));
	cout << "tmpB = \n" << tmpB << endl;
	mat tmpC = S(span(3, 5), span(3, 5));
	std::cout << "tmpC = \n" << tmpC << std::endl;
	std::cout << "invC = \n" << inv(tmpC) << endl;
	mat tmpD = C(span(0, 2), span(0, 2));
	std::cout << "tmpD = \n" << tmpD << std::endl;
	mat tmpE = inv(tmpC) * trans(tmpB);
	
	std::cout << "tmpE = \n" << tmpE << std::endl;

	cx_vec geval;
	cx_mat gevec;
	mat W = inv(tmpD) * (tmpA - tmpB * tmpE);
	eig_gen(geval, gevec, W);
	cout << "W = \n" << W << endl;
	std::cout << "C = \n" << C << std::endl;
	std::cout << "geval = \n" << geval << std::endl;
	std::cout << "gevec = \n" << gevec << std::endl;
	

	uvec I = find(real(geval.diag()) < 1e-8 && geval.diag() != datum::inf);
	cout << "I = \n" << I << endl;
	// Extract eigenvector corresponding to negative eigenvalue
	vec A = real(gevec.col(I(0)));

	cout << "A = \n" << A << endl;
	// % Recover the bottom half...
	cout << "tmpE = \n" << tmpE << endl;
	vec evec_y = -tmpE * A;
	cout << "evec_y = \n" << evec_y << endl;
  	A = join_cols(A, evec_y);
  	std::cout << "A = \n" << A << std::endl;
	#else
	C(0, 2) = -2; C(1, 1) = 1; C(2, 0) = -2;

	cx_vec geval;
	cx_mat gevec;
	eig_pair(geval, gevec, S, C);
	// Todo, here should be modified later
	uvec I = find(real(geval.diag()) < 1e-8 && geval.diag() != datum::inf);
	cout << "I = \n" << I << endl;
	// Extract eigenvector corresponding to negative eigenvalue
	vec A = real(gevec.col(I(0)));
	#endif

	// unnormalize
	vec par(6);
	par(0) = A(0)*sy*sy;
	par(1) = A(1)*sx*sy;
	par(2) = A(2)*sx*sx;
	par(3) = -2*A(0)*sy*sy*mx - A(1)*sx*sy*my + A(3)*sx*sy*sy;
	par(4) = -A(1)*sx*sy*mx - 2*A(2)*sx*sx*my + A(4)*sx*sx*sy;
	par(5) = A(0)*sy*sy*mx*mx + A(1)*sx*sy*mx*my + A(2)*sx*sx*my*my - A(3)*sx*sy*sy*mx - A(4)*sx*sx*sy*my + A(5)*sx*sx*sy*sy;

	std::cout << "par = \n" << par << std::endl;

	// Convert to geometric radii, and centers
	double thetarad = 0.5 * atan2(par(1), par(0) - par(2));
	double cost = cos(thetarad);
	double sint = sin(thetarad);
	double sin_squared = sint * sint;
	double cos_squared = cost *cost;
	double cos_sin = sint * cost;

	std::cout << "thetarad = \n" << thetarad << std::endl;
	std::cout << "cost = \n" << cost << std::endl;
	std::cout << "sint = \n" << sint << std::endl;
	std::cout << "sin_squared = \n" << sin_squared << std::endl;
	std::cout << "cos_squared = \n" << cos_squared << std::endl;
	std::cout << "cos_sin = \n" << cos_sin << std::endl;

	double Ao = par(5);
	double Au =   par(3) * cost + par(4) * sint;
	double Av = - par(3) * sint + par(4) * cost;
	double Auu = par(0) * cos_squared + par(2) * sin_squared + par(1) * cos_sin;
	double Avv = par(0) * sin_squared + par(2) * cos_squared - par(1) * cos_sin;

	std::cout << "Ao = \n" << Ao << std::endl;
	std::cout << "Au = \n" << Au << std::endl;
	std::cout << "Av = \n" << Av << std::endl;
	std::cout << "Auu = \n" << Auu << std::endl;
	std::cout << "Avv = \n" << Avv << std::endl;

	// % ROTATED = [Ao Au Av Auu Avv]
	double tuCentre = - Au/(2*Auu);
	double tvCentre = - Av/(2*Avv);
	double wCentre = Ao - Auu*tuCentre*tuCentre - Avv*tvCentre*tvCentre;

	double uCentre = tuCentre * cost - tvCentre * sint;
	double vCentre = tuCentre * sint + tvCentre * cost;

	double Ru = -wCentre/Auu;
	double Rv = -wCentre/Avv;

	std::cout << "tuCentre = \n" << tuCentre << std::endl;
	std::cout << "tvCentre = \n" << tvCentre << std::endl;
	std::cout << "wCentre = \n" << wCentre << std::endl;
	std::cout << "uCentre = \n" << uCentre << std::endl;
	std::cout << "vCentre = \n" << vCentre << std::endl;
	std::cout << "Ru = \n" << Ru << std::endl;
	std::cout << "Rv = \n" << Rv << std::endl;

	int Su = 0;
	if(Ru > 0)
	{
		Su = 1;
	}
	else if(Ru < 0)
	{
		Su = -1;
	}

	int Sv = 0;
	if(Rv > 0)
	{
		Sv = 1;
	}
	else if(Rv < 0)
	{
		Sv = -1;
	}

	Ru = sqrt(abs(Ru)) * Su;
	Rv = sqrt(abs(Rv)) * Sv;

	std::cout << "Ru = \n" << Ru << std::endl;
	std::cout << "Rv = \n" << Rv << std::endl;

	f.cx = uCentre;
	f.cy = vCentre;
	f.a      = Ru;
	f.b 	 = Rv;
	f.theta  = thetarad;

	std::cout << f.cx << "\n" << f.cy << "\n" << f.a << "\n" << f.b << "\n" << f.theta << std::endl;

	return true;
}
