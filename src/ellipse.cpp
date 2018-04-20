#include <iostream>
#include "ellipse.hpp"

using namespace ecare;

ellipse::ellipse()
{
	generate_points();
	ellipse_fitting();
}

ellipse::ellipse(const mat points)
{
	scatter = points;
	ellipse_fitting();
}

bool ellipse::generate_points()
{
	// Create an ellipse
	vec t = linspace<vec>(0, 2 * datum::pi);
	int Rx = 300;
	int Ry = 200;
	int Cx = 250;
	int Cy = 150;
	double Rotation = datum::pi + 0.1; //Radians
	double NoiseLevel = 0.0; // Will add Gaussian noise of this std.dev. to points

	vec x = Rx * cos(t);
	vec y = Ry * sin(t);

	vec nx = x*cos(Rotation)-y*sin(Rotation) + Cx + randn<vec>(t.size())*NoiseLevel;
	vec ny = x*sin(Rotation)+y*cos(Rotation) + Cy + randn<vec>(t.size())*NoiseLevel;

	scatter = join_rows(nx, ny);
	
	cout << "Given = \n" << Cx << "\t" << Cy << "\t" << Rx << "\t" << Ry << "\t" << Rotation * 180 / datum::pi << endl;
	return true;
}

bool ellipse::normalize(vec &x, vec &y)
{
	x = scatter.col(0);
	y = scatter.col(1);

	mx = mean(x);
	my = mean(y);
	sx = (max(x) - min(x)) / 2;
	sy = (max(x) - min(y)) / 2;

	x = (x - mx) / sx;
	y = (y - my) / sy;

	return true;
}

bool ellipse::construct_features(const vec &x, const vec &y, mat &S)
{
	mat D;
	D = join_rows(x % x, x % y);
	D = join_rows(D, y % y);
	D = join_rows(D, x);
	D = join_rows(D, y);
	D = join_rows(D, ones<vec>(size(x)));
	
	S = D.t() * D;

	return true;
}

bool ellipse::solve_equation(const mat &S, vec &A)
{
	mat C = zeros<mat>(6, 6);
	C(0, 2) = 2.0; C(1, 1) = -1.0; C(2, 0) = 2.0;
	//std::cout << "S = \n" << S << std::endl;

	mat tmpA = S(span(0, 2), span(0, 2));
	mat tmpB = S(span(0, 2), span(3, 5));
	mat tmpC = S(span(3, 5), span(3, 5));
	//cout << "tmpC = \n" << tmpC << endl;
	//cout << "invC = \n" << inv(tmpC) << endl;
	mat tmpD = C(span(0, 2), span(0, 2)); 
	//mat tmpE = inv(tmpC) * trans(tmpB);
	mat tmpE = solve(tmpC, trans(tmpB));

	cx_vec geval;
	cx_mat gevec;
	//eig_gen(geval, gevec, inv(tmpD) * (tmpA - tmpB * tmpE));
	eig_gen(geval, gevec,solve(tmpD, (tmpA - tmpB * tmpE)));

	uvec I = find(real(geval.diag()) < 1e-8 && geval.diag() != datum::inf);
	A = real(gevec.col(I(0)));
	vec evec_y = -tmpE * A;
	A = join_cols(A, evec_y);

	return true;
}

bool ellipse::unnormalize(const vec &A, vec &par)
{
	par << A(0) * sy * sy << endr
	    << A(1) * sx * sy << endr
	    << A(2) * sx * sx << endr
	    << -2 * A(0) * sy * sy * mx - A(1) * sx * sy * my + A(3) * sx * sy * sy << endr
	    << -A(1) * sx * sy * mx - 2 * A(2) * sx * sx * my + A(4) * sx * sx * sy << endr
	    << A(0) * sy * sy * mx * mx + A(1) * sx * sy * mx * my + A(2) * sx * sx * my * my - \
	       A(3) * sx * sy * sy * mx - A(4) * sx * sx * sy * my + A(5) * sx * sx * sy * sy << endr;

	return true;
}

bool ellipse::computer_features(const vec &par)
{
	double thetarad = 0.5 * atan2(par(1), par(0) - par(2));
	double cost = cos(thetarad);
	double sint = sin(thetarad);
	double sin_squared = sint * sint;
	double cos_squared = cost *cost;
	double cos_sin = sint * cost;

	double Ao = par(5);
	double Au =   par(3) * cost + par(4) * sint;
	double Av = - par(3) * sint + par(4) * cost;
	double Auu = par(0) * cos_squared + par(2) * sin_squared + par(1) * cos_sin;
	double Avv = par(0) * sin_squared + par(2) * cos_squared - par(1) * cos_sin;

	double tuCentre = - Au/(2*Auu);
	double tvCentre = - Av/(2*Avv);
	double wCentre = Ao - Auu*tuCentre*tuCentre - Avv*tvCentre*tvCentre;

	double uCentre = tuCentre * cost - tvCentre * sint;
	double vCentre = tuCentre * sint + tvCentre * cost;

	double Ru = -wCentre/Auu;
	double Rv = -wCentre/Avv;

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

	f.cx = uCentre;
	f.cy = vCentre;
	f.a      = Ru;
	f.b 	 = Rv;
	f.theta  = thetarad;

	cout << "Returned = \n" << f.cx << "\t" << f.cy << "\t" << f.a << "\t" << f.b << "\t" << f.theta * 180 / datum::pi << endl;

	return true;
}

bool ellipse::ellipse_fitting()
{
	vec x, y;
	mat D;
	vec A;
	vec par;
	normalize(x, y);
	construct_features(x, y, D);
	solve_equation(D, A);
	unnormalize(A, par);
	computer_features(par);
	return true;
}
