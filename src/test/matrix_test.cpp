#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
	mat A;
	A <<  4.9913e+01 << -5.0074e+01 << -3.2085e-14 << endr
  	  << -5.0074e+01 <<  5.0238e+01 <<  3.2196e-15 << endr
  	  << -3.2085e-14 <<  3.2196e-15 <<  1.0000e+02 << endr;

  	cout << "A = \n" << A << endl;

  	mat B = inv(A);
  	cout << "B = \n" << B << endl;

  	vec u;
	mat v;
	eig_sym(u, v, A);
	cout << "u = \n" << u << endl;
	cout << "v = \n" << v << endl;

	mat	S; 
   	S << 3.7243e+01 << -3.7361e+01 <<  3.7482e+01 <<  4.8177e-01 << -5.0223e-01 <<  4.9913e+01 << endr
  	  <<-3.7361e+01 <<  3.7482e+01 << -3.7604e+01 << -5.0223e-01 <<  5.2224e-01 << -5.0074e+01 << endr
   	  << 3.7482e+01 << -3.7604e+01 <<  3.7728e+01 <<  5.2224e-01 << -5.4180e-01 <<  5.0238e+01 << endr
   	  << 4.8177e-01 << -5.0223e-01 <<  5.2224e-01 <<  4.9913e+01 << -5.0074e+01 << -3.2085e-14 << endr
  	  <<-5.0223e-01 <<  5.2224e-01 << -5.4180e-01 << -5.0074e+01 <<  5.0238e+01 <<  3.2196e-15 << endr
      << 4.9913e+01 << -5.0074e+01 <<  5.0238e+01 << -3.2085e-14 <<  3.2196e-15 <<  1.0000e+02 << endr;

    mat tmpC = S(span(3, 5), span(3, 5));
    cout << "tmpC = \n" << tmpC << endl;
    cout << "InvC = \n" << inv(tmpC) << endl;
	return 0;
}