A = [
	 4.9913e+01  -5.0074e+01  -3.2085e-14;
    -5.0074e+01   5.0238e+01   3.2196e-15;
    -3.2085e-14   3.2196e-15   1.0000e+02; 
    ]

B = inv(A)

[v, u] = eig(A)