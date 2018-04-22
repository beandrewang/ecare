# ecare

The main objective of this project is to implement a falling detector, using the ToF camera data.

## build

### build armadillo-code

```
cd module/armadillo-code
mkdir build
cmake ../ -DCMAKE_INSTALL_PREFIX:PATH=../dev
make & make install
```

You will get directory `dev`, which consist of include and libs ,in armadillo-code,

```
mkdir build
cd build
cmake ../src
make -j8
```

After building, you can get your excutable application `ecare` in `build` directory.

## Human body modeling

Modeling the people skeloton with ellipse, take a reference from [Direct Least Square Fitting of Ellipse, Andrew W. Fitzigibbon](https://github.com/beandrewang/ecare/blob/master/ref/ellipse-specific-fitting.pdf)

We set the feature of an ellipse like bellow:

```
struct feature
{
	double					cx;
	double					cy;
	double					a;
	double 					b;
	double					theta;
};
```

then we can get such conclusion:

* if length of a, b is exchanged, a1 = b2 and a2 = b1, then theta1 - theta2 equal to 90 or 270
* if length of a, b is not changed, a1 = a2 and b1 = b2, then theta1 - theta2 equal to 0 or 180
