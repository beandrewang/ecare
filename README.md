# ecare

The main objective of this project is to implement a falling detector, using the ToF camera data.

## build

## build armadillo-code

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

Modeling the people skeloton with ellipse, take a reference from [paper](https://github.com/beandrewang/ecare/blob/master/src/scripts/ellipse-specific-fitting.pdf)