# ecare

The main objective of this project is to implement a falling detector, using the ToF camera data.

## build

```
mkdir build
cd build
cmake ../src
make -j8
```

After building, you can get your excutable application `ecare` in `build` directory.

## Human body modeling

Modeling the people skeloton with ellipse, take a reference from [paper](https://github.com/beandrewang/ecare/blob/master/src/scripts/ellipse-specific-fitting.pdf)