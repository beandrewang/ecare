#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
#if 0
  mat D;
  D <<  0.9601  <<  -0.8824 <<    0.811 <<  -0.9798 <<    0.9005  <<    1.0 << endr
    <<  0.9686  <<  -0.8837 <<    0.8063  <<  -0.9842 <<    0.8979  <<    1.0 << endr
    <<  0.9487  <<  -0.874  <<    0.8052  <<  -0.974  <<    0.8973  <<    1.0 << endr
    <<  0.9291  <<  -0.8498 <<    0.7773  <<  -0.9639 <<    0.8816  <<    1.0 << endr
    <<  0.9137  <<  -0.833  <<    0.7595  <<  -0.9559 <<    0.8715  <<    1.0 << endr
    <<  0.8649  <<  -0.7989 <<    0.7379  <<  -0.93 <<    0.859 <<    1.0 << endr
    <<  0.8481  <<  -0.7704 <<    0.6998  <<  -0.9209 <<    0.8365  <<    1.0 << endr
    <<  0.7856  <<  -0.7227 <<    0.6648  <<  -0.8863 <<    0.8154  <<    1.0 << endr
    <<  0.7413  <<  -0.6764 <<    0.6173  <<  -0.861  <<    0.7857  <<    1.0 << endr
    <<  0.6674  <<  -0.6207 <<    0.5774  <<  -0.8169 <<    0.7599  <<    1.0 << endr
    <<  0.6415  <<  -0.5792 <<    0.523 <<  -0.8009 <<    0.7232  <<    1.0 << endr
    <<  0.5641  <<  -0.518  <<    0.4757  <<  -0.7511 <<    0.6897  <<    1.0 << endr
    <<  0.4988  <<  -0.4582 <<    0.4209  <<  -0.7062 <<    0.6488  <<    1.0 << endr
    <<  0.4343  <<  -0.4006 <<    0.3696  <<  -0.659  <<    0.6079  <<    1.0 << endr
    <<  0.3867  <<  -0.3536 <<    0.3234  <<  -0.6218 <<    0.5687  <<    1.0 << endr
    <<  0.3201  <<  -0.2954 <<    0.2725  <<  -0.5658 <<    0.522 <<    1.0 << endr
    <<  0.2664  <<  -0.2436 <<    0.2227  <<  -0.5162 <<    0.4719  <<    1.0 << endr
    <<  0.2104  <<  -0.1937 <<    0.1783  <<  -0.4587 <<    0.4223  <<    1.0 << endr
    <<  0.1667  <<  -0.1519 <<    0.1383  <<  -0.4083 <<    0.3719  <<    1.0 << endr
    <<  0.1202  <<  -0.1111 <<    0.1027  <<  -0.3467 <<    0.3205  <<    1.0 << endr
    <<  0.0815  <<  -0.0747 <<    0.0685  <<  -0.2854 <<    0.2617  <<    1.0 << endr
    <<  0.0482  <<  -0.0446 <<    0.0412  <<  -0.2194 <<    0.2031  <<    1.0 << endr
    <<  0.0254  <<  -0.0237 <<    0.0221  <<  -0.1594 <<    0.1486  <<    1.0 << endr
    <<  0.0115  <<  -0.0103 <<    0.0092  <<  -0.1073 <<    0.0961  <<    1.0 << endr
    <<  0.0008  <<  -0.001  <<    0.0013  <<  -0.0281 <<    0.0355  <<    1.0 << endr
    <<  0.0007  <<  -0.0005 <<    0.0004  <<  0.027 <<    -0.0188 <<    1.0 << endr
    <<  0.0078  <<  -0.0071 <<    0.0065  <<  0.0883  <<    -0.0805 <<    1.0 << endr
    <<  0.0234  <<  -0.0219 <<    0.0205  <<  0.1531  <<    -0.1432 <<    1.0 << endr
    <<  0.0441  <<  -0.042  <<    0.04  <<  0.2099  <<    -0.1999 <<    1.0 << endr
    <<  0.0777  <<  -0.0707 <<    0.0643  <<  0.2787  <<    -0.2537 <<    1.0 << endr
    <<  0.1125  <<  -0.1055 <<    0.0989  <<  0.3354  <<    -0.3145 <<    1.0 << endr
    <<  0.1603  <<  -0.1439 <<    0.1291  <<  0.4004  <<    -0.3593 <<    1.0 << endr
    <<  0.2025  <<  -0.1833 <<    0.1658  <<  0.45  <<    -0.4072 <<    1.0 << endr
    <<  0.2606  <<  -0.2399 <<    0.2208  <<  0.5105  <<    -0.4699 <<    1.0 << endr
    <<  0.3075  <<  -0.2849 <<    0.2639  <<  0.5545  <<    -0.5138 <<    1.0 << endr
    <<  0.3693  <<  -0.3427 <<    0.3181  <<  0.6077  <<    -0.564  <<    1.0 << endr
    <<  0.4367  <<  -0.4006 <<    0.3675  <<  0.6608  <<    -0.6062 <<    1.0 << endr
    <<  0.4945  <<  -0.4529 <<    0.4148  <<  0.7032  <<    -0.6441 <<    1.0 << endr
    <<  0.5563  <<  -0.5136 <<    0.4742  <<  0.7459  <<    -0.6886 <<    1.0 << endr
    <<  0.6296  <<  -0.5739 <<    0.5231  <<  0.7935  <<    -0.7232 <<    1.0 << endr
    <<  0.692 <<  -0.6292 <<    0.5722  <<  0.8319  <<    -0.7564 <<    1.0 << endr
    <<  0.7417  <<  -0.6802 <<    0.6239  <<  0.8612  <<    -0.7898 <<    1.0 << endr
    <<  0.7973  <<  -0.7323 <<    0.6725  <<  0.8929  <<    -0.8201 <<    1.0 << endr
    <<  0.8413  <<  -0.775  <<    0.714 <<  0.9172  <<    -0.845  <<    1.0 << endr
    <<  0.8768  <<  -0.8051 <<    0.7393  <<  0.9364  <<    -0.8598 <<    1.0 << endr
    <<  0.9307  <<  -0.8537 <<    0.783 <<  0.9647  <<    -0.8849 <<    1.0 << endr
    <<  0.9536  <<  -0.8793 <<    0.8108  <<  0.9765  <<    -0.9005 <<    1.0 << endr
    <<  0.9941  <<  -0.9065 <<    0.8265  <<  0.9971  <<    -0.9091 <<    1.0 << endr
    <<  0.9998  <<  -0.9151 <<    0.8375  <<  0.9999  <<    -0.9151 <<    1.0 << endr
    <<  1.0247  <<  -0.9367 <<    0.8563  <<  1.0123  <<    -0.9254 <<    1.0 << endr
    <<  1.0068  <<  -0.9224 <<    0.845 <<  1.0034  <<    -0.9192 <<    1.0 << endr
    <<  0.9936  <<  -0.9189 <<    0.8497  <<  0.9968  <<    -0.9218 <<    1.0 << endr
    <<  0.9817  <<  -0.8981 <<    0.8217  <<  0.9908  <<    -0.9065 <<    1.0 << endr
    <<  0.9653  <<  -0.887  <<    0.815 <<  0.9825  <<    -0.9028 <<    1.0 << endr
    <<  0.931 <<  -0.856  <<    0.787 <<  0.9649  <<    -0.8871 <<    1.0 << endr
    <<  0.8921  <<  -0.8221 <<    0.7576  <<  0.9445  <<    -0.8704 <<    1.0 << endr
    <<  0.8399  <<  -0.7765 <<    0.7179  <<  0.9165  <<    -0.8473 <<    1.0 << endr
    <<  0.8055  <<  -0.7356 <<    0.6718  <<  0.8975  <<    -0.8196 <<    1.0 << endr
    <<  0.7539  <<  -0.6912 <<    0.6337  <<  0.8683  <<    -0.7961 <<    1.0 << endr
    <<  0.6765  <<  -0.624  <<    0.5756  <<  0.8225  <<    -0.7587 <<    1.0 << endr
    <<  0.6292  <<  -0.5763 <<    0.5278  <<  0.7932  <<    -0.7265 <<    1.0 << endr
    <<  0.5574  <<  -0.5134 <<    0.473 <<  0.7466  <<    -0.6877 <<    1.0 << endr
    <<  0.4988  <<  -0.4627 <<    0.4291  <<  0.7063  <<    -0.6551 <<    1.0 << endr
    <<  0.4422  <<  -0.4027 <<    0.3668  <<  0.665 <<    -0.6056 <<    1.0 << endr
    <<  0.3793  <<  -0.3454 <<    0.3145  <<  0.6159  <<    -0.5608 <<    1.0 << endr
    <<  0.3114  <<  -0.2883 <<    0.267 <<  0.5581  <<    -0.5167 <<    1.0 << endr
    <<  0.2556  <<  -0.2366 <<    0.219 <<  0.5056  <<    -0.468  <<    1.0 << endr
    <<  0.2048  <<  -0.184  <<    0.1652  <<  0.4526  <<    -0.4065 <<    1.0 << endr
    <<  0.1513  <<  -0.1399 <<    0.1294  <<  0.389 <<    -0.3597 <<    1.0 << endr
    <<  0.116 <<  -0.1049 <<    0.0949  <<  0.3405  <<    -0.308  <<    1.0 << endr
    <<  0.0766  <<  -0.0703 <<    0.0645  <<  0.2768  <<    -0.2539 <<    1.0 << endr
    <<  0.0456  <<  -0.0413 <<    0.0373  <<  0.2136  <<    -0.1932 <<    1.0 << endr
    <<  0.0233  <<  -0.0221 <<    0.0209  <<  0.1528  <<    -0.1447 <<    1.0 << endr
    <<  0.0085  <<  -0.0072 <<    0.0061  <<  0.0923  <<    -0.0781 <<    1.0 << endr
    <<  0.0006  <<  -0.0005 <<    0.0005  <<  0.0248  <<    -0.0218 <<    1.0 << endr
    <<  0.0012  <<  -0.0011 <<    0.0011  <<  -0.0344 <<    0.0324  <<    1.0 << endr
    <<  0.0105  <<  -0.0101 <<    0.0097  <<  -0.1024 <<    0.0984  <<    1.0 << endr
    <<  0.0249  <<  -0.0234 <<    0.022 <<  -0.1578 <<    0.1484  <<    1.0 << endr
    <<  0.05  <<  -0.0459 <<    0.0421  <<  -0.2237 <<    0.2052  <<    1.0 << endr
    <<  0.0801  <<  -0.0723 <<    0.0651  <<  -0.2831 <<    0.2552  <<    1.0 << endr
    <<  0.1195  <<  -0.11 <<    0.1012  <<  -0.3457 <<    0.3182  <<    1.0 << endr
    <<  0.1681  <<  -0.1504 <<    0.1346  <<  -0.41 <<    0.3668  <<    1.0 << endr
    <<  0.2039  <<  -0.1905 <<    0.178 <<  -0.4516 <<    0.4219  <<    1.0 << endr
    <<  0.269 <<  -0.2439 <<    0.2211  <<  -0.5187 <<    0.4702  <<    1.0 << endr
    <<  0.3177  <<  -0.2917 <<    0.2678  <<  -0.5637 <<    0.5175  <<    1.0 << endr
    <<  0.389 <<  -0.3564 <<    0.3266  <<  -0.6237 <<    0.5715  <<    1.0 << endr
    <<  0.4364  <<  -0.4059 <<    0.3774  <<  -0.6606 <<    0.6143  <<    1.0 << endr
    <<  0.5083  <<  -0.4667 <<    0.4286  <<  -0.713  <<    0.6546  <<    1.0 << endr
    <<  0.5675  <<  -0.5182 <<    0.4732  <<  -0.7533 <<    0.6879  <<    1.0 << endr
    <<  0.6361  <<  -0.5757 <<    0.5211  <<  -0.7976 <<    0.7219  <<    1.0 << endr
    <<  0.6897  <<  -0.6334 <<    0.5816  <<  -0.8305 <<    0.7627  <<    1.0 << endr
    <<  0.7428  <<  -0.6835 <<    0.6291  <<  -0.8618 <<    0.7931  <<    1.0 << endr
    <<  0.7779  <<  -0.7218 <<    0.6697  <<  -0.882  <<    0.8184  <<    1.0 << endr
    <<  0.8511  <<  -0.7701 <<    0.6968  <<  -0.9225 <<    0.8348  <<    1.0 << endr
    <<  0.8745  <<  -0.8039 <<    0.7389  <<  -0.9351 <<    0.8596  <<    1.0 << endr
    <<  0.8945  <<  -0.8239 <<    0.7589  <<  -0.9458 <<    0.8712  <<    1.0 << endr
    <<  0.9216  <<  -0.8528 <<    0.789 <<  -0.96 <<    0.8883  <<    1.0 << endr
    <<  0.9621  <<  -0.8793 <<    0.8036  <<  -0.9809 <<    0.8965  <<    1.0 << endr
    <<  0.9542  <<  -0.8829 <<    0.8169  <<  -0.9768 <<    0.9038  <<    1.0 << endr
    <<  0.9756  <<  -0.8879 <<    0.808 <<  -0.9877 <<    0.8989  <<    1.0 << endr;

mat S = D.t() * D;
#else
	mat	S; 
 	S << 3.7243e+01 << -3.7361e+01 <<  3.7482e+01 <<  4.8177e-01 << -5.0223e-01 <<  4.9913e+01 << endr
	  <<-3.7361e+01 <<  3.7482e+01 << -3.7604e+01 << -5.0223e-01 <<  5.2224e-01 << -5.0074e+01 << endr
 	  << 3.7482e+01 << -3.7604e+01 <<  3.7728e+01 <<  5.2224e-01 << -5.4180e-01 <<  5.0238e+01 << endr
 	  << 4.8177e-01 << -5.0223e-01 <<  5.2224e-01 <<  4.9913e+01 << -5.0074e+01 << -3.2085e-14 << endr
	  <<-5.0223e-01 <<  5.2224e-01 << -5.4180e-01 << -5.0074e+01 <<  5.0238e+01 <<  3.2196e-15 << endr
    << 4.9913e+01 << -5.0074e+01 <<  5.0238e+01 << -3.2085e-14 <<  3.2196e-15 <<  1.0000e+02 << endr;
#endif
  cout << "S = \n" << S << endl;
  mat C = S(span(3, 5), span(3, 5));
  cout << "C = \n" << C << endl;
  cout << "cond(C)" << cond(C) << endl;
  cout << "rcond(C)" << rcond(C) << endl;
  cout << "InvC = \n" << inv(C) << endl;
	return 0;
}