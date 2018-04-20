% 
% ellipse_fitting.m
% 
% from paper 'Direct Least Squares Fitting of Ellipse, Andrew W Fitzgibbon, Maurizio Pilu'
% http://www.cnblogs.com/cv-pr/p/4625122.html
% https://stackoverflow.com/questions/21765015/unnormalization-of-ellipse-coefficients-after-direct-ellipse-fitting?answertab=votes#tab-top
%


function coeffs = ellipse_fitting(points)
if nargin == 0
  % Create an ellipse
  t = linspace(0, 2*pi);
  
  Rx = 300;
  Ry = 200;
  Cx = 250;
  Cy = 150;
  Rotation = pi + 0.1; % Radians
  
  NoiseLevel = 10.0; % Will add Gaussian noise of this std.dev. to points
  
  x = Rx * cos(t);
  y = Ry * sin(t);
  nx = x*cos(Rotation)-y*sin(Rotation) + Cx + randn(size(t))*NoiseLevel; 
  ny = x*sin(Rotation)+y*cos(Rotation) + Cy + randn(size(t))*NoiseLevel;
  
  % Clear figure
  clf
  % Draw it
  plot(nx,ny,'o');
  % Show the window
  figure(gcf)
  % Fit it
  params = ellipse_fitting([nx; ny]);
  % Note it may return (Rotation - pi/2) and swapped radii, this is fine.
  Given = round([Cx Cy Rx Ry Rotation*180/pi])
  Returned = round(params.*[1 1 1 1 180/pi])
  
  % Draw the returned ellipse
  t = linspace(0,pi*2);
  x = params(3) * cos(t);
  y = params(4) * sin(t);
  nx = x*cos(params(5))-y*sin(params(5)) + params(1); 
  ny = x*sin(params(5))+y*cos(params(5)) + params(2);
  hold on
  plot(nx,ny,'r-')
  
  return
end

X = points(1, :);
Y = points(2, :);
% normalize data
mx = mean(X);
my = mean(Y);
sx = (max(X) - min(X)) / 2;
sy = (max(Y) - min(Y)) / 2; 

x = (X-mx)/sx;
y = (Y-my)/sy;

% Force to column vectors
x = x(:);
y = y(:);

% Build design matrix
D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ];

% Build scatter matrix
S = D'*D;

% Build 6x6 constraint matrix


% Solve eigensystem
if 0
  C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;
  % Old way, numerically unstable if not implemented in matlab
  [gevec, geval] = eig(S, C);

  % Find the negative eigenvalue
  I = find(real(diag(geval)) < 1e-8 & ~isinf(diag(geval)));
  
  % Extract eigenvector corresponding to negative eigenvalue
  A = real(gevec(:,I));
else
  C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;
  % New way, numerically stabler in C [gevec, geval] = eig(S,C);
  
  % Break into blocks
  tmpA = S(1:3,1:3);
  tmpB = S(1:3,4:6);
  tmpC = S(4:6,4:6);
  tmpD = C(1:3,1:3);
  tmpE = inv(tmpC)*tmpB';
  [evec_x, eval_x] = eig(inv(tmpD) * (tmpA - tmpB*tmpE));
  
  % Find the positive (as det(tmpD) < 0) eigenvalue
  I = find(real(diag(eval_x)) < 1e-8 & ~isinf(diag(eval_x)));
  
  % Extract eigenvector corresponding to negative eigenvalue
  A = real(evec_x(:,I));
  
  % Recover the bottom half...
  evec_y = -tmpE * A;
  A = [A; evec_y];
end
  
% unnormalize
par = [A(1)*sy*sy,   ...
      A(2)*sx*sy,   ...
      A(3)*sx*sx,   ...
      -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
      -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
      A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
      - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
      + A(6)*sx*sx*sy*sy   ...
      ]';


% Convert to geometric radii, and centers

thetarad = 0.5*atan2(par(2),par(1) - par(3));
cost = cos(thetarad);
sint = sin(thetarad);
sin_squared = sint.*sint;
cos_squared = cost.*cost;
cos_sin = sint .* cost;

Ao = par(6);
Au =   par(4) .* cost + par(5) .* sint;
Av = - par(4) .* sint + par(5) .* cost;
Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;

% ROTATED = [Ao Au Av Auu Avv];

tuCentre = - Au./(2.*Auu);
tvCentre = - Av./(2.*Avv);
wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;

uCentre = tuCentre .* cost - tvCentre .* sint;
vCentre = tuCentre .* sint + tvCentre .* cost;

Ru = -wCentre./Auu;
Rv = -wCentre./Avv;
Ru = sqrt(abs(Ru)).*sign(Ru);
Rv = sqrt(abs(Rv)).*sign(Rv);

coeffs = [uCentre, vCentre, Ru, Rv, thetarad];
end