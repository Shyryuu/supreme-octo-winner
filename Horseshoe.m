clear

% Wing definition

b = 12;
b1 = 5.4;

S = 3.3*(1.6+1.2)+5.4*1.6;

width = 0.1;
y = 0:width:12;

nwpoints = length(y);
nwel = nwpoints-1;

ycp = y(1:end-1)+width/2;

ct = 1.2;
cr = 1.6;

c1 = ct+(cr-ct)/(3.3-0)*y(y<3.3);
c2 = cr*ones(1,sum(y>=3.3&y<=8.7));
c3 = cr-(cr-ct)/(3.3)*(y(y>8.7)-8.7);

c = [c1,c2,c3];

xu = c/2;
xd = -c/2;

xcor = c/4;
ycor = y;
cslice = (c(2:end)+c(1:end-1))/2;

xcp = -cslice/4;

twist = (-6*2*pi/360)/6*abs(ycp-6);

alpha0 = DVM(2412);

eps = -alpha0+twist;

n = [sin(eps)',cos(eps)'];

XCP = xcp'*ones(1,nwel);

XCOR1 = (xcor(1:end-1)'*ones(1,nwel))';
XCOR2 = (xcor(2:end)'*ones(1,nwel))';

X1 = XCP-XCOR1;
X2 = XCP-XCOR2;

YCP = ycp'*ones(1,nwel);

YCOR1 = (ycor(1:end-1)'*ones(1,nwel))';
YCOR2 = (ycor(2:end)'*ones(1,nwel))';

Y1 = YCP-YCOR1;
Y2 = YCP-YCOR2;

R1 = (X1.^2+Y1.^2).^0.5;
R2 = (X2.^2+Y2.^2).^0.5;

Resc = X1.*X2+Y1.*Y2;
Rvect = -X1.*Y2+X2.*Y1;

VAB = 1/(4*pi)*(R1+R2)./(R1.*R2.*(R1.*R2+Resc)).*Rvect;

Uesc1 = X1./R1;
Uvect1 = -Y1./R1;

Uesc2 = X2./R2;
Uvect2 = -Y2./R2;

VA = 1/(4*pi)*(1-Uesc1)./Uvect1;
VB = 1/(4*pi)*(1-Uesc2)./Uvect2;

V = VA+VAB-VB;

nz = n(:,2);

NZ = nz*ones(1,nwel);
A = V.*NZ;

arange = [-5 10]*2*pi/360;
apoints = 20;


alphapoints = linspace(arange(1),arange(2),apoints);

%for i=1:length(alphapoints)

  alpha = alphapoints(1);
  alpha = 5;
  Uinf = 1*[-cos(alpha);sin(alpha)];
  b = -n*Uinf;
  
  vort = A\b;
  
  Cl = 2*width/S*sum(vort);

%end




