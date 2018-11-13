function alpha0 = DVM (NACA)

% c = 1
% Uinf = 1

M = floor(NACA/1000);
P = floor((NACA-M*1000)/100);

M = M/100;
P = P/10;

npoints = 1000;

% cosine spacing
dpoints = linspace(pi,0,npoints);
xpoints = 1/2+cos(dpoints)/2;

x1 = xpoints(xpoints < P);
x2 = xpoints(xpoints >= P);

y1 = M/P^2*(2*P*x1-x1.^ 2);
y2 = M/(1-P)^2*(1-2*P+2*P*x2-x2.^2);

zpoints = [y1,y2];

nel = npoints-1;

I = ones(1,nel);

deltax = xpoints(2:nel+1)-xpoints(1:nel);
deltaz = zpoints(2:nel+1)-zpoints(1:nel);
l = (deltax.^2+deltaz.^2).^0.5;

n = [-deltaz./l;deltax./l]';

xcp = xpoints(1:nel)+3/4*(deltax);
xvp = xpoints(1:nel)+1/4*(deltax);

zcp = zpoints(1:nel)+3/4*(deltaz);
zvp = zpoints(1:nel)+1/4*(deltaz);

[XP,XC] = meshgrid(xvp,xcp);
[ZP,ZC] = meshgrid(zvp,zcp);

R2 = (XC-XP).^2+(ZC-ZP).^2;

N1 = n(:,1)*ones(1,nel);
N2 = n(:,2)*ones(1,nel);

A = 1./(2*pi*R2).*((ZC-ZP).*N1-(XC-XP).*N2);

gamma = ones(1,nel)*inv(A)*n;

alpha0 = atan(-gamma(1)/gamma(2));

end



