%------------------------------------------------------------------------------
%
%                          LAMBERBATTIN
%
%   this subroutine solves Lambert's problem using Battins method. The method
%   is developed in Battin (1987). It uses continued fractions to speed the
%   solution and has several parameters that are defined differently than
%   the traditional Gaussian technique.
%
% Inputs:         Description                    Range/Units
%   ro          - IJK Position vector 1          km
%   r           - IJK Position vector 2          km
%   dm          - direction of motion            'pro','retro'
%   Dtsec       - Time between ro and r          s
%
% OutPuts:
%   vo          - IJK Velocity vector            km/s
%   v           - IJK Velocity vector            km/s
%
% Reference:
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill
% , New York; 4th edition(2013).
% 
% Last modified:   2023/12/09   Meysam Mahooti
% 
%------------------------------------------------------------------------------  
function [vo, v] = LAMBERTBATTIN(ro, r, dm, Dtsec)

small = 1.0e-10;
gm = 398600.4415;   % km3/s2
y1 = 0;
magr = norm(r);
magro = norm(ro);
CosDeltaNu= dot(ro,r)/(magro*magr);
rcrossr = cross(ro,r);
magrcrossr = norm(rcrossr);
if strcmp(dm,'pro')
    SinDeltaNu = magrcrossr/(magro*magr);
else
    SinDeltaNu = -magrcrossr/(magro*magr);
end
DNu = atan2(SinDeltaNu,CosDeltaNu);

% the angle needs to be positive to work for the long way
if DNu < 0.0
    DNu = 2.0*pi+DNu;
end

RoR   = magr/magro;
eps   = RoR - 1.0;
tan2w = 0.25*eps*eps / ( sqrt( RoR ) + RoR * ...
                       ( 2.0 + sqrt( RoR ) ) );
rp    = sqrt( magro*magr )*( (cos(DNu*0.25))^2 + tan2w );
if ( DNu < pi )
    L = ( (sin(DNu*0.25))^2 + tan2w ) / ( (sin(DNu*0.25))^2 + tan2w + cos( DNu*0.5 ) );
else
    L = ( (cos(DNu*0.25))^2 + tan2w - cos( DNu*0.5 ) ) / ( (cos(DNu*0.25))^2 + tan2w );
end
m    = gm*Dtsec*Dtsec / ( 8.0*rp*rp*rp );
x    = 10.0;
xn   = L;
chord= sqrt( magro*magro + magr*magr - 2.0*magro*magr*cos( DNu ) );
s    = ( magro + magr + chord )*0.5;
lim1 = sqrt(m/L);

Loops= 1;
while ((abs(xn-x) >= small) && (Loops <= 30))
    x    = xn;    
    tempx= seebatt(x);
    Denom= 1.0 / ( (1.0+2.0*x+L) * (4.0*x + tempx*(3.0+x) ) );
    h1   = ( L+x )^2 * ( 1.0+ 3.0*x + tempx )*Denom;
    h2   = m*( x - L + tempx )*Denom;
    
    % ----------------------- Evaluate CUBIC ------------------
    b = 0.25*27.0*h2 / ( (1.0+h1)^3 );
    if (b < -1.0) % reset the initial condition
        xn = 1.0 - 2.0*L;
    else
        if (y1 > lim1)
            xn = xn * (lim1/y1);
        else
            u = 0.5*b / ( 1.0 + sqrt( 1.0 + b ) );            
            k2 = seebattk(u);
            y = ( ( 1.0+h1 ) / 3.0 )*( 2.0 + sqrt( 1.0+b )/( 1.0+2.0*u*k2*k2 ) );
            xn= sqrt( ( (1.0-L)*0.5 )^2 + m/(y*y) ) - ( 1.0+L )*0.5;
        end
    end
    Loops = Loops + 1;
    y1 = sqrt(m/((L+x)*(1.0+x)));
end

a = gm*Dtsec*Dtsec / (16.0*rp*rp*xn*y*y );

% ------------------ Find Eccentric anomalies -----------------
% ------------------------ Hyperbolic -------------------------
if ( a < -small )
    arg1 = sqrt( s / ( -2.0*a ) );
    arg2 = sqrt( ( s-chord ) / ( -2.0*a ) );
    % ------- Evaluate f and g functions --------
    AlpH = 2.0 * asinh( arg1 );
    BetH = 2.0 * asinh( arg2 );
    DH   = AlpH - BetH;
    F    = 1.0 - (a/magro)*(1.0 - cosh(DH) );
    GDot = 1.0 - (a/magr) *(1.0 - cosh(DH) );
    G    = Dtsec - sqrt(-a*a*a/gm)*(sinh(DH)-DH);
else
    % ------------------------ Elliptical ---------------------
    if ( a > small )
        arg1 = sqrt( s / ( 2.0*a ) );
        arg2 = sqrt( ( s-chord ) / ( 2.0*a ) );
        Sinv = arg2;
        BetE = 2.0*asin(Sinv);
        if ( DNu > pi )
            BetE= -BetE;
        end
        Sinv= arg1;
        am  = s*0.5;
        ae  = pi;
        be  = 2.0*asin( sqrt( (s-chord)/s ) );
        tm  = sqrt(am*am*am/gm)*(ae - (be-sin(be)));
        if ( Dtsec > tm )
            AlpE= 2.0*pi-2.0*asin( Sinv );
        else
            AlpE= 2.0*asin( Sinv );
        end
        DE  = AlpE - BetE;
        F   = 1.0 - (a/magro)*(1.0 - cos(DE) );
        GDot= 1.0 - (a/magr)* (1.0 - cos(DE) );
        G   = Dtsec - sqrt(a*a*a/gm)*(DE - sin(DE));
    else
        % --------------------- Parabolic ---------------------
        arg1 = 0.0;
        arg2 = 0.0;
		error(' a parabolic orbit ');
    end
end

vo = zeros(1,3);
v  = zeros(1,3);
for i= 1:3
    vo(i)= ( r(i) - F*ro(i) )/G;
    v(i) = ( GDot*r(i) - ro(i) )/G;
end

