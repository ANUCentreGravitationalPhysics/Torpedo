% adapted from Calum Torrie's ssmake4p by Ken Strain 3/01
% yaw and roll bugs fixed?  
% 
% 2/14/02, Mark Barton: 
%   corrected quad matrix elements derived by generalizing 
%   Calum's triple matrix elements using Mathematica
%   rename to ssmake4pv2eMB
% 2/24/05, Mark Barton: 
%   new matrix elements with more realistic blades
%   generated from scratch in Mathematica
%   rename to ssmake4pv2eMB2
% 3/30/05, Mark Barton: 
%   adjust sign of certain elements to match
%   conventions in Calum's original code 

global pend
quadopt;
% general redefintion of variables to simplify notation
% and cope with different pendulum designs
g = pend.g;
kcn	=  1/2 * (2*pi*pend.ufcn)^2*pend.mn;
kc1	=  1/2 * (2*pi*pend.ufc1)^2*pend.m1;
kc2	=  1/2 * (2*pi*pend.ufc2)^2*pend.m2;
ln   = pend.ln;
l1   = pend.l1;
l2   = pend.l2;
l3   = pend.l3;
kwn	=pend.Y1*pi*pend.rn^2/ln*pend.nwn/2;
kw1	=pend.Y1*pi*pend.r1^2/l1*pend.nw1/2;	
kw2	=pend.Y2*pi*pend.r2^2/l2*pend.nw2/2;	
kw3	=pend.Y3*pi*pend.r3^2/l3*pend.nw3/2;	

%********************************************************************
% allows choice of 2 wires to set separation to zero
sm = 0; % separation of top wires at structure
sn = 0; % separation of top wires at mass N

s0=pend.su;
s1=pend.su;
s2=pend.si;
s3=pend.si;
s4=pend.sl;
s5=pend.sl;

dm=pend.dm;
dn=pend.dn;
d0=pend.d0;
d1=pend.d1;
d2=pend.d2;
d3=pend.d3;
d4=pend.d4;


%********************************************************************
mn = pend.mn;
m1 = pend.m1;
m2 = pend.m2;
m3 = pend.m3;
Inx = pend.Inx;
I1x = pend.I1x;
I2x = pend.I2x;
I3x = pend.I3x;
Iny = pend.Iny;
I1y = pend.I1y;
I2y = pend.I2y;
I3y = pend.I3y;
Inz = pend.Inz;
I1z = pend.I1z;
I2z = pend.I2z;
I3z = pend.I3z;

mn3  = mn+m1+m2+m3;
m13	= m1+m2+m3;
m23	= m2+m3;

nn0 = pend.nn0;
nn1 = pend.nn1;
n0  = pend.n0;
n1  = pend.n1;
n2  = pend.n2;
n3  = pend.n3;
n4  = pend.n4;
n5  = pend.n5;

%***********************************************************************************
% cosine and sine of the angle the wire makes with the vertical (z)

sin=(nn1-nn0)/ln;			%sin(omegan)
si1=(n1-n0)/l1;			%sin(omega1)
si2=(n3-n2)/l2;			%sin(omega2)
si3=(n5-n4)/l3;			%sin(omega3)	

cn=(ln^2-(nn1-nn0)^2)^0.5/ln;	%cos(omegan)
c1=(l1^2-(n1-n0)^2)^0.5/l1;	%cos(omega1)
c2=(l2^2-(n3-n2)^2)^0.5/l2;	%cos(omega2)
c3=(l3^2-(n5-n4)^2)^0.5/l3;	%cos(omega3)

%********************************************************
bd = pend.bd;
%********************************************************************

kn = kwn;
k1 = kw1;
k2 = kw2;
k3 = kw3;
su=pend.su;
si=pend.si;
sl=pend.sl;

symbexport4 % the output of the Mathematica

FFFLP = diag([-1 1 -1 1 -1 1 -1 1]); % allow for different sign conventions

% LONGITUDINAL-PITCH
AAALP = [...
    zeros(8) eye(8)
    FFFLP'*(-kmlp)\(xmlp-cqxmlp'/qm*cqxmlp)*FFFLP -bd*eye(8)
];

bmlp = FFFLP'*(-kmlp)\(cxsmlp-cqxmlp'/qm*cqsm); % ground displacement inputs: x,y,z,yaw,pitch,roll
bmlp2 = FFFLP'*kmlp\eye(8); % coordinate force inputs: pitchn,xn,pitch1,x1,...

%NOTE: -ve signs in last four entries below are a fudge to match the sign
%of the old matrix elements - needs triple-checking.
BBBLP = [...
    zeros(8,9)
    bmlp(:,1) bmlp2(:,2) bmlp2(:,4) bmlp2(:,6) bmlp2(:,8)  -bmlp2(:,1) -bmlp2(:,3) -bmlp2(:,5) -bmlp2(:,7)
];

CCCLP = [...
 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 %l n
 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 %l 
 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0
 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 %long vel
 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 %pitch vel

 ];
 
 DDDLP = [...
 0 0 0 0 0 0 0 0 0 
 0 0 0 0 0 0 0 0 0 
 0 0 0 0 0 0 0 0 0 
 0 0 0 0 0 0 0 0 0 
 0 0 0 0 0 0 0 0 0 
 0 0 0 0 0 0 0 0 0 
 0 0 0 0 0 0 0 0 0 
 0 0 0 0 0 0 0 0 0 
 0 0 0 0 0 0 0 0 0 
 0 0 0 0 0 0 0 0 0 
 ];

lpe = ss(AAALP,BBBLP,CCCLP,DDDLP);

% TRANSVERSE-ROLL
FFFTR = diag([1 1 1 1 1 1 1 1]); % allow for different sign conventions

AAATR = [...
    zeros(8) eye(8)
    FFFTR'*(-kmtr)\(xmtr-cqxmtr'/qm*cqxmtr)*FFFTR -bd*eye(8)
];

bmtr = FFFTR'*(-kmtr)\(cxsmtr-cqxmtr'/qm*cqsm); % ground displacement inputs: x,y,z,yaw,pitch,roll
bmtr2 = FFFTR'*kmtr\eye(8); % coordinate force inputs: rolln,yn,roll1,y1,...

BBBTR = [...
    zeros(8,6)
    bmtr(:,2) bmtr2(:,2) bmtr2(:,4) bmtr(:,6) bmtr2(:,1)  bmtr2(:,3)
];

CCCTR = [...
    0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0 
    0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0
    0	0 	0	0	0	1	0	0	0	0	0	0	0	0	0	0
    0	0 	0	0	0	0	0	1	0	0	0	0	0	0	0	0
    1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0
    0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0  
    0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0
];
 
DDDTR = [...
    0   0  0  0  0  0
    0   0  0  0  0  0
    0   0  0  0  0  0
    0   0  0  0  0  0   
    0   0  0  0  0  0     
    0   0  0  0  0  0     
    0   0  0  0  0  0     
    0   0  0  0  0  0     
  	0   0  0  0  0  0     
    0   0  0  0  0  0     
];

trer = ss(AAATR,BBBTR,CCCTR,DDDTR);   

%YAW

FFFY = diag([1 1 1 1]); % allow for different sign conventions

AAAY = [...
    zeros(4) eye(4)
    FFFY'*(-kmy)\(xmy-cqxmy'/qm*cqxmy)*FFFY -bd*eye(4)
];
bmy = FFFY'*(-kmy)\(cxsmy-cqxmy'/qm*cqsm); % ground displacement inputs: x,y,z,yaw,pitch,roll
bmy2 = FFFY'*kmy\eye(4); % coordinate force inputs: yawn,yaw1,...

BBBY = [...
    zeros(4,5)
    bmy(:,4) bmy2(:,1) bmy2(:,2) bmy2(:,3) bmy2(:,4)
];

CCCY = [...
    1	0	0	0	0	0 0 0	
   	0	1	0	0	0	0 0 0	
   	0	0	1	0	0	0 0 0
    0	0	0 	1	0	0 0 0
    0   0   0   0   1   0 0 0
];

DDDY = [
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0
    0 0 0 0 0  
];

ye = ss(AAAY,BBBY,CCCY,DDDY);


% VERTICAL
FFFV = diag([1 1 1 1]); % allow for different sign conventions

AAAV = [...
    zeros(4) eye(4)
    FFFV'*(-kmv)\(xmv-cqxmv'/qm*cqxmv)*FFFV -bd*eye(4)
];

bmv = FFFV'*(-kmv)\(cxsmv-cqxmv'/qm*cqsm); % ground displacement inputs: x,y,z,yaw,pitch,roll
bmv2 = FFFV'*kmv\eye(4); % coordinate force inputs: yawn,yaw1,...

BBBV = [...
    zeros(4,3)
    bmv(:,3) bmv2(:,1) bmv2(:,2)
];

CCCV = [
    1	0	0 0 0	0	0	0	
    0	1	0 0 0	0	0	0	
    0	0	1 0 0	0	0	0	
    0	0	0 1 0	0	0	0
    0	0	0 0 1	0	0	0
];
   
DDDV = [
    0 0 0 
    0 0 0 
    0 0 0  
    0 0 0 
    0 0 0
];

ve = ss(AAAV,BBBV,CCCV,DDDV);

