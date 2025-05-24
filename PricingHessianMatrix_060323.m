% 111022-pricing obj concave proof (ori)-not work
%{
syms P U X V z1 z2 z3 z4 alpha2 d c w_p q k w_w T p h

pri_obj = P*(alpha2.*d.*c.*w_p.*q-k.*w_w.*T.*U)./(k.*w_p.*P+c.*w_p.*q)-c.*X-p.*U-h.*V;

result = hessian(pri_obj,[P,U,X,V]);


B=[z1
  z2
  z3
  z4
  ];

semidef = transpose(B)*result*B

semidef_sim = simplify(semidef)
%}

% 061423-pricing obj concave proof (trick) - san jie duan
%{
syms P1 U1 X1 P2 U2 X2 P3 U3 X3 P4 U4 X4 P5 U5 X5 P6 U6 X6 P7 U7 X7 P8 U8 X8 V1 V2 V3 y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 alpha2 d1 c1 d2 c2 d3 c3 w_p q k w_w T p1 h1 p2 h2 p3 h3

pri_obj_3 = P1*(alpha2.*d1.*c1.*w_p.*q-k.*w_w.*T.*P1)./(k.*w_p.*P1+c1.*w_p.*q)-c1.*X1-p1.*U1-h1.*V1+...
  P2*(alpha2.*d2.*c2.*w_p.*q-k.*w_w.*T.*P2)./(k.*w_p.*P2+c2.*w_p.*q)-c2.*X2-p2.*U2-h2.*V2+...
  P3*(alpha2.*d3.*c3.*w_p.*q-k.*w_w.*T.*P3)./(k.*w_p.*P3+c3.*w_p.*q)-c3.*X3-p3.*U3-h3.*V3;

result_3 = hessian(pri_obj_3,[P1,U1,X1,V1,P2,U2,X2,V2,P3,U3,X3,V3]);

B1=[y1
  y2
  y3
  y4
  y5
  y6
  y7
  y8
  y9
  y10
  y11
  y12
  ];

semidef_fak_3 = transpose(B1)*result_3*B1

semidef_sim_fak_3 = simplify(semidef_fak_3)
%}

% 061423-pricing obj concave proof (trick) - er jie duan-yes
%{
syms P1 U1 X1 P2 U2 X2 P3 U3 X3 P4 U4 X4 P5 U5 X5 P6 U6 X6 P7 U7 X7 P8 U8 X8 V1 V2 y1 y2 y3 y4 y5 y6 y7 y8 alpha2 d1 c1 d2 c2 w_p q k w_w T p1 h1 p2 h2

pri_obj_2b2 = P1*(alpha2.*d1.*c1.*w_p.*q-k.*w_w.*T.*P1)./(k.*w_p.*P1+c1.*w_p.*q)-c1.*X1-p1.*U1-h1.*V1+...
  P2*(alpha2.*d2.*c2.*w_p.*q-k.*w_w.*T.*P2)./(k.*w_p.*P2+c2.*w_p.*q)-c2.*X2-p2.*U2-h2.*V2;

result_2b2 = hessian(pri_obj_2b2,[P1,U1,X1,V1,P2,U2,X2,V2]);

B1=[y1
  y2
  y3
  y4
  y5
  y6
  y7
  y8
  ];

semidef_fak_2b2 = transpose(B1)*result_2b2*B1

semidef_sim_fak_2b2 = simplify(semidef_fak_2b2)
%}

% 111022-pricing obj concave proof (trick)-yi jie duan-yes
%{
syms P U X V y1 y2 y3 y4 alpha2 d c w_p q k w_w T p h

pri_obj1 = P*(alpha2.*d.*c.*w_p.*q-k.*w_w.*T.*P)./(k.*w_p.*P+c.*w_p.*q)-c.*X-p.*U-h.*V;

result1 = hessian(pri_obj1,[P,U,X,V]);


B1=[y1
  y2
  y3
  y4
  ];

semidef_fak = transpose(B1)*result1*B1

semidef_sim_fak = simplify(semidef_fak)
%}

% 040123-nonlinear constraint 1 convex proof (trick)-yes
%{
syms P U1 U2 X y1 y2 y3 y4 alpha2 d c w_p q k w_w T p h

nonli1 = U1-U2+(alpha2.*d.*c.*w_p.*q-k.*w_w.*T.*P)./(k.*w_p.*P+c.*w_p.*q)-X;

res_nonli1 = hessian(nonli1,[P,U1,U2,X]);

A1=[y1
  y2
  y3
  y4
  ];

semidef_nonli1_fak = transpose(A1)*res_nonli1*A1

semidef_sim_nonli1_fak = simplify(semidef_nonli1_fak)
%}

% 040123-nonlinear constraint 2 convex proof (trick)-yes
%{
syms P U smX y1 y2 y3 alpha2 d c w_p q k w_w T p h

nonli2 = U-smX+(alpha2.*d.*c.*w_p.*q-k.*w_w.*T.*P)./(k.*w_p.*P+c.*w_p.*q);

res_nonli2 = hessian(nonli2,[P,U,smX]);

C1=[y1
  y2
  y3
  ];

semidef_nonli2_fak = transpose(C1)*res_nonli2*C1

semidef_sim_nonli2_fak = simplify(semidef_nonli2_fak)
%}

% 111022-concave toy proof
%{
syms x1 x2

obj_conc=-x1^2-(x2+2)^2+4;

result=hessian(obj_conc,[x1,x2]);

B=[z1
  z2];

semidef_conc=transpose(B)*result*B;
semidef_conc_simp=simplify(semidef_conc)
%}

% 111022-convex toy proof
%{ 
syms x1 x2

obj_conv=x1^2+x2^2-1;

result_conv=hessian(obj_conv,[x1,x2]);

B=[z1
  z2];

semidef_conv=transpose(B)*result_conv*B;
semidef_conv_simp=simplify(semidef_conv)
%}