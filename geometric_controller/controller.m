function out = controller(u,P)

% input(25*1):desired trajectory and full state feedback, x v R Omega time
% output(4*1): force and moment control input

% process inputs
xd    = u(1:3);
b1d   = u(4:6);

% current state
x     = u(7:9);
v     = u(10:12);
R     = reshape(u(13:21),3,3);
Omega = u(22:24);
t     = u(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% my implementatoin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ex = x - xd;
ev = v - [0; 0; 0]; % because the desire position won't change, so xd_dot = 0
e3 = [0; 0; 1];
xd_ddot = [0; 0; 0]; % because the desire position won't change, so xd_ddot = 0

temp = (-P.kx * ex - P.kv * ev - P.mass * P.gravity * e3 + P.mass * xd_ddot);
b3c = -temp / norm(temp);
temp1 = cross(b3c, b1d);
b1c = -cross(b3c, temp1) / norm(temp1);
b2c = cross(b3c, b1c);

Rc = [b1c(1) b1c(2) b1c(3); b2c(1) b2c(2) b2c(3); b3c(1) b3c(2) b3c(3)];
Rc_dot = [0 0 0; 0 0 0; 0 0 0]; % not sure
Omega_c = Omega - R.'* Rc * vee(Rc.' * Rc_dot); % not sure

eR = vee((Rc.' * R - R.' * Rc)) / 2;
eOmega = Omega - R.' * Rc * Omega_c;

J = [P.Jxx 0 0; 0 P.Jyy 0; 0 0 P.Jzz]; 

f = -(-P.kx * ex - P.kv * ev - P.mass * P.gravity * e3 + P.mass * xd_ddot).'* R * e3;
M = -P.kR * eR - P.kOmega * eOmega + cross(Omega, J * Omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f = 0;
% M = [0; 0; 0];

out = [f;M];
end