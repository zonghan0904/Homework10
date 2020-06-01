function out = trajectory(u,P)

% input(1*1):time
% output(6*1):desired trajectory and desired yaw angle

t = u(end);
b1d = [1 0 0];
xd = [0.4*t 0.4*sin(0.2*pi*t) 0.6*cos(0.2*pi*t)];

out = [xd b1d];
end