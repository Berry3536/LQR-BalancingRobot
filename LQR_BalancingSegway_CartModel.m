%{
* 
* Updated  : 2023-08-22
* Maker    : Janguk Kim
* Filename : LQR_BalancingSegway_CartModel.m
* Purpose  : LQR optimization of Segway Platform
*
%}
format compact;
clc;clear;

%% Segway system variable
M = 0.550;     % mass of the wheels and shaft (kg)
m = 0.345;     % mass of the pendulum (kg)
b = 0.1;       % estimate of viscous friction coefficient (N-m-s)
I = 0.00053084;% moment of inertia of the pendulum (kg*m^2)
g = 9.807;     % acceleration due to gravity (m/s^2)
l = 0.046;     % length to pendulum center of mass (m)


%% Linear state equation from system modeling
% Linear state equation: dx/dt = Ax+Bu, y = Cx+Du

% denominator for the A and B matrices
p = I*(M+m)+M*m*l^2;

% A matrix of dx/dt = Ax+Bu
A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];

% B matrix of dx/dt = Ax+Bu
B = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];

% C matrix of y = Cx+Du
C = [1 0 0 0;
     0 0 1 0];

% D matrix of y = Cx+Du
D = [0;
     0];

%% Fine tuning of weight matrices

% Q matrix of quadratic cost function (state weight matrix)
Q = [10000     0     0     0
     0     0     0     0
     0     0     10000     0
     0     0     0     0];
 
% R matrix of quadratic cost function (input weight matrix)
R = 1;

%% Optimization

% Calculate K gain of optimal control input u*=-Kx
K = lqr(A,B,Q,R);
