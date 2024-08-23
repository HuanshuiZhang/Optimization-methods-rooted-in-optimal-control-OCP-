clc;
clear
addpath('..\OCP')

x0 = 15;
R = 1;
N = 50;
a1 = 1; a2=1; a3 = 1;
[X,FVAL,EXITFLAG,GRAD,HESSIAN] = OCP(@(x)myfun(x,a1),@(x)mygradfun(x,a2),@(x)myhessfun(x,a3),x0,R,N,1e-15,'Method1','on');