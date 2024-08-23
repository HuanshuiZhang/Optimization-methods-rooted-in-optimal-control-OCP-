clc;
clear
close all
addpath('..\OCP')

x0 = [10;5];

N = 100;
paras = [1; 2; 3];

R = 1* eye(2);
[X,FVAL,EXITFLAG,GRAD,HESSIAN] = OCP('myfun','mygradfun','myhessfun',x0,R,N,1e-15,'Method1','on',paras);

R = diag([0.03,0.00008]);
[X,FVAL,EXITFLAG,GRAD,HESSIAN] = OCP('myfun','mygradfun','myhessfun',x0,R,N,1e-15,'Method2','on',paras);
