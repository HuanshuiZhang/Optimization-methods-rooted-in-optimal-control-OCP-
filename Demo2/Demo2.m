clc;
clear
close all
addpath('..\OCP')

x0 = [15;20];

N = 100;
paras = [1; 2; 3];

R = 1* eye(2);
[X,FVAL,EXITFLAG,GRAD,HESSIAN] = OCP('myfun','mygradfun','myhessfun',x0,R,N,1e-15,'Method1','on',paras);

R = 0.1* eye(2);
[X,FVAL,EXITFLAG,GRAD,HESSIAN] = OCP('myfun','mygradfun','myhessfun',x0,R,N,1e-15,'Method2','on',paras);