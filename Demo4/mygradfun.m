function F = mygradfun(x,varagin)
F = [ 2*x(1) - 400*x(1)*(-x(1)^2 + x(2)) - 2;-200*x(1)^2 + 200*x(2)];