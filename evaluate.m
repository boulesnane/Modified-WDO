function  [y]=evaluate(posx)
global param
%Sphere Objective Function
x=0;
for i=1:param.npar

     x=x+(posx(i))^2;
end
y=x;
end