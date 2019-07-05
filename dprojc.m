function z=dprojc(x,y,u)
if abs(x)<=u*abs(y)
    z=1;
else
    z=0;
end
end