function z=projc(x,y,u)
if abs(x)<=u*abs(y)
    z=x;
else
    if x>u*abs(y)
    z=u*abs(y);
    else 
        z=-u*abs(y);
    end
end
end
% function z=projc(x,y,u)
% if y==0
%     z=0;
% else
% if abs(x)<=u*abs(y)
%     z=x/abs(y);
% else
%     if x>u*abs(y)
%     z=u;
%     else 
%         z=-u;
%     end
% end
% end
% end