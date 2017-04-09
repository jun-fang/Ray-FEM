function [p,t] = annulusmesh(x,y,r1,r2,h)
% r1: inner radius
% r2: outer radius

[p,t] = odtmesh2d(@fd,@huniform,h,[-1,-1;1,1],[],0,x,y,r1,r2);

    function s = fd(p,x,y,r1,r2)
    s = (sqrt(sum((p(:,1)-x).^2+(p(:,2)-y).^2,2))-r1)...
        .*(sqrt(sum((p(:,1)-x).^2+(p(:,2)-y).^2,2))-r2);
    end

    function h = huniform(p,varargin)
    h = ones(size(p,1),1);
    end
end