function rays_comp = post_compressor(rays, pct)
%% compress closed-by ray directions to reduce redundancy
% and diminsh ill-conditioning of the assembling matrix

% rays: 1xNray
Nray = length(rays);

if Nray == 1
    rays_comp = rays;
end

if Nray == 2
    dray = abs(rays(1) - rays(2));
    if dray < pct
        rays_comp = mean(rays)/abs(means(rays));
    else
        rays_comp = rays;
    end
end

if Nray == 3
    ray1 = rays;  rays_comp = rays;
    ray2 = rays(2:end);  ray2 = [ray2, rays(1)];
    dray = abs(ray1 - ray2);
    
    idx = find(dray < pct);
    n = length(idx);
    
    if n == 0
        rays_comp = rays;
    elseif n >= Nray - 1
        ave = mean(rays);
        rays_comp = ave/abs(ave);
    else
        ray_mid = ( ray1(idx) + ray2(idx) )/2;
        ray_mid = ray_mid./abs(ray_mid);
        switch idx
            case 1 
                rays_comp = [ray_mid, rays(3)];
            case 2 
                rays_comp = [rays(1), ray_mid];
            case 3 
                rays_comp = [ray_mid, rays(2)];
        end
    end
    
    if length(rays_comp) < Nray
        rays_comp = post_compressor(rays_comp, pct);
    end
end


if Nray == 4
    ray1 = rays;  rays_comp = rays;
    ray2 = rays(2:end);  ray2 = [ray2, rays(1)];
    dray = abs(ray1 - ray2);
    
    
    idx = find(dray < pct);
    n = length(idx);
    
    if n == 0
        rays_comp = rays;
    elseif n >= Nray - 1
        ave = mean(rays);
        rays_comp = ave/abs(ave);
    elseif n == 1        
        ray_mid = ( ray1(idx) + ray2(idx) )/2;
        ray_mid = ray_mid./abs(ray_mid);
        
        switch idx
            case 1 
                rays_comp = [ray_mid, rays(3), rays(4)];
            case 2 
                rays_comp = [rays(1), ray_mid, rays(4)];
            case 3 
                rays_comp = [rays(1), rays(2), ray_mid];
            case 4 
                rays_comp = [ray_mid, rays(2), rays(3)];
        end
    else
        ray_mid = ( ray1(idx) + ray2(idx) )/2;
        ray_mid = ray_mid./abs(ray_mid);
        
        if sum(abs(idx - [1,2])) == 0
            rays_comp = [ray_mid, rays(4)];
        elseif sum(abs(idx - [1,3])) == 0
            rays_comp = ray_mid;
        elseif sum(abs(idx - [1,4])) == 0
            rays_comp = [ray_mid(1), rays(3), ray_mid(2)];
        elseif sum(abs(idx - [2,3])) == 0
            rays_comp = [rays(1), ray_mid];
        elseif sum(abs(idx - [2,4])) == 0
            rays_comp = ray_mid;
        elseif sum(abs(idx - [3,4])) == 0
            rays_comp = [ray_mid(2), rays(2), ray_mid(1)];
        end
    end
    
    
    if length(rays_comp) < Nray
        rays_comp = post_compressor(rays_comp, pct);
    end  
end
        
     
end

