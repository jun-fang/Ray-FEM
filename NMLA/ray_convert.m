function ray = ray_convert(ray_input,opt)
re = real(ray_input);
im = imag(ray_input);
ray = atan2(im,re);
if opt == 1   %% ray angle
    ray = ray + 2*pi*(ray < 0);
elseif opt == 2   %% complex form
    ray = exp(1i*ray);
end
    