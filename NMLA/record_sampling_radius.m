%% record the raduis of sampling circle of NMLA
% varying with two parameters: 1) estimate radius Rest, 2) frequency omega
% the polynomial plays a role as well

high_omega = [40 80 160 320 640]*pi;     % high frequency
low_omega = sqrt(high_omega);            % low frequency
high_wl = 2*pi./high_omega; 
low_wl = 2*pi./low_omega;
high_pml = high_wl;
low_pml = low_wl;


sd = 1/2;     
Rest = sqrt(2)*sd;

high_r = NMLA_radius(high_omega,Rest);
md = sd + high_r + high_pml;
md = ceil(md*10)/10;

Rest = sqrt(2)*md;
low_r = NMLA_radius(low_omega,Rest);
ld = md + low_r + low_pml;
ld = ceil(ld*10)/10;



fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( 'NMLA radius of sampling circle:\n\n');
fprintf( 'High frequency omega/2pi:   ');
fprintf( '&  %1.0f  ',high_omega/(2*pi));
fprintf( '\n\nRadius for high freq:   ');
fprintf( '&  %1.2f  ',high_r);

fprintf( '\n\nMiddle domain size:     ');
fprintf( '&  %1.2f  ',md);

fprintf( '\n\nRadius for low freq:    ');
fprintf( '&  %1.2f  ',low_r);
fprintf( '\n\nLarge domain size:      ');
fprintf( '&  %1.2f  ',ld);

fprintf( ['\n' '-'*ones(1,80) '\n']);



