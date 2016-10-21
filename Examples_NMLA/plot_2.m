cp_omega = [10 20 40 60 80 120 160]*pi;
Rest = 2;
low_r = NMLA_radius(sqrt(cp_omega),Rest)
high_r = NMLA_radius(cp_omega,Rest)
M_a = sm_a + ceil(high_r/0.1)*0.1
L_a = sm_a + M_a + ceil(low_r/0.1)*0.1


%% Show convergence rate with respect to omega
if (0)
    % angle error
    figure(21);
    loglog(rec_omega/(2*pi), rec_ang_err1,'bs-');
    hold on;
    loglog(rec_omega/(2*pi), rec_ang_err2,'r*-');
    axis([0 40 -inf inf])% axis tight;
    xlabel('Frequency \omega/2\pi');
    ylabel('Relative L^2 error');
    legend('Angle Error 1','Angle Error 2','LOCATION','Best');
    
    % NR-FEM error
    figure(22);
    loglog(rec_omega/(2*pi), rec_NR_err1,'bs-');
    hold on;
    loglog(rec_omega/(2*pi), rec_NR_err2,'r*-');
    hold on;
    loglog(rec_omega/(2*pi), rec_ER_err,'ko-');
    axis([0 40 -inf inf]);
    xlabel('Frequency \omega/2\pi');
    ylabel('Relative L^2 error');
    legend('NR-FEM Error 1','NR-FEM Error 2','ER-FEM Error','LOCATION','Best');
    
    
end

%% Show convergence rate with respect to h
if (0)
    figure(5);
    %     loglog(rec_h, rec_NR_err1,'bs-');
    %     hold on;
    loglog(1./rec_h, rec_NR_err2,'r*-');
    hold on;
    loglog(1./rec_h, rec_ER_err,'ko-');
    hold on;
    loglog(1./rec_h, rec_int_err,'g^-');
    axis([50 1000 -inf inf]);
    xlabel('mesh size 1/h');
    ylabel('Relative L^2 error');
    legend('NR-FEM Error 2','ER-FEM Error','Interpoltation Error','LOCATION','Best');
end

