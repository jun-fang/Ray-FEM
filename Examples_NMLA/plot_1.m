%% Show convergence rate with respect to omega
if (1)
    % stability constant
    figure(1);
    plot(rec_omega/(2*pi), rec_sta_con,'*-');
    axis([0 200 0.35 0.45])% axis tight;
    xlabel('Frequency \omega/2\pi');
    ylabel('Stability constant C_{sta}');
    
    % angle error
    figure(2);
    loglog(rec_omega/(2*pi), rec_ang_err1,'bs-');
    hold on;
    loglog(rec_omega/(2*pi), rec_ang_err2,'r*-');
    axis([0 200 -inf inf])% axis tight;
    xlabel('Frequency \omega/2\pi');
    ylabel('Relative L^2 error');
    legend('Angle Error 1','Angle Error 2','LOCATION','Best');
    
    % NR-FEM error
    figure(3);
    loglog(rec_omega/(2*pi), rec_NR_err1,'bs-');
    hold on;
    loglog(rec_omega/(2*pi), rec_NR_err2,'r*-');
    hold on;
    loglog(rec_omega/(2*pi), rec_ER_err,'ko-');
    hold on;
    loglog(rec_omega/(2*pi), rec_int_err,'g^-');
    axis([0 200 -inf inf]);
    xlabel('Frequency \omega/2\pi');
    ylabel('Relative L^2 error');
    legend('NR-FEM Error 1','NR-FEM Error 2','ER-FEM Error','Interpoltation Error','LOCATION','Best');
    
    % optimality constant
    figure(4);
    plot(rec_omega/(2*pi), rec_NR_err2./rec_int_err,'*-');
    axis([0 200 0.35 0.55])% axis tight;
    xlabel('Frequency \omega/2\pi');
    ylabel('Optimality relation C_{opt}');
end

%% Show convergence rate with respect to h
if (1)
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

