function [A,phase]=conducting_sphere(cond_layer1,r_layerbot1)

global R_call ind_period r_layertop mu0


r0k=(1.0+1i)*sqrt(0.5*mu0*cond_layer1*2*pi/ind_period).*r_layertop; % no dimension

r1k=(1.0+1i)*sqrt(0.5*mu0*cond_layer1*2*pi/ind_period).*r_layerbot1; % no dimension

realdif = real(r0k)-real(r1k);
      c2r0k   = 3.0./r0k.^2-1.0;
      c1r0k   =-3.0./r0k;
      c2r1k   = 3.0./r1k.^2-1.0;
      c1r1k   =-3.0./r1k;
      coshyp  =1.0+exp(-2.0*realdif);
      sinhyp  =1.0-exp(-2.0*realdif);
      term1   =sin(realdif).*coshyp+1i* cos(realdif).*sinhyp;
      term2   =cos(realdif).*coshyp-1i*sin(realdif).*sinhyp;
      xi      =(c2r1k.*term2+c1r1k.*term1)./(c2r1k.*term1-c1r1k.*term2);
      ZZ      =-0.5*(c2r0k+c1r0k.*xi)*(r_layertop/R_call).^3;
      ZZ_aprox=1/6*((1.0+1i)*sqrt(0.5*mu0*cond_layer1*2*pi/ind_period)).^2.*r_layerbot1.*(r_layertop-r_layerbot1);
      A       = 2.0*abs(ZZ);
      phase   = atan2(imag(ZZ),real(ZZ))*180/pi;
  
      

