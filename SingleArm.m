% Wavelength- and bandwidth- tunable microring Calculate test(Single Arm)
% Reference: https://github.com/ywj1112/Microring-resonator-design-with-an-agile-and-accurate-simulation-method
%            Shoman, H., Jayatilleka, H., Park, A. H. K., Mistry, A., Jaeger, N. A. F., Shekhar, S., & Chrostowski, L. (2019). Compact wavelength- and bandwidth-tunable microring modulator. Opt Express, 27(19), 26661-26675. doi:10.1364/OE.27.026661
% Author: Zhenyu ZHAO
clear all
lambda = (1525:0.0001:1580)*1e-9;

[Ethru ,Edrop]=RingResonator(lambda, 10e-6 , 1.42*(2*pi*5.9e-6)/4);
figure(1)
hold on 
plot (lambda*1e9, db( abs(Ethru)))

set(gca,'FontSize', 16)
set(gca,'FontName', 'Times New Roman')
ylabel('Relative power (dB)'), xlabel('Wavelength (nm)')


function [Ethru, Edrop]=RingResonator(lambda, r , arm_length)

neff = neff_lambda(lambda);

l1 = pi*r / 2;
l2 = 3*pi*r / 4;
l3 = 3*pi*r / 4;
l_arm = arm_length;
alpha_wg_dB=3;    
alpha_wg=-log(10^(-alpha_wg_dB/10));

A1=exp(-alpha_wg*100*l1);   
A2=exp(-alpha_wg*100*l2);
A3=exp(-alpha_wg*100*l3); 
Ap=exp(-alpha_wg*100*l_arm); 
alpha1 = sqrt(A1);
alpha2 = sqrt(A2);
alpha3 = sqrt(A3);
alphap = sqrt(Ap);


k1=0.1i;
k2=0.1i; 
k3=0.1i; 
k1_conj = -conj(k1);
k2_conj = -conj(k2);
k3_conj = -conj(k3);
t1=sqrt(1-abs(k1)^2);
t2=sqrt(1-abs(k2)^2);
t3=sqrt(1-abs(k3)^2);
t1_conj = conj(t1);
t2_conj = conj(t2);
t3_conj = conj(t3);
theta1=(2*pi./lambda).*neff*(l1);
theta2=(2*pi./lambda).*neff*(l2);
theta3=(2*pi./lambda).*neff*(l3);
%phi1 = (2*pi./lambda).*neff*(l_arm);
phi1 = theta1 + 2*pi/8;


Ethru = (t2*t1*alphap*exp(1i*phi1) + k2*k1_conj*alpha1*exp(1i*theta1)) + ...
    (t2*k1*alphap*exp(1i*phi1) + k2*t1_conj*alpha1*exp(1i*theta1)).* ...
    (t3_conj*t2_conj*k1_conj*alpha1*exp(1i*theta1)*alpha2.*exp(1i*theta2)*alpha3.*exp(1i*theta3) + ...
    t3_conj*k2_conj*t1*alphap*exp(1i*phi1)*alpha2.*exp(1i*theta2)*alpha3.*exp(1i*theta3)) ./ ...
    (1 - t3_conj*t2_conj*t1_conj*alpha1*exp(1i*theta1)*alpha2.*exp(1i*theta2)*alpha3.*exp(1i*theta3)  ...
    - t3_conj*k2_conj*k1*alphap*exp(1i*phi1)*alpha2.*exp(1i*theta2)*alpha3.*exp(1i*theta3));

Edrop= (k3*t2_conj*k1_conj*alpha1*exp(1i*theta1)*alpha2.*exp(1i*theta2) +k3*k2_conj*t1*alphap*exp(1i*phi1)*alpha2.*exp(1i*theta2)) + ...
    (k3*t2_conj*t1_conj*alpha1*exp(1i*theta1)*alpha2.*exp(1i*theta2) +k3*k2_conj*k1*alphap*exp(1i*phi1)*alpha2.*exp(1i*theta2)).* ...
    (t3_conj*t2_conj*k1_conj*alpha1*exp(1i*theta1)*alpha2.*exp(1i*theta2)*alpha3.*exp(1i*theta3) + ...
    t3_conj*k2_conj*t1*alphap*exp(1i*phi1)*alpha2.*exp(1i*theta2)*alpha3.*exp(1i*theta3))./ ...
    (1 - t3_conj*t2_conj*t1_conj*alpha1*exp(1i*theta1)*alpha2.*exp(1i*theta2)*alpha3.*exp(1i*theta3)  ...
    - t3_conj*k2_conj*k1*alphap*exp(1i*phi1)*alpha2.*exp(1i*theta2)*alpha3.*exp(1i*theta3));
end

function [neff]=neff_lambda(lambda)
neff = 2.57 - 0.85*(lambda*1e6-1.55);
end