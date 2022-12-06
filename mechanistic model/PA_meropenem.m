function dydt= PA_meropenem(t,y,parameter)
dydt=zeros(size(y));

r = parameter(1) ; 
r_tilde = parameter(2);
phi = parameter(3);
alpha = parameter(4);
gamma = parameter(5);
delta = parameter(6);
psi = parameter(7);
rho = parameter(8) ;
rho_deg = parameter(9);
T_50 = parameter(10) ;
A_50 = parameter(11) ;

P=y(1);
H=y(2);
S=y(3);
A=y(4);
N=y(5);

R1 = r*N*H -(gamma*A/(T_50+A))*H + delta*S - (rho*A/(A_50+A))*H - phi*H;
R2 = (gamma*A/(T_50+A))*H - delta*S - psi*S;
R3 = -alpha*A - (rho_deg*A/(A_50+A))*H;
R4 = - r_tilde*N*H  ;


dydt = [(R1 + 0.5*R2) R1 R2 R3 R4];
end

