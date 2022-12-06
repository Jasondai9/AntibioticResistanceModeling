import numpy as np
import sys
import matplotlib.pyplot as plt
import math

"""
V. anguillarum bacterial growth in G. mellonella (moth) larvae with TET (Tetracycline) as the antibiotic
https://doi.org/10.1371%2Fjournal.pcbi.1008037

Parameter	Definition	Value
r	Average Replication rate of bacteria	0.4779
m	Co-efficient for the host immune response	0.6772
n	Hill co-efficient in the immune response	0.9193
v	Standard deviation for host heterogeneity	0.0525 (unused)
a1	Maximum kill rate of antibiotic	0.7281
a2	Level of antibiotic giving half max kill rate	0.1910
k	Hill co-efficient in AB induced death.	2.9821
a	Decay rate of antibiotic (half-life = 5.9hrs)	0.1174
Bdead	Bacterial load at which the host dies [38]	109
"""
def simulate_simple_markov(parameters_np, tot_time, tau, init_bac, init_antibiotic):
    r, m, n, v, a1, a2, k, a, Bdead = list(parameters_np)

    # initial values
    B = init_bac
    A = init_antibiotic
    t = 0

    # track number of bacteria over time
    bacteria_list = []
    bacteria_list.append(B)

    while t < tot_time: # for all time

        # maximum cap per organism
        if B >= Bdead:
            B = Bdead
        
        elif B > 0:
            # R+
            rate_of_replication = r*B
            
            # R-
            rate_of_immune_death = m*np.power(B, n)
            rate_of_antibiotic_death = (a1*B*np.power(A, k)) / (np.power(A, k)+np.power(a2, k))
            rate_of_death = rate_of_immune_death + rate_of_antibiotic_death

            # take a step in the markov chain
            B = B + np.random.poisson(tau * rate_of_replication) - np.random.poisson(tau * rate_of_death)
            A = A + np.random.poisson(tau * a * A)
        else:
            B = 0
        t += tau
        bacteria_list.append(B)
    
    return bacteria_list


def simple_markov():
    parameters_np = np.loadtxt(sys.argv[1])
    num_runs = 100
    total_time = 48
    tau = 0.25
    initial_bac = 10e5
    dosages = [0, 0.1, 0.2, 0.3, 0.4]
    plt.figure(figsize=(10, 6))

    # for all initial antibody dosages
    for init_A in dosages:
        # average over num_runs
        all_runs = []
        for _ in range(num_runs):
            bacteria_list = simulate_simple_markov(parameters_np, tot_time=total_time, tau=tau, init_bac=initial_bac, init_antibiotic=init_A)
            all_runs.append(bacteria_list)
        a = np.array(all_runs)
        means = np.mean(a, axis=0)
        stdevs = np.std(a, axis=0)

    
        plt.plot(means, label=f'{np.round(init_A, 1)}')
        plt.fill_between(np.arange(total_time/tau +1), means-stdevs, means+stdevs, alpha=0.1)
    
    plt.title(f'Bacteria population over time')
    plt.legend(title='Antibiotic dosage (mg)')
    plt.ylabel('Population')
    plt.xlabel(f'Time (h)')
    plt.yscale('log')
    plt.xticks(np.arange(stop=total_time/tau +1, step=int(1/tau)*2), np.arange(stop=total_time+1, step=2))
    plt.savefig(f'plots/simple_markov.png')

    return


"""
Matlab code for mechanistic ODE model
P. aeruginosa bacterial growth with β-lactam antibiotic meropenem
https://doi.org/10.1371%2Fjournal.pcbi.1006012

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
"""
def mechanistic_markov():
    pass


def main():
    simple_markov()

if __name__ == "__main__":
    main()