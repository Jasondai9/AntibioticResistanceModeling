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
Bdead	Bacterial load at which the host dies [38]	10e9
"""
def simulate_simple_markov(parameters_np, tot_time, tau, init_bac, init_antibiotic):
    r, m, n, v, a1, a2, k, a, Bdead = list(parameters_np)

    # initial values
    B = init_bac
    A = init_antibiotic
    t = 0

    # track number of bacteria over time
    bacteria_list = []
    antibiotic_list = []
    bacteria_list.append(B)
    antibiotic_list.append(A)

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
        antibiotic_list.append(A)
    
    return bacteria_list, antibiotic_list

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
#Parameter	Description	Units
#alpha	Antibiotic decay rate.	min−1
#r	Growth rate of rod-shaped bacteria.	min−1
#r_tilde 	Decay rate of nutrient due to consumption by rod-shaped bacteria.	OD−1 min−1
#gamma	Transition rate from rod to spherical shape.	min−1
#delta	Transition rate from spherical to rod shape.	min−1
#rho	Death rate of rod-shaped bacteria due to antibiotic.	min−1
#rho_deg 	Decay rate of antibiotic due to irreversible binding to bacteria.	μg ml−1 OD−1 min−1
#phi	Death rate of rod-shaped bacteria.	min−1
#psi	Death rate of spherical bacteria.	min−1
#T50	Antibiotic concentration required for half maximal killing effect	μg ml−1
#A50	Antibiotic concentration required for half maximal transition effect	μg ml−1

def simulate_mechanistic_markov(parameters_np, param_num, tot_time, tau, init_rod_OD, init_antibiotic):

    r, r_tilde, phi, alpha, gamma, delta, psi, rho, rho_deg, T_50, A_50 = list(parameters_np[:, param_num])
    

    # initial values
    H = init_rod_OD
    S = 0
    A = init_antibiotic
    N = 1

    t = 0
    # track number of bacteria over time
    rod_list = []
    sphere_list = []
    antibody_list = []
    nut_proportion = []
    rod_list.append(H)
    sphere_list.append(S)
    antibody_list.append(A)
    nut_proportion.append(N)

    while t < tot_time: # for all time
        # forward euler
        dH = r*N*H -(gamma*A/(T_50+A))*H + delta*S - (rho*A/(A_50+A))*H - phi*H
        dS = (gamma*A/(T_50+A))*H - delta*S - psi*S
        dA = -alpha*A - (rho_deg*A/(A_50+A))*H
        dN = -r_tilde*N*H

        H = H + tau * dH
        S = S + tau * dS
        A = A + tau * dA
        N = N + tau * dN

        t += tau
        rod_list.append(H)
        sphere_list.append(S)
        antibody_list.append(A)
        nut_proportion.append(N)
    
    return rod_list, sphere_list, antibody_list, nut_proportion


def mechanistic_markov():
    parameters_np = np.loadtxt(sys.argv[1])
    total_time = 100000
    tau = 1
    init_rod_OD = 0.17
    dosages = [3,4,10] #[0, 0.1, 0.2, 0.3, 0.4]
    param_num = 1 # 0 or 1, all dead or recover
    

    # for all initial antibody dosages
    for init_A in dosages:
        # average over num_runs
        rod_list, sphere_list, antibody_list, nut_proportion = simulate_mechanistic_markov(parameters_np, param_num, tot_time=total_time, tau=tau, init_rod_OD=init_rod_OD, init_antibiotic=init_A)

        total_OD_list = np.array(rod_list) + np.array(sphere_list) / 2
        total_bac_list = np.array(rod_list) + np.array(sphere_list)
        total_bac_list = total_bac_list * 10e5 / init_rod_OD

        fig, ax = plt.subplots(3,1, sharex=True, figsize=(7,6), gridspec_kw={'height_ratios': [3, 2, 2]})
        ax[0].set_title(f'Bacterial growth over time. Meropenem dosage: {init_A}μg')
        ax[0].plot(rod_list, label=f'rods')
        ax[0].plot(sphere_list, label=f'spheres')
        ax[0].plot(total_OD_list, label=f'combined OD')
        ax[0].legend(title='Cell type')
        ax[0].set_ylabel('Optical density')
        

        ax[1].plot(antibody_list, label=f'antibody μg/ml')
        ax[1].plot(nut_proportion, label=f'nutrient proportion')
        ax[1].legend()
        ax[1].set_ylabel('Concentration')
        
        ax[2].plot(total_bac_list, label=f'bacteria population')
        ax[2].legend()
        ax[2].set_ylabel('Count')
        ax[2].set_xlabel(f'Time (min)')
        plt.savefig(f'plots/mechanistic/mechanistic_markov_dos-{init_A}_param_{param_num}_long.png')
        plt.clf
    return


def simulate_new_markov(parameters_np, param_num, tot_time, tau, init_rod_OD, init_antibiotic, beta, init_res_prop):

    r, r_tilde, phi, alpha, gamma, delta, psi, rho, rho_deg, T_50, A_50= list(parameters_np[:, param_num])
    

    # initial values
    Hs = init_rod_OD * (1 - init_res_prop)
    Hr = init_rod_OD * init_res_prop
    Ss = 0
    Sr = 0
    A = init_antibiotic
    N = 1

    t = 0
    # track number of bacteria over time
    rods_list = []
    spheres_list = []
    rodr_list = []
    spherer_list = []
    antibody_list = []
    nut_proportion = []
    rods_list.append(Hs)
    spheres_list.append(Ss)
    rodr_list.append(Hr)
    spherer_list.append(Sr)
    antibody_list.append(A)
    nut_proportion.append(N)

    while t < tot_time: # for all time
        # forward euler
        dHs = r*N*Hs -(gamma*A/(T_50+A))*Hs + delta*Ss - (rho*A/(A_50+A))*Hs - phi*Hs - beta * Hs * A
        dHr = r*N*Hr -(gamma*A/(T_50+A))*Hr + delta*Sr - phi*Hr + beta * Hr * A
        dSs = (gamma*A/(T_50+A))*Hs - delta*Ss - psi*Ss
        dSr = (gamma*A/(T_50+A))*Hr - delta*Sr - psi*Sr
        dA = -alpha*A - (rho_deg*A/(A_50+A))*Hs
        dN = -r_tilde*N*Hs - r_tilde*N*Hr

        Hs = Hs + tau * dHs
        Hr = Hr + tau * dHr
        Ss = Ss + tau * dSs
        Sr = Sr + tau * dSr
        A = A + tau * dA
        N = N + tau * dN

        t += tau
        rods_list.append(Hs)
        spheres_list.append(Ss)
        rodr_list.append(Hr)
        spherer_list.append(Sr)
        antibody_list.append(A)
        nut_proportion.append(N)
    
    return rods_list, spheres_list,rodr_list,spherer_list, antibody_list, nut_proportion


def new_markov():
    parameters_np = np.loadtxt(sys.argv[1])
    total_time = 10000
    tau = 1
    beta = 10e-6
    init_res_prop = 0.1
    init_rod_OD = 0.17
    dosages = [2] #[0, 0.1, 0.2, 0.3, 0.4]
    param_num = 1 # 0 or 1, all dead or recover
    

    # for all initial antibody dosages
    for init_A in dosages:
        # average over num_runs
        rods_list, spheres_list,rodr_list,spherer_list, antibody_list, nut_proportion = simulate_new_markov(parameters_np, param_num, tot_time=total_time, tau=tau, init_rod_OD=init_rod_OD, init_antibiotic=init_A, beta=beta, init_res_prop=init_res_prop)

        total_OD_list = np.array(rods_list) + np.array(spheres_list) / 2 + np.array(rodr_list) + np.array(spherer_list) / 2
        total_bac_list = np.array(rods_list) + np.array(spheres_list) + np.array(rodr_list) + np.array(spherer_list)
        total_bac_list = total_bac_list * 10e5 / init_rod_OD

        total_resistant_list = (np.array(rodr_list) + np.array(spherer_list)) * 10e5
        total_sensitive_list = (np.array(rods_list) + np.array(spheres_list)) * 10e5

        total_rod_list = (np.array(rodr_list) + np.array(rods_list)) * 10e5
        total_sphere_list = (np.array(spherer_list) + np.array(spheres_list)) * 10e5

        fig, ax = plt.subplots(3,1, sharex=True, figsize=(7,6), gridspec_kw={'height_ratios': [3, 2, 2]})
        ax[0].set_title(f'Bacterial growth over time. Meropenem dosage: {init_A}μg')
        ax[0].plot(rods_list, label=f'susceptible rods')
        ax[0].plot(spheres_list, label=f'susceptible spheres')
        ax[0].plot(rodr_list, label=f'resistant rods')
        ax[0].plot(spherer_list, label=f'resistant spheres')
        ax[0].plot(total_OD_list, label=f'combined OD')
        ax[0].legend(title='Cell type')
        ax[0].set_ylabel('Optical density')
        

        ax[1].plot(antibody_list, label=f'antibody μg/ml')
        ax[1].plot(nut_proportion, label=f'nutrient proportion')
        ax[1].legend()
        ax[1].set_ylabel('Concentration')
        
        ax[2].plot(total_bac_list, label=f'bacteria population')
        ax[2].legend()
        ax[2].set_ylabel('Count')
        ax[2].set_xlabel(f'Time (min)')
        plt.savefig(f'plots/new/mechanistic_markov_dos-{init_A}_param_{param_num}.png')
        plt.clf

        fig, ax = plt.subplots(3,1, sharex=True, figsize=(7,6), gridspec_kw={'height_ratios': [3, 3, 2]})
        ax[0].set_title(f'Total sensitive/resistant and rod/sphere population. Meropenem dosage: {init_A}μg')
        ax[0].plot(total_sensitive_list, label=f'susceptible')
        ax[0].plot(total_resistant_list, label=f'resistant')
        ax[0].legend()
        ax[0].set_ylabel('Count')

        ax[1].plot(total_rod_list, label=f"rods")
        ax[1].plot(total_sphere_list, label=f"spheres")
        ax[1].legend()
        ax[1].set_ylabel('Count')

        ax[2].plot(antibody_list, label=f'antibody μg/ml')
        ax[2].plot(nut_proportion, label=f'nutrient proportion')
        ax[2].legend()
        ax[2].set_ylabel('Concentration')
        
        ax[2].set_xlabel(f'Time (min)')
        plt.savefig(f'plots/new/mechanistic_markov_dos-{init_A}_param_{param_num}_types.png')
        plt.clf

    return

def main():
    #simple_markov()
    #mechanistic_markov()
    new_markov()

if __name__ == "__main__":
    main()