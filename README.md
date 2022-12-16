# AntibioticResistanceModeling

    python simulate.py \
        params.txt

# Parameters

Model 3 and 4 parameters are saved in mechanistic_markov_params.txt

| Parameter |	Description	| Units |
| ------- | -------- | --------- |
| alpha |	Antibiotic decay rate.	| min−1 |
| r |	Growth rate of rod-shaped bacteria.	| min−1 |
| r_tilde | 	Decay rate of nutrient due to consumption by rod-shaped bacteria. | OD−1 min−1 |
| gamma |	Transition rate from rod to spherical shape.	| min−1 |
| delta |	Transition rate from spherical to rod shape.	| min−1 |
| rho |	Death rate of rod-shaped bacteria due to antibiotic.	| min−1 |
| rho_deg | 	Decay rate of antibiotic due to irreversible binding to bacteria. | μg ml−1 OD−1 min−1 |
| phi |	Death rate of rod-shaped bacteria.	| min−1 |
| psi |	Death rate of spherical bacteria.	| min−1 |
| T50 |	Antibiotic concentration required for half maximal killing effect | μg ml−1 |
| A50 |	Antibiotic concentration required for half maximal transition effect | μg ml−1 |
