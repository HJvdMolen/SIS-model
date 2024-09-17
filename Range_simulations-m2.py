"""
Run a range of simulations, exporting μ and σ  over no_of_runs simulations.
E.g. run simulations for {R₀ = 3, 4, .. , 9} * {D = 0, 0.25, 0.5, 0.75, 1}

(C) Henk-Jan van der Molen, 2024-09-28
"""

# Main program
from math import sqrt
import time
import SIS_simulation_m2 as sis

start_time = time.time()

network_size = 5000  # no. of nodes in the network
k = 32  # no. of contacts for each node
n1 = 0.5  # fraction of SubPopulation 1
n2 = 1 - n1  # fraction of SubPopulation 2
assert n1 == n2 == 0.5, f"n₁, n₂ must be 0.5, got n₁ = {n1}, n₂ = {n2}"

directed = False  # False = bidirectional network
no_of_steps = 2000  # no. of time steps in 1 simulation
no_of_runs = 50  # no. of simulations to run

# (network_size, n1, no_of_steps,r0, D)
g = sis.Graph(network_size, n1, no_of_steps, 0, 0)

# SubPopulations #1, #2
num_of_nodes1 = int(network_size * n1)
num_of_nodes2 = network_size - num_of_nodes1
#                                   num_of_nodes,beta,k,delta,no_of_steps,directed)
p1 = sis.SubPopulation(num_of_nodes1, 0, k, 1, no_of_steps, directed)
p2 = sis.SubPopulation(num_of_nodes2, 0, k, 1, no_of_steps, directed)

print("R₀; D; R₁; R₂; v(∞); μ(v); σ(v); #steps; #runs; #nodes; k; network=directed")

for r0 in range(4, 5):  # basic reproduction number = β.k / δ
    g.r0 = r0

    for x in range(5):
        g.D = x / 4  # Diversity [0, 0.25, 0.5, 0.75, 1=max]

        p1.r0 = r0 * (1 + g.D / (4 * n1 * n2))
        p2.r0 = (r0 - n1 * p1.r0) / n2
        g.calc_steady_state(p1, p2)  # calculates v(∞)

        # β: infection chance from 1 infected contact | δ: recovery chance
        p1.beta, p1.delta = sis.calc_beta_delta(p1.r0, k)
        p2.beta, p2.delta = sis.calc_beta_delta(p2.r0, k)

        total = total2 = 0
        for j in range(no_of_runs):
            g.run_simulation(p1, p2)
            total += g.mu
            total2 += g.mu**2

        mu = total / no_of_runs
        sd = sqrt((total2 - total**2 / no_of_runs) / no_of_runs)

        print(
            f"{r0:.2f}; {g.D:.2f}; {p1.r0:.2f}; {p2.r0:.2f}; {g.v_8:.5f}; {mu:.5f}; {sd:.5f}; {no_of_steps}; {no_of_runs}; {network_size}; {k}; {directed}"
        )


run_time = round(time.time() - start_time, 3)
print()
print(f"*** calculation time = {run_time} sec.")
