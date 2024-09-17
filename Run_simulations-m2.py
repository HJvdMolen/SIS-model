"""
Run multiple simulations, filling a list sim [] with v(∞) values.
Using each v(∞) value, calculate the μ and σ  over all simulations.

(C) Henk-Jan van der Molen, 2024-08-28
"""

# Main program
from math import sqrt
from statistics import mean, stdev
import time
import SIS_simulation_m2 as sis

start_time = time.time()

network_size = 5000  # no. of nodes in the network
k = 32  # no. of contacts for each node
n1 = 0.5  # fraction of SubPopulation 1
n2 = 1 - n1  # fraction of SubPopulation 2
directed = False  # False = bidirectional network
assert n1 == n2 == 0.5, f"n₁, n₂ must be 0.5, got n1 = {n1}, n2 = {n2}"

r0 = 4.0  # basic reproduction number = β.k / δ
assert r0 > 1, f"R₀ must be > 1, got {r0}"

D = 1.00  # Diversity index, [0 = min, 1 = max]
assert 0 <= D <= 1, f"0 <= Diversity <= 1, got {D}"

no_of_steps = 2000  # no. of time steps in 1 simulation
no_of_runs = 50  # no. of simulations to run

g = sis.Graph(network_size, n1, no_of_steps, r0, D)

# SubPopulation #1, #2
num_of_nodes1 = int(network_size * n1)
num_of_nodes2 = network_size - num_of_nodes1
r1 = r0 * (1 + D / (4 * n1 * n2))
r2 = (r0 - n1 * r1) / n2

# β: infection chance from 1 infected contact | δ: recovery from infection
beta1, delta1 = sis.calc_beta_delta(r1, k)
beta2, delta2 = sis.calc_beta_delta(r2, k)

p1 = sis.SubPopulation(num_of_nodes1, beta1, k, delta1, no_of_steps, directed)
p2 = sis.SubPopulation(num_of_nodes2, beta2, k, delta2, no_of_steps, directed)
# p1.export_graph()

g.calc_steady_state(p1, p2)  # calculates v(∞)
print(f"*** #runs = {no_of_runs} ***")
print()
print("Run# / v(∞): μ, σ")

total = total2 = 0
for j in range(no_of_runs):
    g.run_simulation(p1, p2)
    print(f"  {j:0=2d}; {g.mu:.5f}; {g.sd:.5f}")

    total += g.mu
    total2 += g.mu**2

mu = total / no_of_runs
sd = sqrt((total2 - total**2 / no_of_runs) / no_of_runs)

run_time = round(time.time() - start_time, 3)

print()
print("R₀; D; R₁; R₂; v(∞); μ(v); σ(v); #steps; #runs; #nodes; k; network=directed")
print(
    f"{r0}; {D}; {r1}; {r2}; {g.v_8:.5f}; {mu:.5f}; {sd:.5f}; {no_of_steps}; {no_of_runs}; {network_size}; {k}; {directed}"
)
print(f"*** calculation time = {run_time} sec.")

for j in range(no_of_steps):
    g.sim_avg[j] /= no_of_runs

# Calculate (μ, σ) from steady state part = right half of list[]
mu = mean(g.sim_avg[no_of_steps // 2 :])
sd = stdev(g.sim_avg[no_of_steps // 2 :])

title = f"SIS simulation: R₀, R₁, R₂ = {g.r0:.2f}, {p1.r0:.2f}, {p2.r0:.2f}; v(∞) = {g.v_8:.3f}"
xlabel = f"Time steps: {no_of_steps} with {network_size} nodes; network directed = {directed}"
ylabel = f"AVG of {no_of_runs} runs: μ(v) = {mu:.5f}, σ(v) ={sd:.5f}"
sis.plot_function(g.sim_avg, title, xlabel, ylabel, g.v_8)
