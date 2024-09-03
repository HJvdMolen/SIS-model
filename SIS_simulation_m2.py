''' Run a single, bottom-up simulation of the heterogeneous SIS model. 
The SIS model used is heterogeneous, with (m =) 2 subpopulations. 

Graph implementation based on: https://stackabuse.com/courses/
graphs-in-python-theory-and-implementation/lessons/representing-graphs-in-code/

*** Formulas *** 
rnd = random sample from U[0..1]

Chance:
node is infected by 1 infected contact  = beta * 1
node is infected by k infected contacts = beta * k
node recovers from infection            = delta
node remains infected                   = 1 - delta

If  node = 0    # Susceptible, chance of (re)infection =
    node = 1 * (rnd < beta * count_infected_contacts)
Else            # Infected, chance to remain infected =
    node = 1 * (rnd > delta)

Author: Henk-Jan van der Molen, 2024-08-29'''

def plot_function(function, title, xlabel, ylabel, v_8):
    import matplotlib.pyplot as plt

    plt.plot(function)
    plt.title(title)
    plt.plot([v_8 for _ in range(10 + len(function))], linestyle='dashed', linewidth=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    return


def openFile(folder, textfile, mode):
    import os, sys
    '''
    file mode:  r = read only (default), 
                w = write only (erases existing file), 
                a = append, 
                r+ = read + write
    '''
    file = os.path.join(folder, textfile)
    try:
        handle = open(file, mode)
        print(f"Opened file {file} in mode {mode}")
    except IOError:
        print(f"Cannot open {file} in mode {mode}!")
        raw_input()
        sys.exit()
    return handle


def calc_beta_delta(ri, k): # calculates β, δ  as big as possible
    assert k % 4 == 0,      f"Variable k must be a multiple of 4, got {k}"
    assert k >= 4,          f"Variable k must be >= 4, got {k}"
    assert ri >= 0,         f"R₁ must be >= 0, got {ri}"

    delta = 1 / 2
    
    while True:
        beta = ri * delta / k
        if beta * k < 1 : break
        delta /= 2
                
    assert 0 < delta <= 1,      f"0 < δ <= 1 expected, got: {delta}"
    assert beta * k <= 1,       f"β.k.v <= 1 expected, got: {beta} * {k}"
    
    #print(f"R₁ = {ri}; k = {k}; beta = {beta:.3f}; delta = {delta:.3f}")
    return beta, delta
    

class Graph:
    # Constructor
    def __init__(self, network_size, n1, no_of_steps, r0, D):
        self.set_vars(network_size, n1, no_of_steps, r0, D)
        return
        

    def set_vars(self, network_size, n1, no_of_steps, r0, D):
        self.network_size = network_size
        self.n1 = n1
        self.n2 = 1 - n1
        self.no_of_steps = no_of_steps
        self.r0 = r0
        self.D  = D
        
        # Calculated later...
        self.v_8 = self.D = 0
        self.mu = self.sd = 0

        # Create simulation + average simulation list
        self.sim_tot = [0 for j in range(self.no_of_steps)]
        self.sim_avg = [0 for j in range(self.no_of_steps)]
        return
        

    def calc_steady_state(self, p1, p2):    # calculate steady state of v
        from math import sqrt
        
        a           = p1.r0 * p2.r0
        if a != 0:
            b       = p1.r0 + p2.r0 - a
            c       = 1 - self.r0
            d       = b * b - 4 * a * c
            self.v_8 = (sqrt(d) - b) / (2 * a)
        elif p1.r0 == 0:
            self.v_8 = self.n2 - 1 / p2.r0
        elif p2.r0 == 0:
            self.v_8 = self.n1 - 1 / p1.r0

        # print(f"*** SIS simulation with {self.no_of_steps} time steps ***")
        # print(f"#nodes = {self.network_size}, n₁ = {self.n1}, k = {p1.k}, D = {self.D:.2f}") 
        # print(f"R₀, R₁, R₂ = {self.r0}, {p1.r0}, {p2.r0}; v(∞) = {self.v_8:.5f}")
        # print()
        return


    # Execute complete simulation
    def run_simulation(self, p1, p2):
        from statistics import mean, stdev
        
        p1.init_simulation()
        p2.init_simulation()
        self.sim_tot = [0 for j in range(self.no_of_steps)]

        for j in range(1, self.no_of_steps):
            p1.simulation_step(p2, j)
            p2.simulation_step(p1, j)

            avg_v = self.n1 * p1.sim[j] + self.n2 * p2.sim[j]
            self.sim_tot[j]  = avg_v
            self.sim_avg[j] += avg_v
                        
        # Calculate (μ, σ) from steady state part = right half of list[]
        self.mu =  mean(self.sim_tot[self.no_of_steps//2:])
        self.sd = stdev(self.sim_tot[self.no_of_steps//2:]) if self.sim_tot[-1] > 0 else 0
        return


class SubPopulation:
    def __init__(self, num_of_nodes, beta, k, delta, no_of_steps, directed):
        self.set_vars(num_of_nodes, beta, k, delta, no_of_steps, directed)

        # Initialize & fill the adjacency list
        self.adj_list = {node: set() for node in range(num_of_nodes)}

        # adj_list is also used for contacts between SubPopulations = /2
        if self.directed:   # directed graph
            contacts_out = self.k // 2
        else:               # undirected = include "mirror" edges
            contacts_out = self.k // 4

        for key in self.adj_list.keys():
            for j in range(contacts_out):
                # Link node(n) to node(n+1, n+2, .. , n+contacts_out)
                self.add_edge(key, (key + 1 + j) % self.num_of_nodes, directed)
        return


    def set_vars(self, num_of_nodes, beta, k, delta, no_of_steps, directed):
        self.num_of_nodes = num_of_nodes
        self.no_of_steps = no_of_steps
        self.directed = directed

        self.beta = beta
        self.k = k
        self.delta = delta
        self.r0 = beta * k / delta
        
        # Create the nodes list
        self.node = [0 for j in range(self.num_of_nodes)]

        # Create Temp copy of self.node[]
        self.temp = [0 for j in range(self.num_of_nodes)]

        # Create simulation list for SubPopulation
        self.sim = [0 for j in range(self.no_of_steps)]
        return 


    def add_edge(self, node1, node2, directed): # Add edge node1 -> node2
        assert node1 != node2, f"node {node1} cannot connect to itself"
        
        self.adj_list[node1].add(node2)
        if self.directed == False:              # also add "mirror" edge
            self.adj_list[node2].add(node1)
        return


    # Export representation of the graph to text file
    def export_graph(self):
        path = "/home/henk-jan/programs/python/SIS-simulation/"
        textfile = "graph.txt"
        export = openFile(path, textfile, 'w')
        print(f"Export the graph to file: {path + textfile}")
        
        for key in self.adj_list.keys():
            line = f"Connections of node {key} : "
            for node in self.adj_list[key]:
                line += str(node) + ", "
                
            line += "\n"
            export.write(line)

        export.close()
        print("... Export file closed.")
        return


    def init_simulation(self):
        import random
        
        # Set all nodes to be Susceptible (== 0)
        for j in range(self.num_of_nodes):
            self.node[j] = 0
            self.temp[j] = 0
        
        if False:    # if True, start with 100% infected
            self.nodes_infected = self.num_of_nodes
            self.node = [1 for j in range(self.num_of_nodes)]
        else:        # Randomly infect 1% of susceptible nodes (0 => 1)
            self.nodes_infected = self.num_of_nodes // 100
            for j in range(self.nodes_infected):
                while True:
                    p = random.randrange(0, self.num_of_nodes)
                    if self.node[p] == 0 : break
                self.node[p] = 1    # node is infected
            
        # Fill in first value in simulation list 
        self.sim = [0 for j in range(self.no_of_steps)]
        self.sim[0] = self.nodes_infected / self.num_of_nodes
        
        return


    # #infected contacts (node[]==1); p1 + p2 have identical adjacency list
    def count_infected_contacts(self, other, key):
        cic = 0
        for nd in self.adj_list[key]:
            cic +=  self.node[nd]   # If node[nd] = 1 -> Infected
            cic += other.node[nd]
        return cic


    # Execute single time step in the simulation for all nodes
    def simulation_step(self, other, step_no):
        import random
        
        self.nodes_infected = 0
        for j in range(self.num_of_nodes):
            rnd = random.uniform(0,1)
            
            if self.node[j] == 0:    # Susceptible node, can be infected 
                contacts = self.count_infected_contacts(other, j)
                assert self.beta*contacts <= 1, f"β.(k.v) > 1, with {self.beta}*{contacts}"
                
                self.temp[j] = 1 * (rnd < self.beta * contacts)
            else:                    # Infected node, can recover 
                self.temp[j] = 1 * (rnd > self.delta)
            self.nodes_infected += self.temp[j]

        self.sim[step_no] = self.nodes_infected / self.num_of_nodes
        self.node = self.temp[:]

        return
 

# Main program
if __name__ == '__main__':
    network_size  = 5000             # no. of nodes in the network
    k             = 16               # no. of contacts for each node
    n1            = 0.5              # fraction of SubPopulation 1
    n2            = 1 - n1           # fraction of SubPopulation 2
    directed      = False            # False = bidirectional network
    assert n1 == n2 == 0.5, f"n₁, n₂ must be 0.5, got n1 = {n1}, n2 = {n2}"
    
    r0            = 5.0              # basic reproduction number R₀
    assert r0 > 1, f"R₀ > 1 expected, got: {r0}"

    D             = 0.75             # Diversity index, [0 = min, 1 = max]
    assert 0 <= D <= 1, f"0 <= Diversity <= 1, got {D}"
    
    # SubPopulation #1 & #2
    num_of_nodes1 = int(network_size * n1)
    num_of_nodes2 = network_size - num_of_nodes1

    r1          = r0 * (1 + D / (4 * n1 * n2))
    r2          = (r0 - n1 * r1) / n2
    
    # β: infection chance from 1 infected contact | δ: recovery from infection
    beta1, delta1 = calc_beta_delta(r1, k)
    beta2, delta2 = calc_beta_delta(r2, k)

    no_of_steps   = 2000             # no. of time steps in 1 simulation

    g             = Graph(network_size, n1, no_of_steps, r0, D)
    p1            = SubPopulation(num_of_nodes1, beta1, k, delta1, no_of_steps, directed)
    p2            = SubPopulation(num_of_nodes2, beta2, k, delta2, no_of_steps, directed)
    # p1.export_graph()
    
    g.calc_steady_state(p1, p2)       # calculates v(∞)
    g.run_simulation(p1, p2)          # execute all time steps 4all nodes
    
    title  = f"SIS simulation: R₀, R₁, R₂ = {g.r0:.2f}, {p1.r0:.2f}, {p2.r0:.2f}; v(∞) = {g.v_8:.3f}"
    xlabel = f"#time steps: {no_of_steps} with {network_size} nodes, k = {k}, directed network: {directed}"
    ylabel = f"μ(v) = {g.mu:.5f}, σ(v) ={g.sd:.5f}"
    plot_function(g.sim_avg,  title, xlabel, ylabel, g.v_8)
