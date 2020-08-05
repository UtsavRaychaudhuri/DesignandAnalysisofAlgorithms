from matplotlib import pyplot as plt
from random import uniform
import tsplib95
import math
import matplotlib
matplotlib.use('Agg')


class TSP_ACS(object):

    def __init__(
            self,
            steps,
            rho=0.1,
            alpha=1,
            beta=3,
            elite_weight=1,
            pheromone_deposit_weight=1.0,
            min_scaling_factor=0.001):
        self.steps = int(steps)
        self.depot_assigned_to_ant = 0
        self.best_distance = 99999999999999
        self.best_route = None
        self.rho = rho
        self.alpha = alpha
        self.beta = beta
        self.elite_weight = elite_weight
        self.pheromone_deposit_weight = pheromone_deposit_weight
        self.min_scaling_factor = min_scaling_factor

    def load_file_and_create_graph(self, filename):
        problem = tsplib95.load(filename)
        G = problem.get_graph()
        nodes = list(problem.get_nodes())
        outer_array = []
        for i, o in enumerate(nodes):
            inner_array = []
            for j, v in enumerate(nodes):
                if i == j:
                    inner_array.append([0, 1])
                elif i > j:
                    inner_array.append(outer_array[j][i])
                else:
                    inner_array.append([G.edges[o, v]['weight'], 1])
            outer_array.append(inner_array)
        self.edges = outer_array
        self.cities = self.ants = len(nodes)
        self.nodes = nodes
        self.graph = problem

    def ant_select_a_city_to_explore(self):
        roulette_wheel = 0
        tabu = [city for city in range(self.cities) if city not in self.route]
        heuristic_total_weight = 0
        for i in tabu:
            heuristic_total_weight += self.edges[self.route[-1]][i][0]
        for i in tabu:
            if self.edges[self.route[-1]
                          ][i][1] == 0 or self.edges[self.route[-1]][i][0] == 0:
                continue
            roulette_wheel += pow(self.edges[self.route[-1]][i][1], self.alpha) * pow(
                (heuristic_total_weight / self.edges[self.route[-1]][i][0]), self.beta)
        random_val = uniform(0, roulette_wheel)
        roulette_wheel_position = 0
        for i in tabu:
            if self.edges[self.route[-1]
                          ][i][1] == 0 or self.edges[self.route[-1]][i][0] == 0:
                continue
            roulette_wheel_position += pow(self.edges[self.route[-1]][i][1], self.alpha) * pow(
                (heuristic_total_weight / self.edges[self.route[-1]][i][0]), self.beta)
            if roulette_wheel_position >= random_val:
                return i

    def ant_find_route(self):
        self.route = [self.depot_assigned_to_ant]
        self.depot_assigned_to_ant += 1
        while len(self.route) < self.cities:
            self.route.append(self.ant_select_a_city_to_explore())
        return self.route

    def add_pheromone(self, tour, distance, weight=1):
        pheromone_to_add = 1 / distance
        for i in range(self.cities - 2):
            self.edges[self.route[i]][self.route[i + 1]
                                      ][1] += weight * pheromone_to_add

    def calculate_route_distance(self):
        self.distance = 0
        for i in range(self.cities - 2):
            self.distance += self.edges[self.route[i]][self.route[i + 1]][0]
        return self.distance

    def acs(self):
        for step in range(self.steps):
            for ant in range(self.ants):
                self.add_pheromone(
                    self.ant_find_route(),
                    self.calculate_route_distance())
                if self.distance < self.best_distance:
                    self.best_distance = self.distance
                    self.best_route = self.route
            self.depot_assigned_to_ant = 0
            for i in range(self.cities - 1):
                for j in range(i + 1, self.cities - 1):
                    self.edges[i][j][1] *= 1 - self.rho

    def elitist(self):
        for _ in range(self.steps):
            for ant in range(self.ants):
                self.add_pheromone(
                    self.ant_find_route(),
                    self.calculate_route_distance())
                if self.distance < self.best_distance:
                    self.best_distance = self.distance
                    self.best_route = self.route
            self.add_pheromone(
                self.best_route,
                self.best_distance,
                weight=self.elite_weight)
            self.depot_assigned_to_ant = 0
            for i in range(self.cities - 1):
                for j in range(i + 1, self.cities):
                    pheromone = self.edges[i][j][1]
                    self.edges[i][j][1] *= 1 - self.rho

    def max_min(self):
        for step in range(self.steps):
            iteration_best_tour = None
            iteration_best_distance = float("inf")
            for ant in range(self.ants):
                self.ant_find_route()
                if self.calculate_route_distance() < iteration_best_distance:
                    iteration_best_tour = self.route
                    iteration_best_distance = self.distance
            self.depot_assigned_to_ant = 0
            if float(step + 1) / float(self.steps) <= 0.75:
                self.add_pheromone(
                    iteration_best_tour,
                    iteration_best_distance)
                max_pheromone = self.pheromone_deposit_weight / iteration_best_distance
            else:
                if iteration_best_distance < self.best_distance:
                    self.best_route = iteration_best_tour
                    self.best_distance = iteration_best_distance
                self.add_pheromone(self.best_route, self.best_distance)
                max_pheromone = self.pheromone_deposit_weight / self.best_distance
            min_pheromone = max_pheromone * self.min_scaling_factor
            for i in range(self.cities - 1):
                for j in range(i + 1, self.cities - 1):
                    self.edges[i][j][1] *= (1.0 - self.rho)
                    if self.edges[i][j][1] > max_pheromone:
                        self.edges[i][j][1] = max_pheromone
                    elif self.edges[i][j][1] < min_pheromone:
                        self.edges[i][j][1] = min_pheromone

    def run(self, mode):
        self.mode = mode
        print('Started : {0}'.format(self.mode))
        if self.mode == '1':
            self.mode = "ACS"
            self.acs()
        elif self.mode == '2':
            self.mode = "Elitist"
            self.elitist()
        else:
            self.mode = "Max-Min"
            self.max_min()
        print('Ended : {0}'.format(self.mode))
        print(self.best_route)
        print(
            'Total distance travelled to complete the tour : {0}\n'.format(
                round(
                    self.best_distance,
                    2)))
        self.plot()

    def plot(self, line_width=1, point_radius=math.sqrt(0.1), annotation_size=1, dpi=120, save=True, name=None):
        x = [self.graph.node_coords[i+1][0] for i in self.best_route]
        x.append(x[0])
        y = [self.graph.node_coords[i+1][1] for i in self.best_route]
        y.append(y[0])
        plt.plot(x, y, linewidth=line_width)
        plt.scatter(x, y, s=math.pi * (point_radius ** 2.0))
        plt.title(self.mode)
        j = 0
        for i in self.best_route:
            plt.annotate(j, self.graph.node_coords[i+1], size=annotation_size)
            j += 1
        if save:
            if name is None:
                name = '{0}.png'.format(self.mode)
            plt.savefig(name, dpi=dpi)
        plt.show()
        plt.gcf().clear()


path = input("Please enter the tsp file and the path to that file")
print(path)
print("1. ACS Without Optimization")
print("2. ACS with Elitist Optimization")
print("3.ACS with MinMax Optimization")
choice = input("Your choice- 1/2/3?")
steps = input("The number of steps you want to run this for")
sol = TSP_ACS(steps)

sol.load_file_and_create_graph(path)
sol.run(choice)
