from operator import itemgetter
import math
import tsplib95
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


class Cristophides(object):
    def run(self, data, name):
        self.generate_cities(data)
        MST, odd_vertices = self.mst()
        self.minimum_weight_matching(MST, self.edges, odd_vertices)
        euler_tour = self.euler_tour(MST)
        path = [euler_tour[0]]
        visited_vertices = set(path)
        distance = 0
        for i in euler_tour[1:]:
            if i not in visited_vertices:
                distance += self.edges[self.data * path[-1] + i][2]
                path.append(i)
                visited_vertices.add(i)
        path.append(path[0])
        # distance += self.edges[self.data * path[-1] + path[0]][2]
        self.path=path
        self.data=data
        print(path)
        print(distance)
        self.plot(name=name)

    def generate_cities(self, filename):
        problem = tsplib95.load(filename)
        G = problem.get_graph()
        nodes = list(problem.get_nodes())
        outer_array = []
        for i, o in enumerate(nodes):
            for j, v in enumerate(nodes):
                if i == j:
                    outer_array.append([o - 1, v - 1, 0])
                else:
                    outer_array.append([o - 1, v - 1, G.edges[o, v]['weight']])
        self.edges = outer_array
        self.data = len(nodes)
        self.graph=problem

    def mst(self):
        weighted_edges = sorted(self.edges, key=itemgetter(2))
        parent = []
        rank = []
        oddness = []
        final_mst = []

        def find_parent(parent, x):
            if parent[x] == x:
                return x
            return find_parent(parent, parent[x])

        def union(parent, rank, x, y):
            parent_x = find_parent(parent, x)
            parent_y = find_parent(parent, y)
            if rank[parent_x] > rank[parent_y]:
                parent[parent_y] = parent_x
            elif rank[parent_x] < rank[parent_y]:
                parent[parent_x] = parent_y
            else:
                parent[parent_x] = parent_y
                rank[parent_y] += 1

        for i in range(self.data):
            parent.append(i)
            rank.append(0)
            oddness.append(0)
        for i in weighted_edges:
            x_parent = find_parent(parent, i[0])
            y_parent = find_parent(parent, i[1])
            if x_parent != y_parent:
                final_mst.append(i)
                union(parent, rank, i[0], i[1])
                oddness[i[0]] += 1
                oddness[i[1]] += 1
        odd_vertices = [x for x in range(len(oddness)) if oddness[x] % 2 == 1]
        return final_mst, odd_vertices

    def minimum_weight_matching(self, MST, graph, odd_vertices):
        while(odd_vertices):
            odd_vertex = odd_vertices.pop()
            length = 999999999999
            for v in odd_vertices:
                if graph[self.data * odd_vertex + v][2] < length:
                    closest = v
                    length = graph[self.data * odd_vertex + v][2]
            MST.append([odd_vertex, closest, length])
            odd_vertices.remove(closest)

    def euler_tour(self, MatchedMSTree):
        # find neigbours
        neighbours = {}
        for edge in MatchedMSTree:
            if edge[0] not in neighbours:
                neighbours[edge[0]] = []

            if edge[1] not in neighbours:
                neighbours[edge[1]] = []

            neighbours[edge[0]].append(edge[1])
            neighbours[edge[1]].append(edge[0])

        tour = []

        def recursively_calculate_euler_tour(neighbours, v):
            for i in neighbours[v]:
                neighbours[v].remove(i)
                neighbours[i].remove(v)
                recursively_calculate_euler_tour(neighbours, i)
            tour.append(v)
        recursively_calculate_euler_tour(neighbours, MatchedMSTree[0][0])
        return tour

    def plot(self, line_width=1, point_radius=math.sqrt(2.0), annotation_size=8, dpi=120, save=True, name=None):
        x = [self.graph.display_data[i+1][0] for i in self.path-1]
        x.append(x[0])
        y = [self.graph.display_data[i+1][1] for i in self.path-1]
        y.append(y[0])
        plt.plot(x, y, linewidth=line_width)
        plt.scatter(x, y, s=math.pi * (point_radius ** 2.0))
        plt.title('Christofides '+self.data)
        j=0
        for i in self.path:
            plt.annotate(j, self.graph.display_data[i+1], size=annotation_size)
            j+=1
        if save:
            if name is None:
                name = '{0}.png'.format('Christofides '+self.data)
            plt.savefig(name, dpi=dpi)
        plt.show()
        plt.gcf().clear()


sol = Cristophides()
sol.run('ALL_tsp/ali535.tsp', 'ali535')
