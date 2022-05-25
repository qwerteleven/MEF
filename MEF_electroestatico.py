import numpy as np

def det_j(nodos):
    n_1, n_2, n_3 = nodos

    J = [
        [n_2[0] - n_1[0], n_3[0] - n_1[0]],
        [n_2[1] - n_1[1], n_3[1] - n_1[1]]
    ]

    det_J = (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]) 

    return det_J


def j_tras_inv(nodos):
    n_1, n_2, n_3 = nodos
    from numpy.linalg import inv

    a = np.array([
        [n_2[0] - n_1[0], n_3[0] - n_1[0]],
        [n_2[1] - n_1[1], n_3[1] - n_1[1]]
    ])

    a = inv(a)
    a = np.transpose(a)

    return a


def condiciones_dirichlet(element):

    for n_nodo in element['n_nodo']:
        G = 10e6
        K[n_nodo - 1, n_nodo - 1] = G
        F[n_nodo - 1] = G * element['condition_value']


def get_k(index_n, index_m, nodos, e = 1):

    determinante = det_j(nodos)
    jacobiano_inv = j_tras_inv(nodos)

    k = e * 1/2 * (jacobiano_inv @ gr[index_n]) @ (jacobiano_inv @ gr[index_m]) * determinante

    return k 


def get_f(index_n, nodos):
    determinante = det_j(nodos)

    f = gr[index_n] @ wi[index_n] * determinante * w[index_n]

    return f


def condiciones_neumann(element):
    import math
    eta = [
        (1 - (1 / math.sqrt(3))) / 2,
        (1 + (1 / math.sqrt(3))) / 2
    ]

    for index, n_nodo in enumerate(element['n_nodo']):
        N = 0
        n_1 = element['nodos'][0]
        n_2 = element['nodos'][1]

        refence = 0.5 * math.sqrt((n_2[0] - n_1[0]) ** 2  + (n_2[1] - n_1[1]) ** 2) * element['condition_value'][index]
        
        N += refence * (1 - eta[0]) 
        N += refence * eta[1]

        F[n_nodo - 1] = N

    return element


def ensamblaje_k_f(element):

    for index_n in range(len(element['n_nodo'])):
        for index_m in range(len(element['n_nodo'])):

            K[element['n_nodo'][index_n] - 1, element['n_nodo'][index_m] - 1] += get_k(index_n,
                                                                                            index_m,
                                                                                            element['nodos'],
                                                                                            element['condition_value'])

        F[element['n_nodo'][index_n] - 1] += get_f(index_n, element['nodos']) 


def get_system(elements):

    for element in elements:
        if element['condition'] == 'E':
            ensamblaje_k_f(element)

    for element in elements:

        if element['condition'] == 'dirichlet':
            condiciones_dirichlet(element)

    for element in elements:
        if element['condition'] == 'neumann':
            condiciones_neumann(element)

    
def main(nodes, conections, elements):

    get_system(elements)
    
    X2 = np.linalg.solve(K, F)


    xs = [value[0]  for key, value in nodes.items()]
    ys = [value[1]  for key, value in nodes.items()]
    
    zs = np.array([value[0] for value in X2])


    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    surf = ax.plot_trisurf(xs, ys, zs, cmap=plt.cm.coolwarm)
    ax.set_title('Solution')
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


def process_nodes(string_nodes):
    nodes = {}
    n_nodes_group_line = 1
    group_first_line = 2

    lines_nodes = string_nodes.strip().split('\n')

    # n_nodes = int(lines_nodes[0].stip().split(' ')[-1])

    while (group_first_line < len(lines_nodes)):
        n_nodes_group = int(lines_nodes[n_nodes_group_line].strip().split(' ')[-1])

        group_node = [int(line.strip()) for line in lines_nodes[group_first_line:group_first_line + n_nodes_group]]

        group_position_node = [
                            [float(value) for value in line.strip().split(' ')] 
                            for line in lines_nodes[group_first_line + n_nodes_group:group_first_line + n_nodes_group * 2]
                        ]

        group_first_line += n_nodes_group * 2 + 1
        n_nodes_group_line = group_first_line - 1

        # print(group_first_line, group_first_line + n_nodes_group)

        for node, position in zip(group_node, group_position_node):
            nodes[str(node)] = position

    return nodes


def process_conections(string_conections):
    
    conections = []
    n_nodes_group_line = 1

    lines_nodes = string_conections.strip().split('\n')

    # n_nodes = int(lines_nodes[0].stip().split(' ')[-1])

    while (n_nodes_group_line < len(lines_nodes)):
        n_nodes_group = int(lines_nodes[n_nodes_group_line].strip().split(' ')[-1])

        group_conections = [
                            [int(value) for value in line.strip().split(' ')][1:]
                            for line in lines_nodes[n_nodes_group_line + 1 : n_nodes_group_line + n_nodes_group]
                        ]

        n_nodes_group_line = n_nodes_group_line + n_nodes_group + 1

        # print(group_first_line, group_first_line + n_nodes_group)

        conections.append(group_conections)

    return conections


def create_elements(nodes, conections, conditions):
    elements = []

    for index_group, group_elements in enumerate(conections):
        for element_nodes in group_elements:

            element = {
                "n_nodo": element_nodes,
                "nodos": [nodes[str(node)] for node in element_nodes],
                "condition": conditions[index_group][0],
                "condition_value": conditions[index_group][1]
            }

            if conditions[index_group][0] == 'neumann':
                import math
                element["condition_value"] = [math.sin(node[1]) for node in element["nodos"]]

            elements.append(element)

    return elements


def read_msh(filename):
    import re

    f = open(filename, 'r')
    text = ''.join([line for line in f])
    f.close()
    
    string_nodes = re.findall("(?s)(?<=\$Nodes).*?(?=\$EndNodes)", text)
    string_conections = re.findall("(?s)(?<=\$Elements).*?(?=\$EndElements)", text)

    nodes = process_nodes(string_nodes[0])
    conections = process_conections(string_conections[0])

    elements = create_elements(nodes, conections, conditions)

    return nodes, conections, elements



conditions = [
        ['dirichlet', 0],

        ['dirichlet', 0],

        ['dirichlet', 0],
        
        ['neumann', 0],

        ['none', 0],
        ['dirichlet', 1],
        ['dirichlet', 1],

        ['E', 2],
        ['E', 1]
]


if __name__ == '__main__':
    filename = './mStrip.msh'
    nodes, conections, elements = read_msh(filename)

    gr = [
        np.array([-1, -1]),
        np.array([1,  0]),
        np.array([0,  1])
    ]

    w = [
        8/9,
        5/9,
        5/9
    ]

    wi = [
        [1/6, 1/6],
        [2/3, 1/6],
        [1/6, 2/3]
    ]

    K = np.zeros((len(nodes), len(nodes)))
    F = np.zeros((len(nodes), 1))
    main(nodes, conections, elements)

