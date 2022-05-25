import numpy as np


def get_ae(n_1, n_2, n_3):
    k = 0.5 * ((n_2[0] - n_1[0]) * (n_3[1] - n_1[1]) - (n_3[0] - n_1[0]) * (n_2[1] - n_1[1]))
    return k

def get_eta(n_1, n_2, n_3, x, y):
    return (1 / (2 * get_ae(n_1, n_2, n_3))) * ((x - n_1[0]) * (n_3[1] - n_1[1]) - (n_3[0] - n_1[0]) * (y - n_1[1]))

def get_nu(n_1, n_2, n_3, x, y):
    return (1 / (2 * get_ae(n_1, n_2, n_3))) * ((n_2[0] - n_1[0]) * (y - n_1[1]) - (x - n_1[0]) * (n_2[1] - n_1[1]))


def u_h_solucion(x, y, solucion):
    u = (solucion[0] * (1 - get_eta(x, y) - get_nu(x, y)) + 
         solucion[1] * get_eta(x, y) + 
         solucion[2] * get_nu(x, y))

    return u

def get_chi(nodo):
    return nodo / get_ae(nodo[0], nodo[1])


def get_funtion_form(element):

    if len(element['nodos']) == 2:
        n_1, n_2 = element['nodos']
        n_3 = [0, 0]

        base_1 = (1 - get_eta(n_1, n_2, n_3, n_1[0], n_1[1]))
        base_2 = get_eta(n_1, n_2, n_3, n_2[0], n_2[1])
        element['forma'] = [base_1, base_2]

    if len(element['nodos']) == 3:
        n_1, n_2, n_3 = element['nodos']

        base_1 = (1 - get_eta(n_1, n_2, n_3, n_1[0], n_1[1]) - get_nu(n_1, n_2, n_3, n_1[0], n_1[1]))
        base_2 = get_eta(n_1, n_2, n_3, n_2[0], n_2[1])
        base_3 = get_nu(n_1, n_2, n_3, n_1[0], n_1[1])

        element['forma'] = [base_1, base_2, base_3]

    return element


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


def condiciones_dirichlet(n_nodos):

    for n_nodo, value in n_nodos:
        G = 10e6
        K[n_nodo, n_nodo] = G
        F[n_nodo] = G * value


def get_k(index_n, index_m, nodos, e = 1):
    
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

    k = 0

    determinante = det_j(nodos)
    jacobiano_inv = j_tras_inv(nodos)

    for index, nodo in enumerate(nodos):
        k = e * (((jacobiano_inv @ gr[index_n]) * (jacobiano_inv @ gr[index_m]) * determinante) @ wi[index]) * w[index]

    # k = e * 1/2 * (jacobiano_inv @ gr[index_n]) @ (jacobiano_inv @ gr[index_m]) * determinante
    return k 


def get_f(index_n, nodos):

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

    determinante = det_j(nodos)
    f = 0

    for _, nodo in enumerate(nodos):
        f += gr[index_n] @ wi[index_n] * determinante * w[index_n]

    return f


def f_neumann(j, ro, angular, forma):
    v = (- forma) / (j * angular * ro)
    return v


def jacobiano(n_1, n_2, n_3):
    return 2 * get_ae(n_1, n_2, n_3)


def ensamblaje_k_f(element):

    for index_n in range(len(element['n_nodo'])):
        for index_m in range(len(element['n_nodo'])):
            if element['condition'][index_m] == 'E':
                K[element['n_nodo'][index_n] - 1, element['n_nodo'][index_m] - 1] += get_k(index_n,
                                                                                            index_m,
                                                                                            element['nodos'],
                                                                                            element['condition_value'][index_m])
            else:
                K[element['n_nodo'][index_n] - 1, element['n_nodo'][index_m] - 1] += get_k(index_n, index_m, element['nodos'])

        F[element['n_nodo'][index_n] - 1] += get_f(index_n, element['nodos']) 


def get_sistem(elements, n_nodos_dir):

    for element in elements:
        ensamblaje_k_f(element)

    condiciones_dirichlet(n_nodos_dir)

    
def main():
    elements = get_elements()
    n_nodos_dir = []

    for element in elements:
        for index in range(len(element['n_nodo'])):
            if element['condition'][index] == 'dirichlet':
                n_nodos_dir.append((element['n_nodo'][index], element['condition_value'][index]))

    n_nodos_dir = list(set(n_nodos_dir))

    get_sistem(elements, n_nodos_dir)
    
    X2 = np.linalg.solve(K, F)

    
    coordinate_nodos = get_coor_nodes()

    xs = [value[0]  for key, value in coordinate_nodos.items()]
    ys = [value[1]  for key, value in coordinate_nodos.items()]
    
    zs = np.array([value[0] for value in X2])


    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    surf = ax.plot_trisurf(xs, ys, zs, cmap=plt.cm.coolwarm)
    ax.set_title('Solution')
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()



def get_conection_nodes():

    conection_nodos = {}

    f = open('./elec_elements.txt', 'r')

    while True:

        metadata_line = f.readline().rstrip("\n")
        if metadata_line == '':
            break

        n_artifacts = int(metadata_line.split(' ')[-1])

        for i in range(n_artifacts):
            conections = f.readline().strip().split(' ')
            if len(conections) == 4:
                conection_nodos[conections[0]] = [int(coor) for coor in conections[1:]]

    f.close()

    return conection_nodos


def get_coor_nodes():
    coordinate_nodos = {}

    f = open('./elec_nodes.txt', 'r')

    while True:

        metadata_line = f.readline().rstrip("\n")
        if metadata_line == '':
            break

        n_artifacts = int(metadata_line.split(' ')[-1])


        n_nodos = []
        for i in range(n_artifacts):
            n_nodos.append(f.readline().rstrip("\n"))

        for i in range(n_artifacts):
            coordinate_nodos[n_nodos[i]] = [float(coor) for coor in f.readline().split(' ')[:-1]]

    f.close()

    return coordinate_nodos


def get_conditions():
    conditions = {}

    f = open('./elec_nodes.txt', 'r')
    index = 0

    meta = [
        ['dirichlet', 0],
        ['dirichlet', 0],
        ['dirichlet', 1],
        
        ['dirichlet', 1],
        ['dirichlet', 1],
        ['dirichlet', 1],
        ['dirichlet', 0],

        ['dirichlet', 0],
        ['dirichlet', 0],

        
        ['dirichlet', 0],
        
        ['dirichlet', 0],

        ['dirichlet', 0],

        
        ['neumann', 0],

        ['dirichlet', 0],

        
        ['dirichlet', 0],
        ['dirichlet', 1],
        ['dirichlet', 1],


        ['E', 2],
        ['E', 1]
    ]

    while True:

        metadata_line = f.readline().rstrip("\n")
        if metadata_line == '':
            break

        n_artifacts = int(metadata_line.split(' ')[-1])

        for i in range(n_artifacts):
            conditions[f.readline().rstrip("\n")] = meta[index]

        for _ in range(n_artifacts):
            f.readline()

        index += 1

    f.close()

    return conditions



def get_elements():
    elements = []
    
    coordinate_nodos = get_coor_nodes()
    conection_nodos = get_conection_nodes()
    conditions = get_conditions()

    for key, value in conection_nodos.items():

        position = []
        condition = []
        conditions_values = []

        for nodo in value:
            position.append(np.array(coordinate_nodos[str(nodo)]))
            condition.append(conditions[str(nodo)][0])
            conditions_values.append(conditions[str(nodo)][1])

        element = {
            "n_nodo": value,
            "nodos": position,
            "condition": condition,
            "condition_value": conditions_values
        }

        element = get_funtion_form(element)

        elements.append(element)

    return elements

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

        n_nodes_group_line = n_nodes_group_line + n_nodes_group

        # print(group_first_line, group_first_line + n_nodes_group)

        conections.append(group_conections)

    return conections


def create_elements(nodes, conections, conditions):
    elements = []

    for index_group, group_elements in enumerate(conections):
        for element_nodes in group_elements:

            element = {
                "n_nodo": element_nodes,
                "nodos": [nodes[str(node)] for node in element_nodes[1:]],
                "condition": conditions[index_group][0],
                "condition_value": conditions[index_group][1]
            }


            element = get_funtion_form(element)
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


    print(elements)

    return elements





conditions = [
        ['dirichlet', 0],
        ['dirichlet', 0],
        ['dirichlet', 1],
        
        ['dirichlet', 1],
        ['dirichlet', 1],
        ['dirichlet', 1],
        ['dirichlet', 0],

        ['dirichlet', 0],
        ['dirichlet', 0],

        
        ['dirichlet', 0],
        
        ['dirichlet', 0],

        ['dirichlet', 0],

        
        ['neumann', 0],



        
        ['dirichlet', 1],
        ['dirichlet', 1],


        ['E', 2],
        ['E', 1]
]


if __name__ == '__main__':
    filename = './mStrip.msh'
    read_msh(filename)


    elements = get_elements()
    K = np.zeros((727, 727))
    F = np.zeros((727, 1))
    # main()

