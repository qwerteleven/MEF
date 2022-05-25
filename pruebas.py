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


def condiciones_dirichlet(K, F, n_nodos):

    for n_nodo, value in n_nodos:
        G = 10e6
        K[n_nodo, n_nodo] = G
        F[n_nodo] = F[n_nodo] * G

    return K, F


def get_k(index_n, index_m, nodos, formas):
    
    lambda_onda = 1
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

    k = 0
    masa = 0

    determinante = det_j(nodos)
    jacobiano_inv = j_tras_inv(nodos)

    for index, nodo in enumerate(nodos):
        k = (((jacobiano_inv @ gr[index_n]) * (jacobiano_inv @ gr[index_m]) * determinante) @ nodo) * w[index]
    
    for index, nodo in enumerate(nodos):
        masa += np.dot(formas[index_n], formas[index_m]) * determinante * w[index] 

    return k - ((2 * np.pi) / lambda_onda) ** 2 * masa


def get_f(index_n, index_m, forma, nodos):
    ro = 1
    angular = 1
    j = 1

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

    determinante = det_j(nodos)
    f = 0

    for _, nodo in enumerate(nodos):
        f += (f_neumann(j, ro, angular, gr[index_n]) * forma[index_m]) * determinante @ nodo * w[index_n]

    return (-j) * ro * angular * f


def f_neumann(j, ro, angular, forma):
    v = (- forma) / (j * angular * ro)
    return v


def jacobiano(n_1, n_2, n_3):
    return 2 * get_ae(n_1, n_2, n_3)


def ensamblaje_k_f(element, K, F):

    for index_n in range(len(element['forma'])):
        for index_m in range(len(element['forma'])):
            K[element['n_nodo'][index_n], element['n_nodo'][index_m]] += get_k(index_n, index_m, element['nodos'], element['forma'])
            F[element['n_nodo'][index_n]] += get_f(index_n, index_m, element['forma'], element['nodos']) 

    return K, F


def get_sistem(elements, n_nodos_dir):

    K = np.zeros((len(elements), len(elements)))
    F = np.zeros((len(elements), 1))

    for element in elements:
        K, F = ensamblaje_k_f(element, K, F)

    K, F = condiciones_dirichlet(K, F, n_nodos_dir)

    return K, F



def main():
    elements = get_elements()
    n_nodos_dir = []

    for element in elements:
        for index in range(len(element['n_nodo'])):
            if element['condition'][index] == 'dirichlet':
                n_nodos_dir.append((element['n_nodo'][index], element['condition_value'][index]))

    n_nodos_dir = list(set(n_nodos_dir))

    print(n_nodos_dir)

    K, F = get_sistem(elements, n_nodos_dir)

    import json

    with open('test.json', 'w') as f:
        lists = K.tolist()
        json_str = json.dumps(lists, indent = 2)
        f.write(json_str)
    
    
    X2 = np.linalg.solve(K, F)


def get_conection_nodes():

    conection_nodos = {}

    f = open('./elements.txt', 'r')

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

    f = open('./nodes.txt', 'r')

    while True:

        metadata_line = f.readline().rstrip("\n")
        if metadata_line == '':
            break

        n_artifacts = int(metadata_line.split(' ')[-1])


        n_nodos = []
        for i in range(n_artifacts):
            n_nodos.append(f.readline().rstrip("\n"))

        for i in range(n_artifacts):
            coordinate_nodos[n_nodos[i]] = [float(coor) for coor in f.readline().split(' ')]

    f.close()

    return coordinate_nodos


def get_adj_nodes():

    adj_nodes = {}
    conection_nodos = get_conection_nodes()

    for n_node in range(1, 958 + 1):

        adj_nodes[str(n_node)] = []
        for key, value in conection_nodos.items():
            for v in value:
                if v == n_node:
                    adj_nodes[str(n_node)] = adj_nodes[str(n_node)] + value

        set_value = set(adj_nodes[str(n_node)])
        set_value.remove(n_node)
        
        adj_nodes[str(n_node)] = list(set_value)

    return adj_nodes


def get_conditions():
    conditions = {}

    f = open('./nodes.txt', 'r')
    index = 0

    meta = [
        ['dirichlet', 0],
        ['dirichlet', 0],
        ['dirichlet', 0],
        ['dirichlet', 0],

        ['neumann', 0],
        ['neumann', 0],
        ['neumann', 0],
        ['neumann', 0],

        ['dirichlet', 0],

        ['neumann', 0],

        ['dirichlet', 0],

        ['dirichlet', 1],


        ['neumann', 0],
        ['neumann', 0],
        ['neumann', 0],
        ['neumann', 0],
        ['none', 0]
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

if __name__ == '__main__':
    main()

