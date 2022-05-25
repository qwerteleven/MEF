import numpy as np

malla = [
    [[0, 0],[1, 0],[0, 1]],
    [[2, 2],[1, 2],[2, 1]],
    [[0, 0],[-1, 0],[0, -1]],
    [[0, 0],[1, 0],[0, 1]],
]

nodos = [
    [0, 0],
    [1, 0],
    [0, 1],
    [1, 2]
]


def get_form_funtions(malla):

    form_funtions = [None] * len(malla)

    B = np.zeros((len(malla[0]), len(malla[0])), int)
    np.fill_diagonal(B, 1)

    for index, element in enumerate(malla):
        A = np.column_stack((np.ones((len(element)), dtype=int), np.array(element)))
        X2 = np.linalg.solve(A, B)
        form_funtions[index] = X2.T

    return form_funtions


form_funtions = get_form_funtions(malla)

print(form_funtions)



first_base_funtion = form_funtions[2]

point = [-0.5, -0.5]
e = 0
for base in first_base_funtion:
    e = base[0] + base[1] * point[0] + base[2] * point[1]
    print('e: ', e)















x = np.linspace(0, 2, 101)
y = np.linspace(0, 1, 101)
xs, ys = np.meshgrid(x, y, sparse=True)
zs = np.sqrt(xs**2 + ys**2)


import matplotlib.pyplot as plt

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

surf = ax.plot_surface(xs, ys, zs, cmap=plt.cm.coolwarm)
ax.set_title('3D line plot geeks for geeks')
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()





# forma  -> a + bx + cy
# x, y -> son el punto dentro del elemento de malla
# problema determinar a, b, c -> 
# solucion sistema de ecuaciones para sacar los coeficientes
# AX = B

# A -> posiciones x, y

# A = [
# [0, 0],
# [0, 1],
# [1, 0]
# ]

# X -> coeficientes

# B -> ecuaciones linealmente independientes
# matriz diagonal
