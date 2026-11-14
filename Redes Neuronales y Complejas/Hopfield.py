import numpy as np
import matplotlib.pyplot as plt

### 1. INICIALIZACIÓN DE PARÁMETROS Y VARIABLES DE CONTROL
np.random.seed(42)  # para reproducibilidad
noise_level = 0.4  # % de ruido
max_iters = 1024 # número máximo de iteraciones para la dinámica de Hopfield

# Matriz 16x16 de la letra A (+1 activo, -1 fondo)
patron_A = np.array([
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1,  1,  1, -1, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1, -1, -1],
    [-1, -1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1, -1, -1],
    [-1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1],
    [-1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1],
    [-1, -1, -1,  1,  1, -1, -1, -1, -1, -1, -1,  1,  1, -1, -1, -1],
    [-1, -1, -1,  1,  1, -1, -1, -1, -1, -1, -1,  1,  1, -1, -1, -1],
    [-1, -1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1, -1, -1],
    [-1, -1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1, -1, -1],
    [-1, -1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
])

xi_A = patron_A.flatten()  # Aplanar la matriz a un vector de 256 componentes

# visualización del patrón:
# plt.imshow(patron_A, cmap='gray')
# plt.show()

###########################################################################################

### 2. Asignación de pesos sinápticos (matriz de pesos, W, simétrica, con ceros en la diagonal)
# Este proceso es el de "memorización" mediante la regla de Hebb

N = int(patron_A.size)
print(f"Número de neuronas: {N}")
J = np.zeros([N,N])

for i in range(N):
    for j in range(N):
        if i == j:
            continue
        else:
            J[i, j] = 1 / 256 * np.sum(xi_A[i] * xi_A[j]) # fórmula simplficada de la regla de Hebb (solo tenemos un patrón)
print("Asignación de pesos sinápticos terminada")


# plt.imshow(J, cmap='viridis')
# plt.show()



###########################################################################################


### 3. Recuperación del patrón a partir de un patrón ruidoso


# Generamos un patrón ruidoso a partir del patrón A
# tipo de ruido: cambiar el signo de un porcentaje de los bits del patrón original

print(f"[INFO] ----- Generando patrón ruidoso con {noise_level*100}% de ruido")
xi_A_noisy = xi_A.copy()
for i in range(len(xi_A)):
    if np.random.rand() < noise_level:
        xi_A_noisy[i] = -xi_A[i]

patron_A_noisy = xi_A_noisy.reshape(patron_A.shape)

# plt.figure()
# plt.subplot(1, 2, 1)
# plt.imshow(patron_A, cmap='gray')
# plt.title('Patrón original A')
# plt.subplot(1, 2, 2)
# plt.imshow(patron_A_noisy, cmap='gray')
# plt.title('Patrón ruidoso A')
# plt.show()


# Ahora aplicamos la dinámica de evolución de Hopfield SOBRE LAS NEURONAS
# los pesos sinápticos se mantienen fijos, los dados por la regla de Hebb

# Dinámica de evolución asíncrona (para asegurar el decrecimiento del funcional de Lyapunov)
def hopfield_dynamics(xi, J, max_iter=max_iters):
    print("[INFO] ----- Iniciando dinámica de Hopfield")
    for iteration in range(max_iter):

        # selección de UNA neurona aleatoria para actualizar en cada iteración
        random_neuron = np.random.randint(0, N)
        i = random_neuron

        xi_prev = xi.copy() # guardamos el estado previo 

        # actualización de la neurona i según la regla de activación
        h_i = np.dot(J[i, :], xi)  # campo local sobre la neurona i
        xi[i] = 1 if h_i >= 0 else -1  # regla de activación sobre la neurona i
        # print(f"Iteración {iteration}: Neurona {i} actualizada de {xi_prev[i]} a {xi[i]}")
        if iteration == max_iter // 4:
            print(f"[INFO] ----- 25 % de las iteraciones alcanzadas ({iteration})")
            xi_25 = xi.copy()
        elif iteration == max_iter // 2:
            print(f"[INFO] ----- Mitad de iteraciones alcanzada ({iteration})")
            xi_50 = xi.copy()
        elif iteration == 3 * max_iter // 4:
            print(f"[INFO] ----- 75 % de las iteraciones alcanzadas ({iteration})")
            xi_75 = xi.copy()

    print(f"[CORRECT] ----- Dinámica de Hopfield finalizada ({max_iter} iteraciones)")
    return xi, xi_25, xi_50, xi_75

hopfield_output, hopfield_output_25, hopfield_output_50, hopfield_output_75 = hopfield_dynamics(xi_A_noisy.copy(), J)



###########################################################################################


### 4. VISUALIZACIÓN DE LOS RESULTADOS

# EVOLUCIÓN TEMPORAL DE LA DINÁMICA DE HOPFIELD
patron_A_recovered    = hopfield_output.reshape(patron_A.shape)
patron_A_recovered_25 = hopfield_output_25.reshape(patron_A.shape)
patron_A_recovered_50 = hopfield_output_50.reshape(patron_A.shape)
patron_A_recovered_75 = hopfield_output_75.reshape(patron_A.shape)

imagenes = [
    patron_A_noisy,
    patron_A_recovered_25,
    patron_A_recovered_50,
    patron_A_recovered_75,
    patron_A_recovered
]

titulos = [
    '0% (Entrada)',
    '25% Iteraciones',
    '50% Iteraciones',
    '75% Iteraciones',
    '100% Convergencia'
]

fig, axes = plt.subplots(1, 5, figsize=(16, 5))

for ax, img, title in zip(axes, imagenes, titulos):
    ax.imshow(img, cmap='gray')
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.axis('off')  # Oculta los ejes numéricos (ticks) para una vista limpia

plt.suptitle("Evolución Temporal de la Red de Hopfield", fontsize=14, fontweight='bold', y=1.05)
plt.tight_layout()
plt.show()


# COMPARACIÓN ORIGINAL-INPUT-OUTPUT
plt.figure(figsize=(16, 5))

plt.subplot(1, 3, 1)
plt.imshow(patron_A, cmap='gray')
plt.title('Patrón original A')
plt.subplot(1, 3, 2)
plt.imshow(patron_A_noisy, cmap='gray')
plt.title('Patrón ruidoso A')
plt.subplot(1, 3, 3)
plt.imshow(patron_A_recovered, cmap='gray')
plt.title('Patrón recuperado A')

plt.tight_layout()
plt.suptitle("Comparación Original-Input-Output", fontsize=14, fontweight='bold', y=1.05)
plt.show()

###########################################################################################
