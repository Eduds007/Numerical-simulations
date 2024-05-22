import matplotlib.pyplot as plt

def plot_values(response, metodo, passo):

    plt.figure(figsize=(10, 6))

    plt.subplot(3, 2, 1)  # Subplot 1
    plt.plot(response[:,0], response[:,1], 'r-',)
    plt.title(f'[{metodo}] Posição massa 1 | h={passo}')

    plt.subplot(3, 2, 2)  # Subplot 2
    plt.plot(response[:,0], response[:,2], 'b-',)
    plt.title(f'[{metodo}] Posição massa 2 | h={passo} ')

    plt.subplot(3, 2, 3)  # Subplot 3
    plt.plot(response[:,0], response[:,2], 'r-',)
    plt.title(f'[{metodo}] Velocidade massa 1 | h={passo} ')

    plt.subplot(3, 2, 4)  # Subplot 4
    plt.plot(response[:,0], response[:,3], 'b-',)
    plt.title(f'[{metodo}] Velocidade massa 2 | h={passo}')

    plt.subplot(3, 2, 5)  # Subplot 5
    plt.plot(response[:,0], response[:,4], 'r-',)
    plt.title(f'[{metodo}] Aceleração massa 1 | h={passo} ')

    plt.subplot(3, 2, 6)  # Subplot 6
    plt.plot(response[:,0], response[:,5], 'b-',)
    plt.title(f'[{metodo}] Aceleração massa 2 | h={passo}')

    plt.tight_layout()  

    plt.show()