import numpy as np

def rk4_leg(f0):
    """
    Implementação do método de Runge-Kutta de quarta ordem para resolver um sistema de equações diferenciais 
    descrevendo o movimento de um pêndulo duplo.

    Parâmetros:
    - f0 (array): vetor contendo as condições iniciais [theta1, theta2, theta1_dot, theta2_dot].

    Retorna:
    - fxy (array): vetor contendo as derivadas das variáveis do sistema.
    """

    # Constantes
    m1 = 0.5
    m2 = 1
    l1 = 1
    l2 = 1.5
    g = 9.8

    # Conversão para variáveis da função
    teta1, teta2, teta11, teta21 = f0

    # Variáveis para cálculo
    num_g1 = (-np.sin(teta1 - teta2) * (m2 * l2 * (teta21 ** 2) + m2 * l1 * (teta11 ** 2) * np.cos(teta1 - teta2))- g * ((m1 + m2) * np.sin(teta1) - m2 * (np.sin(teta2) * np.cos(teta1 - teta2))))
    den_g1 = l1 * (m1 + m2 * ((np.sin(teta1 - teta2)) ** 2))
    g1 = num_g1 / den_g1

    num_g2 = (np.sin(teta1 - teta2) * ((m1 + m2) * l1 * (teta11 ** 2) + m2 * l2 * (teta21 ** 2) * np.cos(teta1 - teta2))+ g * ((m1 + m2) * (np.sin(teta1) * np.cos(teta1 - teta2) - np.sin(teta2))))
    den_g2 = l2 * (m1 + m2 * ((np.sin(teta1 - teta2)) ** 2))
    g2 = num_g2 / den_g2

    # Vetor de saída
    fxy = np.array([teta11, teta21, g1, g2, 0, 0])

    return fxy

def rk4_steps(duration, h, f0):
    """
    Realiza a integração numérica usando o método de Runge-Kutta de quarta ordem.

    Parâmetros:
    - duration (float): tempo total em segundos.
    - h (float): tamanho do passo de integração.
    - f0 (array): vetor contendo as condições iniciais [theta1, theta2, theta1_dot, theta2_dot].

    Retorna:
    - response (array): matriz contendo os valores de tempo e as variáveis do sistema ao longo do tempo.
    """

    time = 0
    steps = int(duration / h)
    response = np.zeros((steps, 7))
    response[0, 1:] = f0
    response[0, 0] = 0

    for i in range(1, steps):
        k1 = rk4_leg(f0)
        f1 = f0 + h / 2 * k1
        k2 = rk4_leg(f1)
        f2 = f0 + h / 2 * k2
        k3 = rk4_leg(f2)
        f3 = f0 + h * k3
        k4 = rk4_leg(f3)

        f0 = f0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

        # Atualiza as derivadas no vetor f0
        f0[4] = k1[2]
        f0[5] = k1[3]

        time += 1
        response[i, 0] = time * h
        response[i, 1:] = f0

    return response
