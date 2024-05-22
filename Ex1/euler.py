import numpy as np

def euler_leg(h, f0):
    """
    Implementação do método de Euler para resolver um sistema de equações diferenciais 
    descrevendo o movimento de um pêndulo duplo.

    Parâmetros:
    - h (float): tamanho do passo de integração.
    - f0 (array): vetor contendo as condições iniciais [theta1, theta2, theta1_dot, theta2_dot].

    Retorna:
    - f_out (array): vetor contendo os valores atualizados das variáveis no tempo t + h.
    """

    # Constantes
    m1 = 0.5
    m2 = 1
    l1 = 1
    l2 = 1.5
    g = 9.8

    # Conversão para variáveis da função
    teta1, teta2, teta11, teta21 = f0

    # Equações de movimento
    num_g1 = (-np.sin(teta1 - teta2) * (m2 * l2 * (teta21 ** 2) + m2 * l1 * (teta11 ** 2) * np.cos(teta1 - teta2))- g * ((m1 + m2) * np.sin(teta1) - m2 * (np.sin(teta2) * np.cos(teta1 - teta2))))
    den_g1 = l1 * (m1 + m2 * ((np.sin(teta1 - teta2)) ** 2))
    g1 = num_g1 / den_g1

    num_g2 = (np.sin(teta1 - teta2) * ((m1 + m2) * l1 * (teta11 ** 2) + m2 * l2 * (teta21 ** 2) * np.cos(teta1 - teta2))+ g * ((m1 + m2) * (np.sin(teta1) * np.cos(teta1 - teta2) - np.sin(teta2))))
    den_g2 = l2 * (m1 + m2 * ((np.sin(teta1 - teta2)) ** 2))
    g2 = num_g2 / den_g2

    # Vetor de saída
    fxy = np.array([teta11, teta21, g1, g2, 0, 0])

    f_out = f0 + h * fxy

    return f_out

def euler_steps(duration, h, f0):
    """
    Realiza a integração numérica usando o método de Euler.

    Parâmetros:
    - duration (float): tempo total de integração.
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
        f0 = euler_leg(h, f0)

        time += 1
        response[i, 0] = time * h
        response[i, 1:] = f0

    return response


