import numpy as np
import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8-whitegrid')

COLOR_M1 = '#d62728'
COLOR_M2 = '#1f77b4'
COLOR_DIFF = '#d62728'


def plot_values(response, metodo, passo):
    """
    Plota deslocamento, velocidade e aceleração angulares do pêndulo duplo.

    Gera uma grade 3x2 (linhas: deslocamento / velocidade / aceleração;
    colunas: massa 1 / massa 2), com uma única curva por painel, para que
    cada trajetória fique legível mesmo quando m1 e m2 têm amplitudes
    parecidas. Ambas as massas aparecem na mesma figura, conforme pedido
    no enunciado.

    Parâmetros:
    - response (array): matriz [tempo, teta1, teta2, teta1_dot, teta2_dot,
      teta1_ddot, teta2_ddot] retornada por euler_steps/rk4_steps.
    - metodo (str): nome do método usado (para o título).
    - passo (float): passo de integração h usado (para o título).
    """

    t = response[:, 0]
    rows = [
        ('Angular displacement', r'$\theta$ [rad]', response[:, 1], response[:, 2]),
        ('Angular velocity', r'$\dot{\theta}$ [rad/s]', response[:, 3], response[:, 4]),
        ('Angular acceleration', r'$\ddot{\theta}$ [rad/s²]', response[:, 5], response[:, 6]),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
    fig.suptitle(f'{metodo} — h = {passo}', fontsize=15, fontweight='bold')

    axes[0, 0].set_title('Mass 1')
    axes[0, 1].set_title('Mass 2')

    for row, (name, ylabel, v1, v2) in enumerate(rows):
        axes[row, 0].plot(t, v1, color=COLOR_M1, lw=1.2)
        axes[row, 0].set_ylabel(f'{name}\n{ylabel}')
        axes[row, 1].plot(t, v2, color=COLOR_M2, lw=1.2)

    for ax in axes.flat:
        ax.margins(x=0)
    for ax in axes[-1, :]:
        ax.set_xlabel('Time [s]')

    fig.tight_layout(rect=[0, 0, 1, 0.96])

    plt.show()


def plot_comparison(t_a, y_a, t_b, y_b, label_a, label_b, ylabel, title, ax=None):
    """
    Compara duas soluções para a mesma variável (ex.: Euler x RK4, ou dois
    passos h diferentes), com uma única linha por solução e a área entre
    as curvas sombreada para evidenciar a diferença ponto a ponto.

    As duas soluções podem ter sido calculadas com passos h diferentes
    (logo, vetores de tempo de tamanhos distintos); nesse caso a série b é
    interpolada sobre os instantes de t_a antes de calcular a diferença.

    Parâmetros:
    - t_a, y_a (array): tempo e valores da variável para a solução A.
    - t_b, y_b (array): tempo e valores da variável para a solução B.
    - label_a, label_b (str): rótulos das duas soluções na legenda.
    - ylabel (str): rótulo do eixo y (com unidade).
    - title (str): título do gráfico.
    - ax (matplotlib.axes.Axes, opcional): eixo onde plotar; se None,
      cria uma figura nova e a exibe.
    """

    own_fig = ax is None
    if own_fig:
        _, ax = plt.subplots(figsize=(9, 4))

    y_b_on_a = np.interp(t_a, t_b, y_b)

    ax.plot(t_a, y_a, color=COLOR_M1, lw=1.2, label=label_a)
    ax.plot(t_b, y_b, color=COLOR_M2, lw=1.2, label=label_b)
    ax.fill_between(t_a, y_a, y_b_on_a, color=COLOR_DIFF, alpha=0.15, label='Difference')

    ax.set_title(title)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(ylabel)
    ax.margins(x=0)
    ax.legend(loc='upper left', framealpha=0.9)

    if own_fig:
        plt.tight_layout()
        plt.show()

    return ax
