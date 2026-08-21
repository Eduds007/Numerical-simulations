import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8-whitegrid')

COLOR_M1 = '#d62728'
COLOR_M2 = '#1f77b4'


def plot_values(response, metodo, passo):
    """
    Plota deslocamento, velocidade e aceleração angulares do pêndulo duplo.

    Gera três gráficos empilhados (um por página, um debaixo do outro),
    cada um contendo as curvas de m1 e m2 sobrepostas, conforme pedido
    no enunciado.

    Parâmetros:
    - response (array): matriz [tempo, teta1, teta2, teta1_dot, teta2_dot,
      teta1_ddot, teta2_ddot] retornada por euler_steps/rk4_steps.
    - metodo (str): nome do método usado (para o título).
    - passo (float): passo de integração h usado (para o título).
    """

    t = response[:, 0]
    theta1, theta2 = response[:, 1], response[:, 2]
    theta1_dot, theta2_dot = response[:, 3], response[:, 4]
    theta1_ddot, theta2_ddot = response[:, 5], response[:, 6]

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    fig.suptitle(f'{metodo} — h = {passo}', fontsize=15, fontweight='bold')

    axes[0].plot(t, theta1, color=COLOR_M1, lw=1.3, label=r'$\theta_1$ (m1)')
    axes[0].plot(t, theta2, color=COLOR_M2, lw=1.3, label=r'$\theta_2$ (m2)')
    axes[0].set_title('Angular displacement')
    axes[0].set_ylabel(r'$\theta$ [rad]')

    axes[1].plot(t, theta1_dot, color=COLOR_M1, lw=1.3, label=r'$\dot{\theta}_1$ (m1)')
    axes[1].plot(t, theta2_dot, color=COLOR_M2, lw=1.3, label=r'$\dot{\theta}_2$ (m2)')
    axes[1].set_title('Angular velocity')
    axes[1].set_ylabel(r'$\dot{\theta}$ [rad/s]')

    axes[2].plot(t, theta1_ddot, color=COLOR_M1, lw=1.3, label=r'$\ddot{\theta}_1$ (m1)')
    axes[2].plot(t, theta2_ddot, color=COLOR_M2, lw=1.3, label=r'$\ddot{\theta}_2$ (m2)')
    axes[2].set_title('Angular acceleration')
    axes[2].set_ylabel(r'$\ddot{\theta}$ [rad/s²]')
    axes[2].set_xlabel('Time [s]')

    for ax in axes:
        ax.legend(loc='upper right', framealpha=0.9)
        ax.margins(x=0)

    fig.tight_layout(rect=[0, 0, 1, 0.96])

    plt.show()
