import numpy as np

def rotate(ponto,  angle): #função que rotaciona ponto e o ajusta para ponto a extrema esquerda

    angle = angle*np.pi/180

    rotate = np.array([
      [np.cos(angle), np.sin(angle)],
      [-np.sin(angle), np.cos(angle)]
  ])


    ponto.x, ponto.y = rotate@np.array([ponto.x, ponto.y])

    ponto.x = ponto.x +5.5
    ponto.y = ponto.y +4

    return ponto

def vetor_normal(p1, p2):
    # Coordenadas dos pontos
    x1, y1 = p1.x, p1.y
    x2, y2 = p2.x, p2.y

    # Vetor direcional (x2 - x1, y2 - y1)
    dx = x2 - x1
    dy = y2 - y1

    # Vetor normal (dy, -dx)
    normal = [-dy, dx]
    normal = normal/np.linalg.norm(normal)
    return normal

def dist_ponto_segmento(ponto, a, b):
    """
    Distância do ponto ao segmento de reta a-b, usada em update_normal
    para achar a aresta da placa geometricamente mais próxima de cada
    ponto de fluido (em vez de decidir por prioridade fixa de direção).
    """
    ax, ay = a.x, a.y
    bx, by = b.x, b.y
    abx, aby = bx - ax, by - ay
    comprimento2 = abx**2 + aby**2
    t = ((ponto.x - ax)*abx + (ponto.y - ay)*aby) / comprimento2
    t = max(0.0, min(1.0, t))
    cx, cy = ax + t*abx, ay + t*aby
    return np.hypot(ponto.x - cx, ponto.y - cy)