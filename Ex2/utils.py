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