import numpy as np
import matplotlib.pyplot as plt
from utils import *



class Point:
    def __init__(self, i, j, u_inf, delta, uy, corrente, p_init) -> None:

        self.x = i*delta
        self.y = j*delta
        self.is_rect = False
        self.corrente = corrente
        self.ux = u_inf
        self.uy =uy
        self.pressure = p_init
        self.normal = [0,0]


    def __str__(self):
        return f"({self.x}, {self.y})"

class Retangle(Point):
  def __init__(self, x1, y1,angle, u_inf, delta, uy, corrente, p_init, h, w ):
    
    delta = 1
    self.p2 = rotate(Point(x1,y1+w,u_inf, delta, uy, corrente, p_init), angle)   
    self.p3 = rotate(Point(x1+h,y1+w,u_inf, delta, uy, corrente, p_init), angle)
    self.p1 = rotate(Point(x1,y1,u_inf, delta, uy, corrente, p_init), angle)
    self.p4 = rotate(Point(x1+h,y1,u_inf, delta, uy, corrente, p_init), angle)

    self.perimeter = 2*h + 2*w
    self.angle = angle

  def is_inside(self, ponto):
    def is_left(p1, p2, p3):
        return (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y)

    def is_in_triangle(p, p1, p2, p3):
        b1 = is_left(p, p1, p2) <= 0.0
        b2 = is_left(p, p2, p3) <= 0.0
        b3 = is_left(p, p3, p1) <= 0.0
        return (b1 == b2) and (b2 == b3)

    # Verifica se o ponto está dentro do retângulo
    return is_in_triangle(ponto, self.p1, self.p2, self.p3) or is_in_triangle(ponto, self.p1, self.p3, self.p4)

class Mesh():
  def __init__(self, rect, u_inf, delta, tolerance, lamb, h, corrente_inicial, uy, p_init, iteracao, rho, g):
    self.linhas = int(8/delta)+1
    self.colunas = int(11/delta)+1
    self.rect = rect
    self.tolerance = tolerance
    self.delta = delta
    self.u_inf = u_inf
    self.lamb = lamb
    self.iter = iteracao
    self.ux = u_inf
    self.uy = uy
    self.corrente_inicial = corrente_inicial
    self.p_init = p_init
    self.rho = rho
    self.g = g

    self.matriz = np.empty((self.linhas, self.colunas), dtype=object)
    psi = np.tile(np.linspace(self.ux*h,-self.ux*h, self.linhas), (self.colunas, 1)).T

    # Preencher a matriz com objetos da classe Ponto
    for i in range(self.linhas):
        for j in range(self.colunas):
            self.matriz[i, j] = Point(j,self.linhas-1-i, self.u_inf, delta, self.uy, self.corrente_inicial, self.p_init)
    self.corrente =  np.array([[ponto.corrente for ponto in linha] for linha in self.matriz])

  def update_mask(self):
    shape = self.matriz.shape
    for i in range(1,shape[0]-1):
      for j in range(1,shape[1]-1):
          if self.rect.is_inside(self.matriz[i, j]):
            self.matriz[i,j].is_rect = True

  def update_velocidade(self):
    shape = self.matriz.shape
    for i in range(1,shape[0]-1):

        for j in range(1,shape[1]-1):
          if self.rect.is_inside(self.matriz[i, j]):
            self.matriz[i,j].ux = 0
            self.matriz[i,j].uy = 0
            
          else:
            self.matriz[i,j].uy = (self.matriz[i,j+1].corrente - self.matriz[i,j-1].corrente)/(2*self.delta)
            self.matriz[i,j].ux = (self.matriz[i-1,j].corrente - self.matriz[i+1,j].corrente)/(2*self.delta)

  def update_pressure(self):
    shape = self.matriz.shape
    for i in range(1,shape[0]-1):

        for j in range(1,shape[1]-1):

          if self.rect.is_inside(self.matriz[i, j]):

            self.matriz[i, j].pressure = self.p_init
          else:
            self.matriz[i, j].pressure = self.rho*self.g*(self.matriz[0, 0].pressure/(self.rho*self.g) + (self.matriz[0, 0].ux)**2/(2*self.g) - (self.matriz[i, j].ux**2 + self.matriz[i, j].uy**2)/(2*self.g))


  def update_corrente(self):

    shape = self.matriz.shape
    print(shape)

    for k in range(self.iter):
      aux = self.matriz.copy()
      erro = []
      for i in range(0,shape[0]):

        for j in range(0,shape[1]):
          if self.rect.is_inside(self.matriz[i, j]):
            self.matriz[i,j].corrente = 0

          elif (not (i==0 or i==shape[0]-1)) and (j==0): #Entrada
            #print(i,j)
            #print('entrada')
            novo = self.lamb*(2*aux[i,j+1].corrente +  aux[i+1,j].corrente +  aux[i-1,j].corrente)/4 + (1-self.lamb)*aux[i,j].corrente
            
            self.matriz[i,j].corrente = novo
            

          elif (not(i==0 or i==shape[0]-1)) and (j==shape[1]-1): #Saída
            #print(i,j)
            #print('saída')
            novo = self.lamb*(2*aux[i,j-1].corrente +  aux[i+1,j].corrente +  aux[i-1,j].corrente)/4 + (1-self.lamb)*aux[i,j].corrente
            try:
              erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
            except ZeroDivisionError:
              pass
            self.matriz[i,j].corrente = novo
            pass

          elif (i==0) and (not(j==0 or j==shape[1]-1)): #Superior
            #print(i,j)
            #print('superior')
            novo = self.lamb*(aux[i,j+1].corrente +  aux[i,j-1].corrente +  2*aux[i+1,j].corrente + 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
            try:
              erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
            except ZeroDivisionError:
              pass
            self.matriz[i,j].corrente = novo

          elif (i==shape[0]-1) and (not(j==0 or j==shape[1]-1)): #Inferior
            #print(i,j)
            #print('inferior')
            novo = self.lamb*(aux[i,j+1].corrente +  aux[i,j-1].corrente +  2*aux[i-1,j].corrente - 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
            try:
              erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
            except ZeroDivisionError:
              pass
            self.matriz[i,j].corrente = novo

          elif (i==0 or i==shape[0]-1) and (j==0 or j==shape[1]-1): #Canto

            print(i,j)
            #print('canto')
            if (i==0 and j==0):
              #print('canto superior esquerdo')
              novo = self.lamb*(2*aux[i+1,j].corrente +  2*aux[i,j+1].corrente + 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo

            elif (i==shape[0]-1 and j==0): #canto inferior esquerdo
              #print('canto inferior esquerdo')
              novo = self.lamb*(2*aux[i-1,j].corrente +  2*aux[i,j+1].corrente - 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo

            elif (i==0 and j==shape[1]-1): #canto superior direito
              #print('canto superior direito')
              novo = self.lamb*(2*aux[i+1,j].corrente +  2*aux[i,j-1].corrente + 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo
            elif (i==shape[0]-1 and j==shape[1]-1): #canto inferior direito
              #print('canto inferior direito')
              novo = self.lamb*(2*aux[i-1,j].corrente +  2*aux[i,j-1].corrente - 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo
            #print('canto')

            pass
          #elif (i==0 or i==shape[1]) and (j==0 or i==shape[0]): #Cantos da placa
            #pass
          #elif (self.rect.is_inside(self.matriz[i+1, j]) or self.rect.is_inside(self.matriz[i-1, j]) or self.rect.is_inside(self.matriz[i, j+1]) or self.rect.is_inside(self.matriz[i, j-1])): #borda irregular
            #pass



          else:
            #print(i,j)
            #print(f'antes: { self.matriz[i,j].corrente}')
            #print(i,j)
            novo = self.lamb*(aux[i+1,j].corrente +  aux[i-1,j].corrente +  aux[i,j+1].corrente +  aux[i,j-1].corrente)/4 + (1-self.lamb)*aux[i,j].corrente
            #print(novo)
            try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
            except ZeroDivisionError:
                pass
            self.matriz[i,j].corrente = novo
            #print(f'depois: { self.matriz[i,j].corrente}')
      self.corrente =  np.array([[ponto.corrente for ponto in linha] for linha in self.matriz])
      print(f'Iteração: {k}')
      print(f'Erro: {np.max(erro)}')
      print()
      if np.max(erro) <= self.tolerance:
        break

  def update_normal(self):
    shape = self.matriz.shape
    if self.rect.angle == 15:
      for i in range(1,shape[0]-1):
        for j in range(1,shape[1]-1):
          if not self.matriz[i,j].is_rect:
            #cima
            if self.matriz[i+1,j].is_rect:
              #print(vetor_normal(rect.p4, rect.p1))
              self.matriz[i,j].normal = vetor_normal(self.rect.p2, self.rect.p3)
            #baixo
            elif self.matriz[i-1,j].is_rect:
               self.matriz[i,j].normal = vetor_normal(self.rect.p4, self.rect.p1)
            #esquerda
            elif self.matriz[i,j+1].is_rect:
              self.matriz[i,j].normal = vetor_normal(self.rect.p1, self.rect.p2)
            #dieita
            elif self.matriz[i,j-1].is_rect:
              self.matriz[i,j].normal = vetor_normal(self.rect.p3, self.rect.p4)
    elif self.rect.angle == 90:
      for i in range(1,shape[0]-1):
        for j in range(1,shape[1]-1):
          if not self.matriz[i,j].is_rect:
            #cima
            if self.matriz[i+1,j].is_rect and (not self.matriz[i,j-1].is_rect):
              self.matriz[i,j].normal = vetor_normal(self.rect.p1, self.rect.p2)
            #esquerda
            elif self.matriz[i,j+1].is_rect:
              #print(vetor_normal(rect.p4, rect.p1))
              self.matriz[i,j].normal = vetor_normal(self.rect.p4, self.rect.p1)
            #direita
            elif self.matriz[i,j-1].is_rect and (not self.matriz[i-1,j].is_rect):
              self.matriz[i,j].normal = vetor_normal(self.rect.p2, self.rect.p3)
            #baixo
            elif self.matriz[i-1,j].is_rect:
               self.matriz[i,j].normal = vetor_normal(self.rect.p3, self.rect.p4)

  
  def plot_corrente(self):
    contour_levels = np.linspace(-120, 120, 40)
    contour_levels = contour_levels[contour_levels != 0]
    plt.contour(self.corrente,cmap='cool',levels=contour_levels)
    plt.colorbar(label='Psi')
    plt.title('Corrente')
    plt.xlabel('Eixo X')
    plt.ylabel('Eixo Y')
    plt.gca().invert_yaxis()
    plt.show()

  def plot_velocidade(self):
    V_x = np.array([[ponto.ux for ponto in linha] for linha in self.matriz])
    V_y = np.array([[ponto.uy for ponto in linha] for linha in self.matriz])

    # Criar uma grade de pontos
    x = np.arange(0, V_x.shape[1])
    y = np.arange(0, V_y.shape[0])
    X, Y = np.meshgrid(x, y)

    # Plotar o gráfico de vetores
    plt.figure()
    plt.quiver(X, Y, V_x, V_y)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Campo de Velocidades')
    magnitude = np.sqrt(V_x**2 + V_y**2)
    plt.imshow(magnitude, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', interpolation='bilinear', alpha=0.8)
    plt.colorbar(label='Magnitude')
    plt.gca().invert_yaxis()
    plt.show()

  def plot_pressure(self):
    p = np.array([[ponto.pressure for ponto in linha] for linha in self.matriz])
    x = np.arange(0, p.shape[1])
    y = np.arange(0, p.shape[0])
    X, Y = np.meshgrid(x, y)

    # Encontrar a posição do menor valor de pressão
    min_pos = np.unravel_index(np.argmin(p), p.shape)
    min_pressure = p[min_pos]

    # Plotar o gráfico de contorno
    plt.contourf(X, Y, p, cmap='viridis')
    plt.colorbar(label='Pressão')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Distribuição de Pressão')

    # Inverter o eixo y
    plt.gca().invert_yaxis()

    # Destacar o ponto com o menor valor de pressão
    plt.scatter(min_pos[1], min_pos[0], color='red', s=100, edgecolor='black', label=f'Menor Pressão: {min_pressure:.2f}')
    plt.legend()

    # Mostrar o gráfico
    plt.show()

  def plot_mask(self):
    teste = np.array([[ponto.is_rect for ponto in linha] for linha in self.matriz])
    plt.imshow(teste, cmap='gray', interpolation='nearest')
    plt.title('Máscara do retangulo')
    plt.show()
  
  def plot_normal(self):
    N_x = np.array([[ponto.normal[0] for ponto in linha] for linha in self.matriz])
    N_y = np.array([[ponto.normal[1] for ponto in linha] for linha in self.matriz])

    # Criar uma grade de pontos
    x = np.arange(0, N_x.shape[1])
    y = np.arange(0, N_y.shape[0])
    X, Y = np.meshgrid(x, y)

    # Plotar o gráfico de vetores
    plt.figure()
    plt.quiver(X, Y, N_x, N_y,  scale=20)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Vetor Normal')
    plt.gca().invert_yaxis()
    plt.show()
  
  def get_arrasto(self):
    arrasto = np.sum(np.array([[ponto.pressure*ponto.normal[0] for ponto in linha] for linha in self.matriz]))/self.rect.perimeter
    return arrasto
  
  def get_sust(self):
    sust = np.sum(np.array([[ponto.pressure*ponto.normal[1] for ponto in linha] for linha in self.matriz]))/self.rect.perimeter
    return sust