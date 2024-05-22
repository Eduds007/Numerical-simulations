import numpy as np
import matplotlib.pyplot as plt
from utils import *



class Point:
    """
    Representa um ponto em uma malha de simulação de aerofólio.

    Atributos:
        x (float): Coordenada x do ponto.
        y (float): Coordenada y do ponto.
        is_rect (bool): Indica se o ponto está dentro da região retangular.
        corrente (float): Valor da corrente no ponto.
        ux (float): Componente x da velocidade.
        uy (float): Componente y da velocidade.
        pressure (float): Pressão no ponto.
        normal (list): Vetor normal no ponto.
    """
        
    def __init__(self, i, j, u_inf, delta, uy, corrente, p_init) -> None:
        
        """
        Inicializa um ponto na malha de simulação.

        Args:
            i (int): Índice i do ponto na malha.
            j (int): Índice j do ponto na malha.
            u_inf (float): Velocidade do fluxo livre.
            delta (float): Distância entre os pontos na malha.
            uy (float): Componente y da velocidade.
            corrente (float): Valor da corrente no ponto.
            p_init (float): Pressão inicial no ponto.
        """

        self.x = i*delta
        self.y = j*delta
        self.is_rect = False
        self.corrente = corrente
        self.ux = u_inf
        self.uy =uy
        self.pressure = p_init
        self.normal = [0,0]
        self.is_a = False
        self.is_b = False


class Retangle(Point):
  """
  Representa um retângulo definido por quatro pontos, com rotação aplicada.

  Atributos:
      p1 (Point): Primeiro ponto do retângulo.
      p2 (Point): Segundo ponto do retângulo.
      p3 (Point): Terceiro ponto do retângulo.
      p4 (Point): Quarto ponto do retângulo.
      perimeter (float): Perímetro do retângulo.
      angle (float): Ângulo de rotação do retângulo.
  """
  def __init__(self, x1, y1,angle, u_inf, delta, uy, corrente, p_init, h, w ):
    """
        Inicializa um retângulo com quatro pontos e aplica rotação.

        Args:
            x1 (float): Coordenada x do ponto inferior esquerdo do retângulo.
            y1 (float): Coordenada y do ponto inferior esquerdo do retângulo.
            angle (float): Ângulo de rotação do retângulo.
            u_inf (float): Velocidade do fluxo livre.
            delta (float): Distância entre os pontos na malha.
            uy (float): Componente y da velocidade.
            corrente (float): Valor da corrente no ponto.
            p_init (float): Pressão inicial no ponto.
            h (float): Altura do retângulo.
            w (float): Largura do retângulo.
    """
    
    delta = 1
    self.p2 = rotate(Point(x1,y1+w,u_inf, delta, uy, corrente, p_init), angle)   
    self.p3 = rotate(Point(x1+h,y1+w,u_inf, delta, uy, corrente, p_init), angle)
    self.p1 = rotate(Point(x1,y1,u_inf, delta, uy, corrente, p_init), angle)
    self.p4 = rotate(Point(x1+h,y1,u_inf, delta, uy, corrente, p_init), angle)

    self.perimeter = 2*h + 2*w
    self.angle = angle

  def is_inside(self, ponto):
    """
        Verifica se um ponto está dentro do retângulo.

        Args:
            ponto (Point): O ponto a ser verificado.

        Returns:
            bool: True se o ponto estiver dentro do retângulo, False caso contrário.
    """
    def is_left(p1, p2, p3):
        return (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y)

    def is_in_triangle(p, p1, p2, p3):
        b1 = is_left(p, p1, p2) <= 0.0
        b2 = is_left(p, p2, p3) <= 0.0
        b3 = is_left(p, p3, p1) <= 0.0
        return (b1 == b2) and (b2 == b3)

    return is_in_triangle(ponto, self.p1, self.p2, self.p3) or is_in_triangle(ponto, self.p1, self.p3, self.p4)

class Mesh():
  """
    Representa uma malha de simulação para um aerofólio.

    Atributos:
        linhas (int): Número de linhas na malha.
        colunas (int): Número de colunas na malha.
        rect (Retangle): Objeto Retangle representando o retângulo na malha.
        tolerance (float): Tolerância para o critério de parada.
        delta (float): Distância entre os pontos na malha.
        u_inf (float): Velocidade do fluxo livre.
        lamb (float): Parâmetro de relaxação.
        iter (int): Número de iterações para atualização da função corrente.
        ux (float): Componente x da velocidade.
        uy (float): Componente y da velocidade.
        corrente_inicial (float): Valor inicial da função corrente.
        p_init (float): Pressão inicial.
        rho (float): Densidade do fluido.
        g (float): Aceleração devido à gravidade.
        matriz (np.ndarray): Matriz de objetos Point representando a malha.
        corrente (np.ndarray): Matriz com os valores da função corrente.
  """
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

    """Atualiza a máscara que indica se os pontos estão dentro do retângulo."""

    shape = self.matriz.shape
    for i in range(1,shape[0]-1):
      for j in range(1,shape[1]-1):
          if self.rect.is_inside(self.matriz[i, j]):
            self.matriz[i,j].is_rect = True

  def update_velocidade(self):

    """Atualiza as componentes da velocidade em cada ponto da malha."""

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

    """Atualiza a pressão em cada ponto da malha."""

    shape = self.matriz.shape
    for i in range(1,shape[0]-1):

        for j in range(1,shape[1]-1):

          if self.rect.is_inside(self.matriz[i, j]):

            self.matriz[i, j].pressure = self.p_init
          else:
            self.matriz[i, j].pressure = self.rho*self.g*(self.matriz[0, 0].pressure/(self.rho*self.g) + (self.matriz[0, 0].ux)**2/(2*self.g) - (self.matriz[i, j].ux**2 + self.matriz[i, j].uy**2)/(2*self.g))


  def update_corrente(self):

    """
    Atualiza a função corrente em cada ponto da malha 
    utilizando o método de sobrerelaxação com as equações
    discretas para cada caso do sistema
    """

    shape = self.matriz.shape


    for k in range(self.iter):
      aux = self.matriz.copy()
      erro = []
      for i in range(0,shape[0]):

        for j in range(0,shape[1]):
          if self.rect.is_inside(self.matriz[i, j]):
            self.matriz[i,j].corrente = 0

          elif (not (i==0 or i==shape[0]-1)) and (j==0): #Entrada

            novo = self.lamb*(2*aux[i,j+1].corrente +  aux[i+1,j].corrente +  aux[i-1,j].corrente)/4 + (1-self.lamb)*aux[i,j].corrente
            
            self.matriz[i,j].corrente = novo
            

          elif (not(i==0 or i==shape[0]-1)) and (j==shape[1]-1): #Saída

            novo = self.lamb*(2*aux[i,j-1].corrente +  aux[i+1,j].corrente +  aux[i-1,j].corrente)/4 + (1-self.lamb)*aux[i,j].corrente
            try:
              erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
            except ZeroDivisionError:
              pass
            self.matriz[i,j].corrente = novo
            pass

          elif (i==0) and (not(j==0 or j==shape[1]-1)): #Superior

            novo = self.lamb*(aux[i,j+1].corrente +  aux[i,j-1].corrente +  2*aux[i+1,j].corrente + 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
            try:
              erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
            except ZeroDivisionError:
              pass
            self.matriz[i,j].corrente = novo

          elif (i==shape[0]-1) and (not(j==0 or j==shape[1]-1)): #Inferior

            novo = self.lamb*(aux[i,j+1].corrente +  aux[i,j-1].corrente +  2*aux[i-1,j].corrente - 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
            try:
              erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
            except ZeroDivisionError:
              pass
            self.matriz[i,j].corrente = novo

          elif (i==0 or i==shape[0]-1) and (j==0 or j==shape[1]-1): #Canto

            if (i==0 and j==0):

              novo = self.lamb*(2*aux[i+1,j].corrente +  2*aux[i,j+1].corrente + 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo

            elif (i==shape[0]-1 and j==0): #canto inferior esquerdo

              novo = self.lamb*(2*aux[i-1,j].corrente +  2*aux[i,j+1].corrente - 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo

            elif (i==0 and j==shape[1]-1): #canto superior direito

              novo = self.lamb*(2*aux[i+1,j].corrente +  2*aux[i,j-1].corrente + 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo
            elif (i==shape[0]-1 and j==shape[1]-1): #canto inferior direito

              novo = self.lamb*(2*aux[i-1,j].corrente +  2*aux[i,j-1].corrente - 2*self.delta*self.u_inf)/4 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo


            pass
          elif self.matriz[i,j].is_a: #canto A
              novo = self.lamb*(aux[i-1,j].corrente +  aux[i,j+1].corrente)/2 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo
          elif self.matriz[i,j].is_b: #canto B
              novo = self.lamb*(aux[i+1,j].corrente +  aux[i,j+1].corrente)/2 + (1-self.lamb)*aux[i,j].corrente
              try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
              except ZeroDivisionError:
                pass
              self.matriz[i,j].corrente = novo



          else:

            novo = self.lamb*(aux[i+1,j].corrente +  aux[i-1,j].corrente +  aux[i,j+1].corrente +  aux[i,j-1].corrente)/4 + (1-self.lamb)*aux[i,j].corrente

            try:
                erro.append(np.absolute((novo-aux[i,j].corrente)/novo))
            except ZeroDivisionError:
                pass
            self.matriz[i,j].corrente = novo

      self.corrente =  np.array([[ponto.corrente for ponto in linha] for linha in self.matriz])

      if np.max(erro) <= self.tolerance:
        break

  def update_normal(self):

    """Atualiza os vetores normais na malha."""

    shape = self.matriz.shape
    if self.rect.angle == 15:
      for i in range(1,shape[0]-1):
        for j in range(1,shape[1]-1):
          if not self.matriz[i,j].is_rect:
            #cima
            if self.matriz[i+1,j].is_rect:

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

              self.matriz[i,j].normal = vetor_normal(self.rect.p4, self.rect.p1)
            #direita
            elif self.matriz[i,j-1].is_rect and (not self.matriz[i-1,j].is_rect):
              self.matriz[i,j].normal = vetor_normal(self.rect.p2, self.rect.p3)
            #baixo
            elif self.matriz[i-1,j].is_rect:
               self.matriz[i,j].normal = vetor_normal(self.rect.p3, self.rect.p4)

  
  def plot_corrente(self):

    """Realiza o plot do gráfico da corrente"""

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

    """Realiza o plot do gráfico da velocidade"""

    V_x = np.array([[ponto.ux for ponto in linha] for linha in self.matriz])
    V_y = np.array([[ponto.uy for ponto in linha] for linha in self.matriz])

    x = np.arange(0, V_x.shape[1])
    y = np.arange(0, V_y.shape[0])
    X, Y = np.meshgrid(x, y)

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


    """Realiza o plot do gráfico da pressão"""

    p = np.array([[ponto.pressure for ponto in linha] for linha in self.matriz])
    x = np.arange(0, p.shape[1])
    y = np.arange(0, p.shape[0])
    X, Y = np.meshgrid(x, y)


    min_pos = np.unravel_index(np.argmin(p), p.shape)
    min_pressure = p[min_pos]

    plt.contourf(X, Y, p, cmap='viridis')
    plt.colorbar(label='Pressão')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Distribuição de Pressão')


    plt.gca().invert_yaxis()

    # Destacar o ponto com o menor valor de pressão
    plt.scatter(min_pos[1], min_pos[0], color='red', s=100, edgecolor='black', label=f'Menor Pressão: {min_pressure:.2f}')
    plt.legend()


    plt.show()

  def plot_mask(self):

    """Realiza o plot do gráfico do posicionamento do Aerofólio"""

    teste = np.array([[ponto.is_rect for ponto in linha] for linha in self.matriz])
    plt.imshow(teste, cmap='gray', interpolation='nearest')
    plt.title('Máscara do retangulo')
    plt.show()
  
  def plot_normal(self):
    N_x = np.array([[ponto.normal[0] for ponto in linha] for linha in self.matriz])
    N_y = np.array([[ponto.normal[1] for ponto in linha] for linha in self.matriz])

    x = np.arange(0, N_x.shape[1])
    y = np.arange(0, N_y.shape[0])
    X, Y = np.meshgrid(x, y)

    plt.figure()
    plt.quiver(X, Y, N_x, N_y,  scale=20)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Vetor Normal')
    plt.gca().invert_yaxis()
    plt.show()
  
  def get_arrasto(self):

    """Obter força de arrasto"""

    arrasto = -np.sum(np.array([[ponto.pressure*ponto.normal[0] for ponto in linha] for linha in self.matriz]))/self.rect.perimeter
    return arrasto
  
  def get_sust(self):

    """Obter força de sustentação"""

    sust = np.sum(np.array([[ponto.pressure*ponto.normal[1] for ponto in linha] for linha in self.matriz]))/self.rect.perimeter
    return sust