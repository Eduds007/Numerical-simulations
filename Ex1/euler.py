import numpy as np

def euler_leg( h, f0):

  #constantes
  m1 = 0.5
  m2 = 1
  l1 = 1
  l2 = 1.5
  g = 9.8

  #conversão para variveis da função

  teta1 = f0[0]
  teta2 = f0[1]
  teta11 = f0[2]
  teta21 = f0[3]

  
  num_g1 = (-np.sin(teta1-teta2)*(m2*l2*(teta21**2)+m2*l1*(teta11**2)*np.cos(teta1-teta2))-g*((m1 +m2)*np.sin(teta1) - m2*(np.sin(teta2)*np.cos(teta1-teta2))))
  den_g1 = l1*(m1+m2*((np.sin(teta1-teta2))**2))
  g1 = num_g1/den_g1


  num_g2 = (np.sin(teta1-teta2)*((m1+m2)*l1*(teta11**2)+m2*l2*(teta21**2)*np.cos(teta1-teta2))+g*((m1 +m2)*(np.sin(teta1)*np.cos(teta1-teta2)-np.sin(teta2))))
  den_g2 = l2*(m1+m2*((np.sin(teta1-teta2))**2))
  g2 = num_g2/den_g2


  fxy = np.array([
    teta11,
    teta21,
    g1,
    g2,
    0,
    0
  ])

  f0[4] = g1
  f0[5] = g2

  f_out = f0 + h*fxy

  return f_out

def euler_steps( duration, h, f0 ):
  
  time = 0
  steps = int(duration/h)

  response = np.zeros((steps, 7))
  response[0,1:] = f0
  response[0,0] = 0

  for i in range(1,steps):

    f0 = euler_leg(h, f0)

    time += 1
    response[i,0] = time*h
    response[i, 1:] = f0

  return response

