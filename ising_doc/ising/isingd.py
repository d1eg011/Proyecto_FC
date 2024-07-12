"""Los objetos y las funciones presentadas permiten llevar a cabo el modelo de Ising mediante el algoritmo de Metrópolis y el método de Montecarlo para un material sometido a distintas temperaturas.

Se incluye

- `class Spins`: Los objetos de esta clase están asociados a grillas cuadradas con espines en cada uno de sus puntos. Los métodos permiten calcular propiedades como la energía y la magnetización, además de someter un spin al algoritmo de metrópolis.

- `simulate_Ising`: Permite someter a una grilla al modelo de Ising aplicando múltiples veces el método de Montecarlo para un rango de temperaturas especificado.
"""

import numpy as np

class Spins:

    """
    Los objetos de esta clase corresponden a grillas cuadradas con espines hacia arriba o hacia abajo en cada uno de sus puntos.

    Attributes:
        lattice (list): Lista de listas que define la grilla de la configuración.
        size (int): Tamaño de la grilla cuadrada.

    Methods:
        calculate_energy(self): Calcula la energía de la grilla.
        calculate_magnetization(self): Calcula la magnetización de la grilla.
        update(self, i, j, T, rand): Actualiza la grilla con base en el algoritmo de Metrópolis.

    """
    
    def __init__(self, init_state, size):
        """Constructor de las grillas de espines. Cuenta con los atributos lattice (arreglo correspondiente a la grilla) y size (tamaño de la grilla).
        
        Examples:
        
        Se crea una grilla de tamaño 2 con espines hacia abajo y se obtiene su tamaño.

            >>> spins_t = Spins(-1, 2)
            >>> spins_t.size
            2

        Args:
            init_state (int o float): Valor que indica el estado inicial de la grilla (-1: espines hacia abajo, 0: espines aleatorios, cualquier otro: espines hacia arriba).
            size (int): Tamaño de la grilla cuadrada.
                     
        """

        if init_state == 0:
            self.lattice = np.random.choice([-1, 1], size=(size, size))
        elif init_state == -1:
            self.lattice = np.full((size, size), -1)
        else:
            self.lattice = np.ones((size, size))
        self.size = size
        
    def calculate_energy(self):
    
    	"""Cálculo de la energı́a respecto a los vecinos más cercanos en la grilla.
		
        Examples:

        Se crea una grilla de tamaño 2 con espines hacia abajo y se calcula su energía.
        
	    >>> spins_t = Spins(-1, 2)
            >>> spins_t.calculate_energy()
	    -8.0
	
	Returns:
	    (float): Suma de interacción con los cuatro espines más cercanos por el parámetro energético de interacción -J.
		
	"""
    	
    	neighbors_sum = (np.roll(self.lattice, 1, axis=1)
    	    + np.roll(self.lattice, -1, axis=1) + np.roll(self.lattice, 1, axis=0) + np.roll(self.lattice, -1, axis=0)) 
    	
    	return -0.5 * (self.lattice * neighbors_sum).sum()
    
    def calculate_magnetization(self):
        
        """Cálculo de la magnetización según los espines de una grilla.
        Examples:

        Se crea una grilla de tamaño 2 y se calcula su magnetización.

            >>> spins_t = Spins(-1, 2)
            >>> spins_t.calculate_magnetization()
            -4

        Returns:
            (int): Magnetización de la configuración.
        """    
        
        return np.sum(self.lattice)       
        
    def update(self, i, j, T, rand):
        """Cálculo del cambio de energı́a entre los espines más cercanos en la grilla y cálculo del cambio de la magnetizació        n con base en el algoritmo de metrópolis.
		
        Examples:

        Se crea una grilla de tamaño 4 con espines hacia abajo y se actualiza el spin en la posición (2, 2) bajo una tempera        tura de 2 y una probabilidad aleatoriqa -0.1.

	    >>>	spin_t = Spins(-1, 4)
            >>> spin_t.update(2, 2, 2, 1)
            (0, 0)

	Args:		
            i (int): Fila (posición en x).
			j (int): Columna (posición en y).
            T (int or float): Temperatura.
            rand (int or float): Número aleatorio (pseudo-aleatorio) con el que se ejecuta el algoritmo de metrópolis. Este             último contiene comparaciones del argumento con el logaritmo natural de la probabilidad dada por la distribución de Boltzmann, por lo que generalmente se usan valores negativos.

	Returns:
		dE (int): Si el valor aleatorio dado es menor a la probabilidad dada por la distribución de Boltzmann, da el valor del cambio de la energı́a respecto a los espines más cercanos. Si no, retorna 0.
 		dM (int): Si el valor aleatorio dado es menor a la probabilidad dada por la distribución de Boltzmann, da el valor del cambio de la magnetización respecto a los espines más cercanos. Si no, retorna 0.
            
	"""
	
        dE = (
        2
        * self.lattice[i, j]
        * (
            self.lattice[(i + 1) % self.size, j]
            + self.lattice[(i - 1) % self.size, j]
            + self.lattice[i, (j + 1) % self.size]
            + self.lattice[i, (j - 1) % self.size]
        )
        )
        if rand < - dE / (T):
            self.lattice[i, j] *= -1
            dM = 2 * self.lattice[i, j]
            return dE, dM

        return 0, 0     
	
        
def simulate_Ising(Ti, Tf, steps, size, mcsteps, thermsteps, init_state=1, J=1, kb=1):
    """Ejecución del modelo de Ising para un intervalo de temperaturas especificado a partir de una configuración inicial y la generación de números psuedo-aleatorios (método de Monte Carlo). Primero se lleva a cabo la termalización y luego se guardan los valores según los steps de Montecarlo.
    
    Examples:
    
    Se somete una grilla de tamaño 3 con espines hacia arriba al modelo de Ising con 5 steps de Montecarlo, 2 steps de terma    lización y un intervalo de temperaturas que va de 0.5 a 5 con 3 steps. El resultado varía con la semilla, la cual es tomada de la hora del sistema.

        >>> simulate_Ising(Ti=0.5, Tf=5, steps=3, size=3, mcsteps=5, thermsteps=2, init_state=1)
        (array([[-18., -18., -18., -18., -18.], [-18., -18., -18., -18., -18.], [-18., -18., -18., -18., -18.]]),
 		array([[9., 9., 9., 9., 9.], [9., 9., 9., 9., 9.], [9., 9., 9., 9., 9.]]),
 		array([0.5 , 2.75, 5.  ]))

    Args:
        Ti (int or float): Temperatura inicial.
        Tf (int or float): Temperatura final.
        steps (int): Cantidad de temperaturas en el intervalo indicado.
        size (int): Tamaño de la grilla cuadrada.
        mcsteps (int): Pasos del método de Montecarlo.
        thermsteps (int): Pasos de termalización.
        init_state (int or float): Valor con el que se inicializa la grilla de espines según el constructor.
        J (int or float, optional): Propiedad magnética del material (fijado como 1 por defecto). 
        kb(int or float, optional): Constante de Boltzmann (fijado como 1 por defecto).

    Returns: 
        energies (list): Lista de listas con las energías de la configuración para cada temperatura.
        magnetizations (list): Lista con las magnetizaciones de la configuración para cada temperatura.
        temperatures (list): Lista con las temperaturas empleadas.

    """

    spins = Spins(init_state, size)
    E = spins.calculate_energy()
    M = spins.calculate_magnetization()

    temperatures = np.linspace(Ti, Tf, steps)
    energies = np.zeros((steps, mcsteps))
    magnetizations = np.zeros((steps, mcsteps))

    rand_pos = np.random.randint(size, size=((steps, mcsteps + thermsteps, 2)))
    rands = np.log(np.random.uniform(size=(steps, mcsteps + thermsteps)))

    for k in range(steps):

        for l in range(mcsteps + thermsteps):
            i, j = rand_pos[k, l]

            dE, dM = spins.update(i, j, temperatures[k], rands[k, l])
            E += dE
            M += dM

            if l >= thermsteps:
                energies[k, l-thermsteps] = E
                magnetizations[k, l-thermsteps] = M


    return energies, magnetizations, temperatures
