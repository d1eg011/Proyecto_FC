# Información de las funciones 

## Calculate energy

Para esta función, utiliza la funcionalidad de `numpy.roll`, para sumar sobre los vecinos más cercanos de la siguiente forma:

- Se le dan los argumentos del array de espines, se le indica que el cambio de posiciones en el arreglo, en este caso va a ser solamente en una unidad (1 o -1), respecto a uno de los dos ejes posibles 0 y 1.
- Y devuleven el nuevo arreglo de espines tras aplicar el cambio de posiciones.

Esto se ejemplifica con:

		>>> espines = np.array([[1,-1,1],[-1,1,1]])
		>>>	np.roll(espines, 1, axis= 1)
		array([[1,1,-1],[1,-1,1]])
		>>> np.roll(espines,-1,axis = 0)
		array([[-1,1,1],[1,-1,1]])

- Se aplica para las cuatro posiciones de los espines cercanos (arriba, abajo,izquierda y derecha) y se suman.
- Y al final se suman sobre todas las contribuciones de cada espin y se multiplica por el parámetro de interación J.

## Calculate dE

En este caso se analiza solamnete el cambio de energı́a del espín que cambia con respecto a los cuatro espines cercanos ya mencionados, dicho cambio de energía va a estar basado en:

$$
dE = 2 \hat{S}^{z}_{i,j} \left[ \hat{S}^{z}_{i,j+1} + \hat{S}^{z}_{i,j-1} + \hat{S}^{z}_{i+1,j} + \hat{S}^{z}_{i-1,j} \right]
$$

Donde los $\hat{S}^{z}_{i,j}$ espines en una posición de la grilla ${i,j}$.
