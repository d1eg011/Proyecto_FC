# Proyecto_FC

Proyecto final de Física Computacional: modelo de Ising Clásico bidimensional en una grilla de NxN espines (Monte Carlo)
Grupo: 04 

En la carpeta del src, se encuetran la implemetación del modelo basada en OOP en python en el archivo de `ising.py`, y las implemetaciones de las clases en C++ en serie y paralelo, para la aplicación del modelo y la creación de la simulación.

En la carpeta de Plots, se encuentran las gráficas obtenidas apartir de los códigos en C++.

El jupyter-notebook `test.ipynb` contiene las gráficas del código de python.

El jupyter'notebook `proyecto?fc.ipynb` contiene el informe del proyecto con sus resultados y análisis.

Por último, la carpeta de `ising_doc` tiene los documentos utilizados para generar la documentación del código de python, que se encuentra en: `https://d1eg011.github.io/Proyecto_FC/` 

## Instrucciones y dependencias del código de Ising en C++ 
Para poder utilizar el modulo gnuplot-iostream para hacer plots en C++ es necesario tener instalada la libreria boost de C++, se puede install con el siguiente comando 

`sudo apt-get install libboost-all-dev`

luego para correr el código debe hacerse de la siguiente manera 

* `g++ -c *.cpp`
* `g++ -Wall run_ising.o spin.o -o filename.x -lboost_iostreams -lboost_system -lboost_filesystem -lboost_program_options`

El código run_ising puede ser reemplazado por otros archivos main.
