# Proyecto_FC

Proyecto final de Física Computacional: modelo de Ising bidimensional en una grilla de NxN espines (Monte Carlo)

Grupo: 04 (?)

La carpeta "Apoyo" contiene información útil para comprender cómo funciona el modelo de Ising (ising_1.pdf, ising_2.pdf, ising_3.pdf) y una aplicación del método de Montecarlo para obtener una aproximación del valor de pi (A1.py y A2.py).

Video sobre la implementación del hamiltoniano en Python: https://youtu.be/6W881fbHGlU?si=7u56C31-E7Gd3HeS

## Instrucciones y dependencias del código de Ising en C++ 
Para poder utilizar el modulo gnuplot-iostream para hacer plots en C++ es necesario tener instalada la libreria boost de C++, se puede install con el siguiente comando 

`sudo apt-get install libboost-all-dev`

luego para correr el código debe hacerse de la siguiente manera 

* `g++ -c *.cpp`
* `g++ -Wall run_ising.o spin.o -o filename.x -lboost_iostreams -lboost_system -lboost_filesystem -lboost_program_options`

El código run_ising puede ser reemplazado por otros archivos main.
