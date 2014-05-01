#define CERO 0.00001
#include <iostream>
#include "Matrix.h"
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

class minimosCuadrados{
private:
		double * equis; //Arreglo de los valores de X
		double * ye; //Arreglo de los valores de Y
		double Sx; //Sumatoria de los valores de X
		double Sy; //Sumatoria de los valores de Y
		Matrix * coeficientesNoResueltos; /*Matriz donde se guarda la información
		que será utilizada en Gauss-Jordan para obtener resultado*/
		double * coeficientes; //Coeficientes del polinomio resultante
		int tama; //Tamaño de la cantidad de puntos que se tiene
		int tamaCoefs; //Tamaño de el arreglo que satisface el polinomio
public:
		minimosCuadrados(); //Constructor por default de la clase, pide un archivo e inicializa el menú
		minimosCuadrados(int tam); //Constructor que recibe un tamaño para pedir los puntos
		minimosCuadrados(string archivo); //Constructor que genera las tablas desde un archivo
		~minimosCuadrados(); //Destructor
		void generarSumatorias(); //Genera las sumatorias de X & Y
		void resolverMinimosCuadrados(int ord); //Resuelve por mínimos cuadrados según el orden específicado
		void gaussJordan(); //Resuelve por Gauss-Jordan la matriz para obtener los coeficientes
		double sXOrden(int orden); //Regresa la sumatoria de X^orden
		double sYOrden(int orden); //Regresa la sumatoria de Y^orden
		double sXYOrden(int ordenX, int ordenY); //Regresa la sumatoria de X^ordenX * Y^ordenY
		double sXYOrden(int orden); //Regresa la sumatoria de X^orden * Y^orden
		double sumatoria(double * a, int n); //Regresa la sumatoria de un arreglo de tamaño n
		void printCoeficientes(); //Imprime el polinomio resultante
		double hConstanteAll();  //Checa si las diferencias son finitas o no
		int ordenOptimo(); //Revisa si el orden es el optimo o no
		void print(); //Imprime los datos
		void usar(); //Menú para el usuario
		void readFile();//Permite leer un archivo

};
