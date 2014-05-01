#pragma once
#include <iostream>
#include <string>

using namespace std;

class Matrix{
private:
	int m;
		//Numero de Filas
	int n;
		//Numero de Columnas
	double **e;
		//Datos de la matriz
public:
	Matrix();
		//Constructor por default
	Matrix(int m, int n);
		//Constructor que define numero de filas y columnas
	Matrix(int order);
		//Constructor que genera una matriz cuadrada de orden "order"
	Matrix(int m, int n, double **e);
		//Constructor que genera una matriz de orden m x n con los datos e
	Matrix(int order, double **e);
		//Constructor que genera una matriz cuadrada de orden "order" con los datos e
	Matrix(Matrix *b);
		//Constructor que copia una matriz
	Matrix(string filename);
		//Constructor que genera una matriz a partir de una archivo
	~Matrix();
		//Destructor
	double getElement(int i, int j);
		//Regresa el elemento en la posición i,j
	int getRows();
		//Regresa la cantidad de filas
	int getColumns();
		//Regresa la cantidad de columnas
	void setElement(int i, int j, double value);
		//Inserta el elemento value en la posición i,j

	void setIdentity();
		//Modifica una matriz para volverla identidad
	void setNull();
		//Modifica la matriz para convertirla en nula
	//bool resolverPorGaussJordan(int n, Matrix *ab, Matrix *&sol);
	bool resolverPorGaussJordan(Matrix *&sol);
		//Resuelve por Gauss-Jordan una matriz regresandola a través de una referencia de memoria
	double * resolverPorGaussJordan(int n, Matrix * ampliada);
		//Resuelve por Gauss-Jordan y regresa los coeficientes en un arreglo de doubles
	Matrix * sum(Matrix *b);
		//Devuelve el resultado de la matriz del objeto más la mátriz b
	void print();
		//Imprime la matriz
	bool isIdentity(double r);
		//Regresa si la matriz es identidad
	Matrix * getTranspose();
		//Regresa la traspuesta de una matriz
private:
	double ** createMatrix(int m, int n);
		//Metodo que nos permite crear una matriz de tamaño m, n
};
