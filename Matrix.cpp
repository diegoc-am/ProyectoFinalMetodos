#include <iostream>
#include <cmath>
#include "Matrix.h"
#include <fstream>
#include <string>

using namespace std;

Matrix::Matrix() //matriz nula [0]
{
	this->m = 1;
	this->n = 1;
	this->e = this->createMatrix(1, 1);
	this->e[0][0] = 0.0;
}
Matrix::Matrix(int m, int n) //matriz ordnen m x n
{
	this->m = m;
	this->n = n;
	this->e = this->createMatrix(m, n);
}
Matrix::Matrix(int order) // matriz orden "order"
{
	this->m = order;
	this->n = order;
	this->e = this->createMatrix(m, n);
}
Matrix::Matrix(int m, int n, double **e) //matriz orden m x n y elementos por matriz dinamica e
{
	this->m = m;
	this->n = n;
	this->e = this->createMatrix(m, n);
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			this->e[i][j] = e[i][j];
		}
	}
}
Matrix::Matrix(string fileName){
	ifstream ifs;
	ifs.open(fileName.c_str());
	int mm,nn;
	ifs >> mm;
	ifs >> nn;
	this->m = mm;
	this->n = nn;
	this->e = this->createMatrix(mm, nn);
	int temp;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ifs >> temp;
			this->e[i][j] = temp;
		}
	}
}
Matrix::Matrix(int order, double **e) //matriz cudadrada y elemento de matriz dinamica
{
	this->m = order;
	this->n = order;
	this->e = this->createMatrix(order, order);
	for(int i = 0; i < order; i++)
	{
		for(int j = 0; j < order; j++)
		{
			this->e[i][j] = e[i][j];
		}
	}
}
Matrix::Matrix(Matrix *b) //matriz copia
{
	this->m = b->getRows();
	this->n = b->getColumns();
	this->e = this->createMatrix(m, n);
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			this->e[i][j] = b->getElement(i, j);
		}
	}
}
Matrix::~Matrix() //destructor matriz
{
	for(int i = 0; i < m; i++)
	{
		delete [] e[i];
	}
	delete []e;
}
double Matrix::getElement(int i, int j) //metodos analizadores o get
{
	return this->e[i][j];
}
int Matrix::getRows()
{
	return this->m;
}
int Matrix::getColumns()
{
	return this->n;
}
void Matrix::setElement(int i, int j, double value) //metodos modificadores o set
{
	this->e[i][j] = value;
}
//pre: m == n
void Matrix::setIdentity()
{
	for(int i = 0; i < this->m; i++)
	{
		for(int j = 0; j < this->n; j++)
		{
			if(i == j)
			{
				this->e[i][j] = 1.0;
			}
			else
			{
				this->e[i][j] = 0.0;
			}
		}
	}
}
void Matrix::setNull()
{
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			this->e[i][j] = 0.0;
		}
	}
}
//pre: mismo orden
Matrix * Matrix::sum(Matrix *b) //suna de matrices
{
	int m = b->getRows();
	int n = b->getColumns();
	Matrix *c = new Matrix(m, n);
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			c->setElement(i, j, this->getElement(i, j) + b->getElement(i, j));
		}
	}
	return c;
}
void Matrix::print() //imprime matriz en salida std
{
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			cout << this->e[i][j] << " ";
		}
		cout << endl;
	}
}
double ** Matrix::createMatrix(int m, int n)
{
	double **a = new double *[m];
	for(int i = 0; i < m; i++)
	{
		a[i] = new double[n];
	}
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			a[i][j] = 0.0;
		}
	}
	return a;
}

bool Matrix::isIdentity(double err)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i != j){
				if (fabs(this->e[i][j]) > err){
					return false;
				}
			}
			else{
				if ((this->e[i][j] < 1.0 - err) && (this->e[i][j] > 1.0 + err)){
					return false;
				}
			}
		}
	}
	return true;
}

Matrix * Matrix::getTranspose(){

	int m = this->getColumns(); //numero de filas de matriz transpuesta
	int n = this->getRows(); // numero de columnas de matriz transpuesta
	Matrix *t = new Matrix(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			t->setElement(i, j, this->getElement(j, i));
		}
	}
	return t;
}

//M�todo invocado para Matrices ampliadas recibe la matriz sobre la que se guardar� el arreglo de soluci�n.
//La operaci�n se realiza sobre la matriz que contiene el objeto Matrix
bool Matrix::resolverPorGaussJordan(Matrix *&sol){
	double CERO = 0.000001;
	bool haysol = false;
	double piv;
	double det=1;
	Matrix *t = new Matrix(this);
	//cout << "Usted ha introducido la siguiente matriz ampliada: " << endl;
	//t->print();
	for (int i = 0; i < n-1; i++)
	{
		piv = t->getElement(i,i);
		det = det*piv;
		if (piv == 0 || (abs(piv)) <= CERO){
			for (int k = 0; k < m; k++){
				if (t->e[k][i] != 0 || (abs(piv))>CERO){
					double *temp = t->e[k];
					t->e[k] = t->e[i];
					t->e[i] = temp;
				}
			//cout << t->e[k][i];
			}
		}
		else if (piv != 1){
			for (int k = 0; k < n; k++){
					t->e[i][k] = t->e[i][k]/piv;
			}
		}
			for (int j = 0; j < m; j++)
			{
				if (abs(t->e[j][i]) > CERO && j!=i){
					double reciproco =t->e[j][i] * (-1);
					for (int k = 0; k < n; k++){
						t->e[j][k] = t->e[i][k]*reciproco + t->e[j][k];
					}
				}
			}
			
	}
	//cout << "Tabla Gauss Jordan resuelta" << endl;
	//t->print();
	sol = new Matrix(1, m);
	for (int i = 0; i < m; i++){
		sol->e[0][i] = t->e[i][n-1];
	}
	//cout << "Fin" << endl;
	return haysol;
}

double  * Matrix::resolverPorGaussJordan(int n, Matrix * ampliada){
	double * raices = new double[n];
	int i, j, k;
	double factor;
	for (k = 0; k < n; k++) {
		factor = ampliada->e[k][k];
		cout << "Pivote: " << factor << endl;
		for (j = k; j < n + 1; j++) {
			ampliada->e[k][j] = ampliada->e[k][j] / factor;
		}
		for (i = 0; i < ampliada->n -1 ; i++){
			if (i != k) {
				factor = ampliada->e[i][k];
				for (j = 0; j < n + 1; j++) {
					ampliada->e[i][j] = ampliada->e[i][j] - factor * ampliada->e[k][j];
				}
			}
		}
		ampliada->print();
	}
	cout << "Raices: ";
	for (int ra = 0; ra < n; ra++){
		raices[ra] = ampliada->e[ra][n];
		cout << raices[ra] << ", ";
	}
	cout << endl;
	return raices;
}
