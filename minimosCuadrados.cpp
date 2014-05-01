#define CERO 0.00001
#include <iostream>
#include "Matrix.h"
#include "minimosCuadrados.h"
#include <vector>
#include <fstream>

using namespace std;


		minimosCuadrados::minimosCuadrados(){
			this->equis = 0;
			this->ye = 0;
			this->Sx=0;
			this->Sy=0;
			this->coeficientesNoResueltos = 0;
			this->coeficientes = 0;
			this->tamaCoefs = 0;
			this->tama = 0;
			this->readFile();
			this->usar();
		}


		minimosCuadrados::minimosCuadrados(int tam){
			this->tama = tam;
			this->equis = new double[this->tama];
			this->ye = new double[this->tama];

			for (int i = 0; i < this->tama; i++){
				cout << "Dame el valor para x" << i << ": ";
				cin >> equis[i];
				cout << "Dame el valor para y" << i << ": ";
				cin >> ye[i];
			}
			this->Sx = 0;
			this->Sy = 0;
			this->generarSumatorias();
			this->tamaCoefs=0;
			this->coeficientes = new double[this->tamaCoefs];
			this->coeficientesNoResueltos = new Matrix();
		}
		minimosCuadrados::minimosCuadrados(string archivo){

			//cout << "Leyendo Archivo..." << endl;
			vector<double> x;
			vector<double> y;
			double v;
			ifstream ifs;
			ifs.open(archivo.c_str());
			while (!ifs.eof()){
				ifs >> v;
				x.push_back(v);
				ifs >> v;
				y.push_back(v);
			}
			ifs.close();
			this->tama = x.size();
			this->equis = new double[this->tama];
			this->ye = new double[this->tama];
			for (int i = 0; i < this->tama; i++){
				equis[i] = x.at(i);
				ye[i] = y.at(i);
			}
			this->Sx = 0;
			this->Sy = 0;
			this->generarSumatorias();
			this->tamaCoefs=0;
			this->coeficientes = new double[this->tamaCoefs];
			this->coeficientesNoResueltos = new Matrix();
		}

		minimosCuadrados::~minimosCuadrados(){
			delete[]equis;
			delete[]ye;
			delete[]coeficientes;
			delete coeficientesNoResueltos;
		}

		void minimosCuadrados::generarSumatorias(){
			for (int i = 0 ; i < this-> tama ; i++){
				this->Sx += this->equis[i];
				this->Sy += this->ye[i];
			}
		}

		void minimosCuadrados::resolverMinimosCuadrados(int ord){
			int m = ord + 1;
			int n = ord + 2;
			this->coeficientesNoResueltos = new Matrix(m, n);
			//this->coeficientesNoResueltos->print();
			double val = 0;
			for(int i = 0; i < m; i++){
				for (int j = 0; j < n; j++){
					if (j == n - 1){
						val = this->sXYOrden(i, 1);
					}
					else{
						val = this->sXOrden(i + j);
					}
					this->coeficientesNoResueltos->setElement(i, j, val);
				}
			}
			//this->tamaCoefs = this->coeficientesNoResueltos->getRows();
			cout << "Matriz ampliada de coeficientes no resueltos: " << endl;
			this->coeficientesNoResueltos->print();
		}

		void minimosCuadrados::gaussJordan(){
			Matrix *resultado = new Matrix();
			this->coeficientesNoResueltos->resolverPorGaussJordan(resultado);
			//resultado->print();
			this->tamaCoefs = resultado->getColumns();
			this->coeficientes = new double[tamaCoefs];
			for (int i = 0; i < this->tamaCoefs; i++){
				this->coeficientes[i] = resultado->getElement(0,i);
			}
		}


		double minimosCuadrados::sXOrden(int orden){
			double * temp = new double[this->tama];
			for (int i = 0; i < this->tama ; i++){
				temp[i] = pow(this->equis[i],orden);
			}
			double s = sumatoria(temp, this->tama);
			delete[]temp;
			return s;
		}

		double minimosCuadrados::sYOrden(int orden){
			double * temp = new double[this->tama];
			for (int i = 0; i < this->tama; i++){
				temp[i] = pow(this->ye[i], orden);
			}
			double s = sumatoria(temp, this->tama);
			delete[]temp;
			return s;
		}

		double minimosCuadrados::sXYOrden(int ordenX, int ordenY){
			double * tempX = new double[this->tama];
			double * tempY = new double[this->tama];
			double * temp = new double[this->tama];
			for (int i = 0; i< this->tama; i++){
				tempX[i] = pow(this->equis[i], ordenX);
				tempY[i] = pow(this->ye[i], ordenY);
				temp[i] = tempX[i] * tempY[i];
			}
			double s = sumatoria(temp, this->tama);
			delete[]tempX;
			delete[]tempY;
			delete[]temp;
			return s;
		}

		double minimosCuadrados::sXYOrden(int orden){
			return this->sXYOrden(orden,orden);
		}

		double minimosCuadrados::sumatoria(double * a, int n){
			double s=0;
			for (int i = 0; i < n ; i++){
				s += a[i];
			}
			return s;
		}



		void minimosCuadrados::printCoeficientes(){
			cout << "Resultado: " << endl;
			cout  << this->coeficientes[0];
			for (int i = 1 ; i < this->tamaCoefs; i++){
				if(fabs(this->coeficientes[i]>CERO)){
					if(this->coeficientes[i]==1){
						cout  << " + " << "x^" << i;
					}
					else{
						cout  << " + "<< this->coeficientes[i] << "x^" << i;
					}
				}
			}
			cout << endl;
		}

		double minimosCuadrados::hConstanteAll(){
			double h = this->equis[1] - this->equis[0];
			double oldh;
			//errNorm = fabs(pi2 - pi1)*100.0 / pi2;
			for (int i = 2; i < this->tama; i++){
				oldh = h;
				h = this->equis[i] - this->equis[i - 1];
				if (h - oldh>CERO){
					return -1;
				}
			}
			return h;
		}

		int minimosCuadrados::ordenOptimo(){ //No usar!
			int opt;
			opt = this->hConstanteAll();
			while(this->hConstanteAll()==-1){
				opt = this->hConstanteAll();
			}
			return opt;
		}

		void minimosCuadrados::print(){
			cout << "X | Y" << endl;
			for (int i = 0; i< this->tama ; i++){
				cout << this->equis[i] <<" | " << this->ye[i] << endl;
			}
			cout << "------" << endl;
			cout << "Sx: " << this->Sx << endl;
			cout << "Sy: " << this->Sy << endl;
		}

		void minimosCuadrados::usar(){
					cout << "Mínimos Cuadrados" << endl;
					int meh = 0;
					do{
						cout << "#############################################################################" << endl;
						cout << "1. Resolver" << endl;
						cout << "2. Leer otro archivo" << endl;
						cout << "3. Imprimir Datos" << endl;
						cout << "0. Salir" << endl;
						cout << endl;
						cin >> meh;
						double orden=-1;
						switch (meh){
						case 1:
							cout << "Resolviendo puntos, ingrese el grado que quiere resolver(-1 para automático):";
							cin >> orden;
							if(orden == -1){
								orden = 5;//this->ordenOptimo();

							}
							this->resolverMinimosCuadrados(orden);
							this->gaussJordan();
							this->printCoeficientes();
							break;
						case 2:
							this->readFile();
							break;
						case 3:
							this->print();
							break;
						case 0:
							cout << "Saliendo" << endl;
							break;
						default:
							cout << "Elige una opcion valida" << endl;
							break;
						}
						cout << endl;
					} while (meh != 0);
				}

		void minimosCuadrados::readFile(){
			string archivo;
			cout << "Introduce el nombre del archivo: " << endl;
			cin >> archivo;
			vector<double> x;
			vector<double> y;
			double v;
			ifstream ifs;
			ifs.open(archivo.c_str());
			while (!ifs.eof()){
				ifs >> v;
				x.push_back(v);
				ifs >> v;
				y.push_back(v);
			}
			ifs.close();
			this->tama = x.size();
			this->equis = new double[this->tama];
			this->ye = new double[this->tama];
			for (int i = 0; i < this->tama; i++){
				equis[i] = x.at(i);
				ye[i] = y.at(i);
			}
			this->Sx = 0;
			this->Sy = 0;
			this->generarSumatorias();
			this->tamaCoefs=0;
			this->coeficientes = new double[this->tamaCoefs];
			this->coeficientesNoResueltos = new Matrix();
		}
