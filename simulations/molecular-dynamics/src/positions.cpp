#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
using namespace std;
#include "constantes.h"


// N nombre de particules
//rho = nombre de particule par unité de volume
double initialisePositions(vector<vector<double>>& R, int& N, double&  rho)
{
    int M = 1;
    while (M * M * M < N) M++;
    M--;
    
    // largeur du cube calcule a partir de N et rho
    double L = pow((M*M*M) / rho, 1.0/3);

    double a = L / M;           // constante de la maille

    cout<<"Initialisation des positions"<<endl;
    cout<<"L = "<<L<<endl;
    cout<<"MMM = "<<M*M*M<<endl;
    cout<<"Old N = "<<N<<endl;
    N = M*M*M; //corection de la valeur de N avec reference
    cout<<"New N = "<<N<<endl;
    cout<<"a = "<<a<<endl;
    rho = N/(L*L*L);

    R = vector<vector<double>>(3, vector<double> (N,0));

    int n = 0;
    for (int x = 0; x < M; x++)
    for (int y = 0; y < M; y++)
    for (int z = 0; z < M; z++)
    {
           R[0][n] = (x + 0.25) * a;
           R[1][n] = (y + 0.25) * a;
           R[2][n] = (z + 0.25) * a;
           n++;
     }
	
    return L;
}


