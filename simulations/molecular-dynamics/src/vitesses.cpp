#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
using namespace std;
#include "constantes.h"


double gauss ()	//distribution gaussienne
{
     static bool dispo = false;
     static double ga;
     double fac, rsq, v1, v2;
     if (!dispo) {
          do {
               v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
               v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
               rsq = v1 * v1 + v2 * v2;
          } while (rsq >= 1.0 || rsq == 0.0);
          fac = sqrt(-2.0 * log(rsq) / rsq);
          ga = v1 * fac;
          dispo = true;
          return v2*fac;
     } else {
          dispo = false;
          return ga;
     }
}

//temperature instantanne du systeme calculée avec V
double calculeTemperatureTnstantanee(vector<vector<double>>& V, int N, const double m)
{
	double s = 0;
	for (int i = 0; i < N; i++)
	for (int k = 0; k < 3; k++)
            s +=  V[k][i] * V[k][i];
    

	double kt=KbEnAu;
	return (s*m) /(3 * (N - 1)*kt); //T en K avec masse et kt en u.a
}

//rescale des vitesses
void rescaleVitesses(vector<vector<double>>& V, int N, double T, const double m)
{
    double Ti = calculeTemperatureTnstantanee(V, N, m);

    double l = sqrt( T / Ti );// N-1 car par de translation de l'ensemble
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            V[i][n] *= l;
}

void initaliseVitesses(vector<vector<double>>& V, int N, double T, const double m)
{
    cout<<"Initialisation des vitesses"<<endl;
	cout<<" "<<endl;
	cout<<" "<<endl;
	V = vector<vector<double>>(3, vector<double> (N,0));

	for (int n = 0; n < N; n++)
	for (int i = 0; i < 3; i++) V[i][n] = gauss();
    // Adjuste les vitesses pour avoir une vitesse du centre de masse = 0;
	double VCM[3] = {0, 0, 0};
	for (int n = 0; n < N; n++)
	for (int i = 0; i < 3; i++)
            VCM[i] += V[i][n];

	for (int i = 0; i < 3; i++) VCM[i] /= N;

	for (int n = 0; n < N; n++)
	for (int i = 0; i < 3; i++)
		V[i][n] -= VCM[i];

    	// Rescale les vitesses à la temperatre désirée
	rescaleVitesses(V, N, T, m);
}


void saveVitesses(vector<vector<double>>& V)
{
	//ecriture des vitesses initiales dans un fichier .txt
	//unité en u.a
	ofstream file("vitesses_init.txt");
	int N = V[0].size();
    	for (int i = 0; i < N; i++)
	{
		file<<V[0][i]<<" "<<V[1][i]<<" "<<V[2][i]<<endl;
	}
	file<<"vx"<<" "<<"vy"<<" "<<"vz"<<endl;	
	file.close();
}


void ApplyBerendsen(vector<vector<double>>& V, int N, double T, double tau, double m)
{
	// ce qui revient a appliquer tout les tau pas
	double Ti =  calculeTemperatureTnstantanee(V,N,m);
	double scale = sqrt(fabs(1.0 + tau*(T/Ti-1.0)));
	for ( int i = 0; i < N; i++)
		for ( int j = 0; j < 3; j++)
			V[j][i] *= scale;
}

