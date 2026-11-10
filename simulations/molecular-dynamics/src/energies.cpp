#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
using namespace std;
#include "utiles.h"

double calculEnergiePotentielle(vector<vector<double>>& R, const double L, const double epsilon, const double sigma)
{
	//definitions des distances entre l'atomes i et j
	double xij;
	double yij;
	double zij;


	double rijs;   //distance r au carré
	double rij;		// distance r

	//cste du Potentiel de Lennard-Jones
	double cstA=4*epsilon*pow(sigma,6.0);
	double cstB=4*epsilon*pow(sigma,12.0);
	double Ep=0;
	int N = R[0].size();

	
	for(int i=0; i<N; i++)
		{
			for(int j=i+1; j<N; j++) // j > i
			{
				//calcul distance entre atome i et j
				xij=R[0][i]-R[0][j];
				yij=R[1][i]-R[1][j];
				zij=R[2][i]-R[2][j];

				//verification pour le caclul de Ep
				xij=getMinDiff(xij,L);
				yij=getMinDiff(yij,L);
				zij=getMinDiff(zij,L);

				//calcul de la distance r au carré
				rijs=xij*xij+yij*yij+zij*zij;

				//si rijs vaut 0 on lui attribue une toute petite valeur pour ne pas diviser par 0
				if(rijs==0)
					{rijs=1e-15;}

				rij=sqrt(rijs); //racine de rijs

				//calcul de l'énergie potentielle
				double rij6 = 1/pow(rij,6.0);
				Ep += cstB*rij6*rij6-cstA*rij6;

			}
		}
	return Ep;
}

double calculEnergieCinetique(vector<vector<double>>& V, const double m)
{
	double Ec = 0;
	int N = V[0].size();
	
	for(int i=0; i<N; i++)
		Ec += V[0][i]*V[0][i] + V[1][i]*V[1][i] + V[2][i]*V[2][i];

	Ec *= m/2;

	return Ec;
}



double calculPression(const double rho, const int N, const double m, vector<vector<double>>& V, vector<vector<double>>& R, const double sigma, const double epsilon, const double  L)
{
	double Ec=calculEnergieCinetique(V,m);

	double W;

	for(size_t i=0;i<R[0].size();i++)
        {
                for(size_t j=0;j<R[0].size();j++)
                {
                   if(i!=j)
                   {
                        double  dr[3];
                        double rij = calculDif(R, i, j, dr, L);
                        double sr6 = pow(sigma/rij,6.0);
                        double l = 24*epsilon*(2*sr6*sr6-sr6)/(rij*rij);
                        for(int c=0;c<3;c++)
                                W += l*dr[c]*dr[c];     
                   }
                }
        }

	return (rho/(3*N))*(2*Ec+(W/2));
}
