#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
using namespace std;
#include "constantes.h"
#include "utiles.h"
#include "energies.h"
#include "vitesses.h"


void CalculAccelerationV2(vector<vector<double>>& A, vector<vector<double>>& R, const double L, const double m, const double epsilon, const double sigma)
{
       for(size_t i=0;i<R[0].size();i++)
        {
                for(int c=0;c<3;c++) A[c][i] =0;    //initialisation de l'acceleration 

                for(size_t j=0;j<R[0].size();j++)
                {
                   if(i!=j)
                   {
                        double  dr[3];
                        double rij = calculDif(R, i, j, dr, L);
                        double sr6 = pow(sigma/rij,6.0);
                        double l = 24*epsilon*(2*sr6*sr6-sr6)/(rij*rij);
                        for(int c=0;c<3;c++)
                                A[c][i] += l*dr[c]/m;      //définition de l'accélération 
                   }
                }
        }
}

void CalculAcceleration(vector<vector<double>>& A, vector<vector<double>>& R, const double L, const double m, const double epsilon, const double sigma)
{
	int N = R[0].size();
	//definitions des distances entre l'atomes i et j
	double xij;
	double yij;
	double zij;


	double rijs;   //distance r au carré
	double rij;		// distance r

	//cste du Potentiel de Lennard-Jones
	double cstA=4*epsilon*6*pow(sigma,6.0);
	double cstB=4*epsilon*12*pow(sigma,12.0);

	
	for(int i=0; i<N; i++)
	{
		double sommeEpX=0.0;
		double sommeEpY=0.0;
		double sommeEpZ=0.0;

		
		for(int j=0; j<N; j++)
		{
			//si j==i on continue (meme particule)
			if(j==i)
				{continue;}
				
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

			//calcul des dérivées des energies potentiels pour avoir la force selon x,y et z
			double rij2 = 1/(rij*rij);
			double rij6 = rij2*rij2*rij2;
			double rr = cstA*rij6*rij2 - cstB*rij6*rij6*rij2;
			sommeEpX += rr*xij;
			sommeEpY += rr*yij;
			sommeEpZ += rr*zij;
		}
		//caclul des accélérations et affectation dans les vecteurs
		A[0][i]=-sommeEpX/(m);
		A[1][i]=-sommeEpY/(m);
		A[2][i]=-sommeEpZ/(m);
	}
}


void Euler(vector<vector<double>>& A, vector<vector<double>>& R, vector<vector<double>>& V, const double tpas, 
 const double L, const double m, const double sigma, const double epsilon)
{
	int N = R[0].size();
	for(int i=0; i<N; i++)
	{
		//calcul des vitesses
		for(int k=0;k<3;k++)
		{
        	V[k][i] += A[k][i]*tpas;
			
		//attribution des positions vérifiées dans les vecteurs x,y et z
		    R[k][i] += V[k][i]*tpas;
        }
	}
	applyPeriodicite(R,L);// remettre les particules dans la boite
	//CalculAcceleration(A,R,L,m,epsilon,sigma);// A t+dt
	CalculAccelerationV2(A,R,L,m,epsilon,sigma);// A t+dt
}



void Verlet(vector<vector<double>>& A, vector<vector<double>>& R, vector<vector<double>>& V, const double tpas,
 const double L, const double m, const double sigma, const double epsilon)
{
	int N = R[0].size();
	for(int i=0; i<N; i++)
		{
			//calcul des positions au temps t+1
			for(int k=0;k<3;k++)
				R[k][i] += V[k][i]*tpas + 0.5*A[k][i]*tpas*tpas;
		}
	applyPeriodicite(R,L);// remettre les particules dans la boite
// On ajoute d'abord le terme avec l'acce à t
	for(int i=0; i<N; i++)
		{
			//calcul des vitesses au temps t+1
			for(int k=0;k<3;k++)
				V[k][i] += 0.5*A[k][i]*tpas;
		}

	//calcul de l'accerelation au temps t+dt
	//CalculAcceleration(A,R,L,m,epsilon,sigma);
	CalculAccelerationV2(A,R,L,m,epsilon,sigma);

// On ajoute d'abord le terme avec l'acce à t+dt
	for(int i=0; i<N; i++)
		{
			//calcul des vitesses au temps t+1
			for(int k=0;k<3;k++)
				V[k][i] += 0.5*A[k][i]*tpas;
		}
}


void Leapfrog(vector<vector<double>>& A, vector<vector<double>>& R, vector<vector<double>>& V, const double tpas,
 const double L, const double m, const double sigma, const double epsilon)
{
	int N = R[0].size();
	for(int i=0; i<N; i++)
	{
		//calcul des vitesses au temps t+1/2
		for(int k=0;k<3;k++)
		{	
            V[k][i] += 0.5*A[k][i]*tpas;
			
		//calcul des positions au temps t+1 
		    R[k][i] += V[k][i]*tpas;
        }
	}
	
	applyPeriodicite(R,L);// remettre les particules dans la boite
	//CalculAcceleration(A,R,L,m,epsilon,sigma);// A t+dt
	CalculAccelerationV2(A,R,L,m,epsilon,sigma);// A t+dt
	
	for(int i=0; i<N; i++)
	{
		//calcul des vitesses au temps t+1
		for(int k=0;k<3;k++)
			V[k][i] += 0.5*A[k][i]*tpas;
	}
}


void RK4(vector<vector<double>>& A, vector<vector<double>>& R, vector<vector<double>>& V, const double tpas,
 const double L, const double m, const double sigma, const double epsilon)
{
	int N = R[0].size();

    vector<vector<double>> R1(3, vector<double> (N,0));
    vector<vector<double>> R2(3, vector<double> (N,0));
    vector<vector<double>> R3(3, vector<double> (N,0));
    vector<vector<double>> R4(3, vector<double> (N,0));

    vector<vector<double>> V1(3, vector<double> (N,0));
    vector<vector<double>> V2(3, vector<double> (N,0));
    vector<vector<double>> V3(3, vector<double> (N,0));
    vector<vector<double>> V4(3, vector<double> (N,0));

    vector<vector<double>> R_temp(3, vector<double> (N,0));
    vector<vector<double>> V_temp(3, vector<double> (N,0));

	for(int i=0; i<N; i++)
	{
		for(int k=0;k<3;k++)
		{	
            V1[k][i] = A[k][i]*tpas;
			R1[k][i] = V[k][i]*tpas;
        }
	}
	
    for(int i=0; i<N; i++)
	{
		for(int k=0;k<3;k++)
		{	
            V_temp[k][i] = V[k][i] + V1[k][i]/2.0;
		    R_temp[k][i] = R[k][i] + R1[k][i]/2.0;
        }
	}
    
    applyPeriodicite(R_temp,L);
    //CalculAcceleration(A,R,L,m,epsilon,sigma);// A t+dt
	CalculAccelerationV2(A,R_temp,L,m,epsilon,sigma);// A t+dt
	
    for(int i=0; i<N; i++)
	{
		for(int k=0;k<3;k++)
		{	
            V2[k][i] = A[k][i]*tpas;
			R2[k][i] = V_temp[k][i]*tpas;
        }
	}

    for(int i=0; i<N; i++)
	{
		for(int k=0;k<3;k++)
		{	
            V_temp[k][i] = V[k][i] + V2[k][i]/2.0;
		    R_temp[k][i] = R[k][i] + R2[k][i]/2.0;
        }
	}

    applyPeriodicite(R_temp,L);
    //CalculAcceleration(A,R,L,m,epsilon,sigma);// A t+dt
	CalculAccelerationV2(A,R_temp,L,m,epsilon,sigma);// A t+dt
	
    for(int i=0; i<N; i++)
	{
		for(int k=0;k<3;k++)
		{	
            V3[k][i] = A[k][i]*tpas;
			R3[k][i] = V_temp[k][i]*tpas;
        }
	}

    for(int i=0; i<N; i++)
	{
		for(int k=0;k<3;k++)
		{	
            V_temp[k][i] = V[k][i] + V3[k][i];
		    R_temp[k][i] = R[k][i] + R3[k][i];
        }
	}

    applyPeriodicite(R_temp,L);
    //CalculAcceleration(A,R,L,m,epsilon,sigma);// A t+dt
	CalculAccelerationV2(A,R_temp,L,m,epsilon,sigma);// A t+dt
	
    for(int i=0; i<N; i++)
	{
		for(int k=0;k<3;k++)
		{	
            V4[k][i] = A[k][i]*tpas;
			R4[k][i] = V_temp[k][i]*tpas;
        }
	}

    for(int i=0; i<N; i++)
	{
		for(int k=0;k<3;k++)
		{	
            V[k][i] += (V1[k][i] + 2*V2[k][i] + 2*V3[k][i] + V4[k][i])/6.0;
			R[k][i] += (R1[k][i] + 2*R2[k][i] + 2*R3[k][i] + R4[k][i])/6.0;
        }
	}

	applyPeriodicite(R,L);
    //CalculAcceleration(A,R,L,m,epsilon,sigma);// A t+dt
	CalculAccelerationV2(A,R,L,m,epsilon,sigma);// A t+dt
}

