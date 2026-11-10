#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
using namespace std;
#include "data.h"
#include "positions.h"
#include "vitesses.h"
#include "dynamique.h"

int main()
{
	//constante du fichier
	double epsilon;
	double sigma;
	double m;
	int N;
	double T;
	double rho;
	double tfin;
	double tpas;
	int method;
	int savenSteps;
	double tthermal;
	string conservation;
	int M;
	int M_rdf;
	double cstecouplage;
	
	//vector des fichiers de résulats
	vector<ofstream> ofstream_vector;
	
	//vector des parametres de la conservation NV
	vector<double> NV_param (7,0);


	LectureCI(epsilon, sigma, m, N, T, rho, tpas, tfin, method, savenSteps, tthermal, conservation, "CI_def.txt", M, M_rdf, cstecouplage, NV_param); 
	int N_original=N;	//valeur de N original
	vector<vector<double>> R;	//vecteur de 3 vecteurs pour x,y et z
	double L = initialisePositions(R, N, rho);// N, rho sont probablement modifiés
	AffichageCI( epsilon, sigma,  m, N,  N_original,  T,  rho, tpas,  tfin,  method,  savenSteps, tthermal, conservation, cstecouplage, NV_param);

	//vecteur de 3 vecteurs pour vx,vy et vz
	vector<vector<double>> V;
	initaliseVitesses(V, N, T, m);	//initialisation des vitesses

	if(conservation=="NVE" || conservation=="NVT")
		{
			ofstream_vector=init_result(method);
			dynamiqueNVET(R, V, L,  m, sigma, epsilon, tpas, tfin, method, savenSteps, tthermal, conservation, T, rho, M, M_rdf, cstecouplage, NV_param, ofstream_vector);
		}
	if(conservation=="NV")
		{
			ofstream_vector=init_resultNV(method);
			dynamiqueNV(R, V, L, m, sigma, epsilon, tpas, tfin, method, savenSteps, tthermal, conservation, T, rho, cstecouplage, NV_param, ofstream_vector);
		}
}
