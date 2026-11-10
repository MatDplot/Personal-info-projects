#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<string>
using namespace std;
#include "constantes.h"
#include "chemin.h"


const string Method_type[4]={"Euler","Verlet","Leapfrog","RK4"};


void LectureCI(double& epsilon, double& sigma, double& m, int& N, double& T, double& rho, 
	double& tpas, double& tfin, int &method, int& savenSteps, double& tthermal, string& conservation, string fname, int& M, int& M_rdf, double& constcouplage, vector<double>& NV_param)
{
	double NV_Tmin, NV_Tmax, NV_Tpas, NV_temps, NV_plateau;
	int repetitions;
	bool film;
	
	ifstream CI(fname);
	CI>>epsilon>>sigma>>m>>N>>rho>>T>>method>>tpas>>tfin>>tthermal>>savenSteps>>M>>M_rdf>>constcouplage>>conservation>>NV_Tmin>>NV_Tmax>>NV_Tpas>>NV_temps>>NV_plateau>>repetitions>>film;

	//on effectue les changement d'unitées
	epsilon /= HartreeToeV; //eV to Hartree
	sigma *=10/BohrToAng; //nm to Angstrom to u.a
	m /= AuTou; //u to u.a
	// rho est en Ang^-3 
	rho *= BohrToAng*BohrToAng*BohrToAng;
	tpas *= 1e-15;// fs => au
	tpas /= AuTosecond;
	tfin *= 1e-12; // ps => au
	tfin /= AuTosecond;
	tthermal *= 1e-12; // ps => au
	tthermal /= AuTosecond;
	constcouplage *= 1e-12; //ps => au
	constcouplage /= AuTosecond;

	NV_temps *= 1e-12; //ps => au
	NV_temps /= AuTosecond;
	NV_plateau *= 1e-12; //ps => au
	NV_plateau /= AuTosecond;

	//attribution des paramètres pour la méthode NV dans un vecteur
	NV_param[0]=NV_Tmin;
	NV_param[1]=NV_Tmax;
	NV_param[2]=NV_Tpas;
	NV_param[3]=NV_temps;
	NV_param[4]=NV_plateau;
	NV_param[5]=repetitions;
	NV_param[6]=film;
}

void AffichageCI(double epsilon, double sigma, double m, int N, int N_original, double T, double rho, 
	double tpas, double tfin, int method, int savenSteps, double tthermal, string conservation, double constcouplage, vector<double>& NV_param)
{
	//affichage des constantes du fichier CI_def.txt
	cout<<" "<<endl;
	cout<<" "<<endl;
	cout<<"Variables en unités atomiques "<<endl;
	cout<<"epsilon="<<epsilon<<endl;
	cout<<"sigma="<<sigma<<endl;
	cout<<"m="<<m<<endl;
	cout<<"N="<<N<<endl;
	cout<<"N original="<<N_original<<endl;
	cout<<"rho="<<rho<<endl;
	cout<<"T(K)="<<T<<endl;
	if(method==0)
		cout<<"Integration par la methode d'Euler"<<endl;
	if(method==1)
		cout<<"Integration par la methode de Verlet-Vitesse"<<endl;
	if(method==2)
		cout<<"Integration par la methode de Leapfrog"<<endl;
	if(method==3)
		cout<<"Integration par la methode RK4"<<endl;
	cout<<"tpas="<<tpas<<endl;
	cout<<"tfin="<<tfin<<endl;
	cout<<"tthermal="<<tthermal<<endl;
	cout<<"Constante de couplage="<<constcouplage<<endl;
	cout<<"nombre d'iterations="<<tfin/tpas<<endl;
	cout<<"save tout les "<<savenSteps<<" pas "<<endl;
	cout<<"Quantités conservées="<<conservation<<endl;
	
	if(conservation=="NV")
		{
			cout<<"T_min(K)="<<NV_param[0]<<endl;
			cout<<"T_max(K)="<<NV_param[1]<<endl;
			cout<<"pas T (K)="<<NV_param[2]<<endl;
			cout<<"dynamique de ="<<NV_param[3]<<endl;
			cout<<"plateau de ="<<NV_param[4]<<endl;
			cout<<"répétitions ="<<NV_param[5]<<endl;
			cout<<"conservation des positions ="<<NV_param[6]<<endl;
		}

	

	double a = pow(rho,-1.0/3);
	if(a<sigma/2)
	{
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cout<<"ATENTION  distance entre les atomes est trop petite par rapport à la distance  d'equilibre " <<endl;
		cout<<"a="<<a<<" sigma = "<<sigma<<endl;
		cout<<"La densité est trop grande"<<endl;
		cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	}
	cout<<" "<<endl;
	cout<<" "<<endl;

}

vector<ofstream> init_result(int type) //ouverture des fichiers pour les méthodes NVE et NVT
{
	//création des fichiers pour les résultats
	system("mkdir -p Resultats/Energies");
	system("mkdir -p Resultats/Positions");
	system("mkdir -p Resultats/Vitesses");
	system("mkdir -p Resultats/Temperatures");
	system("mkdir -p Resultats/Histogrammes/RDF");
	system("mkdir -p Resultats/Histogrammes/Vitesses");

	//titre des fichier en fonction de la méthode


	string fxyz_name=GEO;
	string fv_name=VIT;
	string fe_name=ENE;
	string fT_name=TEM;
	

    fxyz_name.append(Method_type[type]);
	fxyz_name.append(".xyz");

	fv_name.append(Method_type[type]);
	fv_name.append(".txt");

	fe_name.append(Method_type[type]);
	fe_name.append(".txt");

	vector<ofstream> ofstream_vector(4);
	ofstream_vector[0].open(fxyz_name);
	ofstream_vector[1].open(fT_name);
	ofstream_vector[2].open(fv_name);
	ofstream_vector[3].open(fe_name);
	

	return ofstream_vector;
}

vector<ofstream> init_resultNV(int type) //ouverture des fichiers pour la méthode NV
{
	//création des fichiers pour les résultats
	system("mkdir -p Resultats/NV");

	//titre des fichier en fonction de la méthode
	string fxyz_name=GEONV;
	string fT_name=TEMNV;
	

    fxyz_name.append(Method_type[type]);
	fxyz_name.append(".xyz");


	vector<ofstream> ofstream_vector(2);
	ofstream_vector[0].open(fxyz_name);
	ofstream_vector[1].open(fT_name);

	return ofstream_vector;
}
