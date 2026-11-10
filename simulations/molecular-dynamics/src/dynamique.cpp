#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
using namespace std;
#include "equa_diff.h"
#include "constantes.h"
#include "utiles.h"
#include "energies.h"
#include "vitesses.h"

#define fxyz ofstream_vector[0]
#define fT ofstream_vector[1]
#define fv ofstream_vector[2]
#define fe ofstream_vector[3]



void (*MethodeIntegration[4]) (vector<vector<double>>& , vector<vector<double>>& , vector<vector<double>>& , const double , const double , const double , const double , const double )={Euler, Verlet, Leapfrog, RK4};
//pointeur de fonction 

void addPositionsToFile(ofstream& f, vector<vector<double>>& R, double L)
{
	f<<R[0].size()+8<<endl;
	f<<"cfc"<<endl;
	for(size_t i=0;i<R[0].size();i++)
		f<<"Ar " <<R[0][i]*BohrToAng <<" "<<R[1][i]*BohrToAng <<" "<<R[2][i]*BohrToAng<<endl;
	// Les bords
	double l = L*BohrToAng;
	f<<"Xe "<<0.0<<" "<<0.0<<" "<<0.0<<endl;
	f<<"Xe "<<l<<" "<<0.0<<" "<<0.0<<endl;
	f<<"Xe "<<l<<" "<<l<<" "<<0.0<<endl;
	f<<"Xe "<<0<<" "<<l<<" "<<0.0<<endl;
	f<<"Xe "<<0.0<<" "<<0.0<<" "<<l<<endl;
	f<<"Xe "<<l<<" "<<0.0<<" "<<l<<endl;
	f<<"Xe "<<l<<" "<<l<<" "<<l<<endl;
	f<<"Xe "<<0<<" "<<l<<" "<<l<<endl;
}

void addVitesses(ofstream& f, vector<vector<double>>& V)
{
	for(size_t i=0;i<V[0].size();i++)
		f<<V[0][i] <<" " <<V[1][i] <<" "<<V[2][i]<<endl;// en unite atomique
}


//thermalisation du systeme
void Thermalisation(double tinit, double tthermal, double tpas, int method, vector<vector<double>>& A, vector<vector<double>>& R, vector<vector<double>>& V, double L, double m,
double sigma, double epsilon, string conservation, ofstream& f, vector<double>& NV_param, double constcouplage, int N, double T)
{
	for(double t=tinit;t<tthermal;t+=tpas)
	{
		cout<<t/tthermal*100<<"\r %";
		
		//résolution des equations du mouvement avec la méthode d'Euler
		//if(method==1)
		//	Euler(A,R,V,tpas,L, m, sigma, epsilon);

		//résolution des equations du mouvement avec la méthode de Verlet	
		//if(method==2)
		//	Verlet(A,R,V,tpas,L,m,sigma,epsilon);

		//résolution des equations du mouvement avec la méthode de Leapfrog	
		//if(method==3)
		//	Leapfrog(A,R,V,tpas,L,m,sigma,epsilon);
		
		//résolution des equations du mouvement avec la méthode de RK4	
		//if(method==4)
		//	RK4(A,R,V,tpas,L,m,sigma,epsilon);

		MethodeIntegration[method](A,R,V,tpas,L,m,sigma,epsilon); //pointeur de fonction permet d'éviter les lignes ci-dessus
				
		if(conservation=="NVE" || conservation=="NVT")
		{
			rescaleVitesses(V,N,T,m);
			f<<(-tthermal+t)*AuTosecond/1e-12<<" "<<calculeTemperatureTnstantanee(V,N,m)<<endl;
		}
		if(conservation=="NV")
		{
			rescaleVitesses(V,N,NV_param[0],m);	
		}
	}
}


// R, V et les autres données sont à definir dans main
void dynamiqueNVET(vector<vector<double>>& R, vector<vector<double>>& V, double L, double m, double sigma, double epsilon, double tpas, 
	double tfin, int method, int savenSteps, double tthermal, string conservation, double T, const double rho, int M, int M_rdf, double constcouplage, vector<double>& NV_param, vector<ofstream>& ofstream_vector)
{
	cout<<"Je suis dans dynamique "<<endl;
	int N = R[0].size(); //nombre de particules
	vector<vector<double>> A(3, vector<double> (N,0));


	//paramètres de temps
	double tinit=0.0;
	int compteurTemps=0; //compteur pour les sauvegardes des positions et vitesses
	

	fT<<"# t          T	"<<endl;

	//initialisation histogramme des vitesses vx, vy et vz
	double amin, h;
	vector<double> H=InitHistogramme(V, M, amin, h);

	//initialisation histogramme du module des vitesses
	double aminvs, hvs;
	vector<double> Hvs=InitHistoVsquare(V, M, aminvs, hvs);

	//initialisation de la RDF
	double amin_rdf, h_rdf;
	amin_rdf=sigma*BohrToAng*0;
	vector<double> H_rdf=InitRDF(M_rdf,L,amin_rdf,h_rdf);


	//CalculAcceleration(A,R,L,m,epsilon,sigma);
	CalculAccelerationV2(A,R,L,m,epsilon,sigma);

	cout<<" "<<endl;
	cout<<"Thermalisation"<<endl;
	Thermalisation(tinit,tthermal,tpas,method,A,R,V,L,m,sigma,epsilon, conservation, fT, NV_param, constcouplage, N, T);



	
	cout<<" "<<endl;
	cout<<"Dynamique moléculaire"<<endl;

	fe<<"# t          Ep               Ec               Etot       (Etot-E0)/E0*100"<<endl;		//en tete du fichier des energies
	//CalculAcceleration(A,R,L,m,epsilon,sigma);
	CalculAccelerationV2(A,R,L,m,epsilon,sigma);
	double E0=calculEnergieCinetique(V, m)+calculEnergiePotentielle(R, L, epsilon,sigma); //Energie mecan à t = 0

	//initialisation des compteurs
	double sommePression=0.0;
	double sommeTemperature=0.0;
	int nombreIteration=0;


	for(double t=tinit;t<tfin;t+=tpas)
	{
		cout<<t/tfin*100<<"\r %"; //affichage pourcentage
					
		if(compteurTemps%savenSteps==0) //sauvegarde tout les savenSteps
		{
			addPositionsToFile(fxyz, R, L);
			addVitesses(fv, V);

			//affichage de Ep, Ec et Em
			double sommeEc=calculEnergieCinetique(V, m); //Energie cinétique du systeme 
			double sommeEp=calculEnergiePotentielle(R, L, epsilon,sigma); //Energie potentielle du systeme 
			double E = sommeEp+sommeEc;
			fe<<t*AuTosecond/1e-12<<" "<<sommeEp*HartreeToeV<<"	 "<<sommeEc*HartreeToeV<<" "<<E*HartreeToeV<<" "<<(E-E0)/E0*100<<endl;

		}
			
		compteurTemps++;
			
		fT<<t*AuTosecond/1e-12<<" "<<calculeTemperatureTnstantanee(V,N,m)<<endl;
			
		sommePression += calculPression(rho, N, m, V, R, sigma, epsilon, L);
		sommeTemperature += calculeTemperatureTnstantanee(V,N,m);
		nombreIteration ++;

		//résolution des equations du mouvement choisie
		MethodeIntegration[method](A,R,V,tpas,L,m,sigma,epsilon);


		if(conservation=="NVT") //themostat si NVT
		{
			ApplyBerendsen(V,N,T,tpas/constcouplage, m);
		}

		//mise à jour des histogrammes
		AddHistogramme(V, H, amin, h);
		AddRDF(R, H_rdf, h_rdf, L, amin_rdf);
		AddHistoVsquare(V, Hvs, aminvs, hvs);
	}
	

	//afichages des valeurs

	cout<<" "<<endl;
	cout<<" "<<endl;

	sommePression=sommePression/nombreIteration;
	cout<<"Pression = "<<sommePression<<" (u.a)"<<endl;
	cout<<"Pression = "<<sommePression*2.9421015697e13<<" (Pa)"<<endl;

	cout<<" "<<endl;
	

	sommeTemperature=sommeTemperature/nombreIteration;
	cout<<"Température = "<<sommeTemperature<<" (K)"<<endl;


	fxyz.close();
	fv.close();
	fe.close();
	fT.close();

	cout<<" "<<endl;

	int n = int(tfin/tpas);
	NormaliseHistogramme(H, n, R[0].size());
	EnregistreHistogramme(H,amin, h,"Resultats/Histogrammes/Vitesses/histo_vitesses.txt");// histogramme des vitesses
	EnregistreHistoTH(H,amin, h,"Resultats/Histogrammes/Vitesses/histo_vitesses_th.txt", m, T, 1);
	

	NormaliseHistogramme(Hvs, n, R[0].size());
	EnregistreHistogramme(Hvs,aminvs, hvs,"Resultats/Histogrammes/Vitesses/histo_modules_vitesses.txt");// histogramme des modules des vitesses
	EnregistreHistoTH(Hvs,aminvs, hvs,"Resultats/Histogrammes/Vitesses/histo_modules_vitesses_th.txt", m, T, 2);
	
	

	NormaliseRDF(H_rdf, n, R[0].size(), h_rdf, rho, amin_rdf);
	EnregistreHistogramme(H_rdf,amin_rdf,h_rdf,"Resultats/Histogrammes/RDF/histo_RDF.txt");
	
	cout<<"Pour visualiser les résulats executer les différents fichier .bash dans les dossiers résultats"<<endl;
}





void dynamiqueNV(vector<vector<double>>& R, vector<vector<double>>& V, double L, double m, double sigma, double epsilon, double tpas, 
	double tfin, int method, int savenSteps, double tthermal, string conservation, double T, const double rho, double constcouplage, vector<double>& NV_param, vector<ofstream>& ofstream_vector)
{
	cout<<"Je suis dans dynamique "<<endl;
	int N = R[0].size();
	vector<vector<double>> A(3, vector<double> (N,0));

	//paramètres de temps
	double tinit=0.0;
	int compteurTemps=0; //compteur pour les sauvegardes des positions et vitesses

	
	fT<<"# t          T	"<<endl;



	//CalculAcceleration(A,R,L,m,epsilon,sigma);
	CalculAccelerationV2(A,R,L,m,epsilon,sigma);

	cout<<" "<<endl;
	cout<<"Thermalisation"<<endl;
	Thermalisation(tinit,tthermal,tpas,method,A,R,V,L,m,sigma,epsilon, conservation, fT, NV_param, constcouplage, N, T);


	cout<<" "<<endl;
	cout<<"Dynamique moléculaire"<<endl;


	//création d'un vecteur contenant toutes les températures à parcourir selon les paramètres rentrés dans le fichier
	vector<double> NV_Tvect(1,NV_param[0]);
	for(double NV_T=NV_param[0]+NV_param[2]; NV_T<NV_param[1]; NV_T+=NV_param[2])
	{
		NV_Tvect.push_back(NV_T);
	}
	for(double NV_T=NV_param[1]; NV_T>=NV_param[0]; NV_T-=NV_param[2])
	{
		NV_Tvect.push_back(NV_T);
	}
			
	double Tempfin=0.0;	//temp à la fin d'une répétition

	for(int j=0; j<NV_param[5];j++)
	{
		cout<<"répétition "<<j+1<<endl; //nombre de répétitions effectuées
		
		int nombreIteration=0;
		int NV_compteurTemps;
		bool Tmax = false; //si plateau = true
		double TempTmax=0.0; //temps du plateau
		
		//on parcours toutes les températures possibles
		for(size_t i=0;i<NV_Tvect.size();i++)
		{
			cout<<i/double(NV_Tvect.size())*100<<"\r %";
				
				
			if(Tmax==true) //plateau attein
			{
				TempTmax=NV_compteurTemps*tpas; //temp du plateau
				Tmax=false;
			}
				

			//si on arrive à la valeur max on réalise une dynamique de 5ps au lieu de celle rentrée par l'utilisateur dans le fichier CI_def
			double NV_tfin=NV_param[3];
			if(NV_Tvect[i]==NV_param[1])
			{
				NV_tfin= NV_param[4];
				Tmax=true;
			}					
				

			NV_compteurTemps=0;
			for(double t=tinit;t<NV_tfin;t+=tpas)
			{
				NV_compteurTemps++;
				if(compteurTemps%savenSteps==0) 
				{
					if(NV_param[6]==true)
						addPositionsToFile(fxyz, R, L);
						
						
					if(NV_Tvect[i]==NV_param[1])
						fT<<(i*NV_param[3]+t+Tempfin)*AuTosecond/1e-12<<" "<<calculeTemperatureTnstantanee(V,N,m)<<endl;
					else
						fT<<(tpas*NV_compteurTemps+i*NV_param[3]+TempTmax+Tempfin)*AuTosecond/1e-12<<" "<<calculeTemperatureTnstantanee(V,N,m)<<endl;
				}
				
				compteurTemps++;
				nombreIteration ++;

				//résolution des equations du mouvement choisie
				MethodeIntegration[method](A,R,V,tpas,L,m,sigma,epsilon);

				ApplyBerendsen(V,N,NV_Tvect[i],tpas/constcouplage, m);
			}
		}
		Tempfin+=tpas*NV_compteurTemps+NV_Tvect.size()*NV_param[3]+TempTmax;
	}

	addPositionsToFile(fxyz, R, L); //sauvegarde de la position finale

	cout<<" "<<endl;
	cout<<" "<<endl;

	fxyz.close();
	fT.close();

	cout<<"Pour visualiser executer le fichier .bash dans le dossiers résultats puis NV"<<endl;
	cout<<"Executer Jmol pour visualiser la configuration"<<endl;
}

