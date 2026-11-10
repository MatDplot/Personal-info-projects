#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
using namespace std;
#include "constantes.h"

double getMinDiff(double v, const double L) //condition pour le calcul de l'acceleration
{
	if(v>L/2.0) return v-L;
	if(v<-L/2.0) return v+L;
	return v;
}

void applyPeriodicite(vector<vector<double>>& R, double L) //conditions de périodicités
{
	for(size_t i=0;i<R[0].size();i++)
		for(int k=0;k<3;k++)
		{
			if(R[k][i]>L) 
				R[k][i] -= int(R[k][i]/L)*L;
			if(R[k][i]<0) 
				R[k][i] += (1-int(R[k][i]/L))*L;
		}
}

double min(const vector <vector<double> > &V) //valeur min d'un vector de vector avec V de taille MxN
{
    double min=V[0][0]; //changer la valeur de min
    for (size_t n = 0; n < V[0].size(); n++)
    {
        for (size_t i = 0; i < V.size(); i++)
        {       
		if (V[i][n]<min) {min=V[i][n];};     
	}
        
    }
    return min ;
}

double minVect(const vector<double> &V) //valeur min d'un vector
{
    double min=V[0]; //changer la valeur de min
    
    for (size_t i = 0; i < V.size(); i++)
        {       
		  if (V[i]<min) {min=V[i];};     
	      }
        
    return min ;
}

double max(const vector <vector<double> > &V) //valeur max d'un vector de vector avec V de taille MxN

{
       double max=V[0][0] ;            
       for (size_t n = 0; n < V[0].size(); n++)
       {
            for (size_t i = 0; i < V.size(); i++)
            {       
		if (V[i][n]>max) {max=V[i][n];};
            }
        
       }
       return max ;
}

double maxVect(const vector<double> &V) //valeur min d'un vector
{
  double max=V[0] ;            
  for (size_t i = 0; i < V.size(); i++)
  {       
		if (V[i]>max) {max=V[i];};
  }
        
  return max ;
}

double DistributionTh_V(const double m, const double T, double v, int S) //distribution theorique des vitesses
{
	double kt=KbEnAu;
  double factor= sqrt(m/(2*M_PI*kt*T));
    
  if(S==2)  //S==2 pour module des vitesses
    {
      factor=factor*factor*factor;
      factor=v*v*4*M_PI*factor;
    }

	return factor*exp((-m*v*v)/(2*kt*T));
}

vector<double> InitHistogramme(const vector <vector<double> >& V, int M, double& a, double& h)
{
  vector<double> H(M,0); //on crée le vecteur H de M élements et on les initialise tous à 0
  a=min(V);
  double b =max(V);
  a *= 2;
  b *= 2;
  h= (b-a)/M;                //on déclare le pas en fonction de maximum et le minimum
  return H;
}

vector<double> InitHistoVsquare(const vector <vector<double> >& V, int M, double& a, double& h)
{
  vector<double> H(M,0); //on crée le vecteur H de M élements et on les initialise tous à 0

  double b1=maxVect(V[0]);
  double b2=maxVect(V[1]);
  double b3=maxVect(V[2]);

  a=0.0;
  double b=sqrt(b1*b1+b2*b2+b3*b3);
  //b valeur du module des vitesses les plus grandes de V
  
  h= (b-a)/M;                //on déclare le pas en fonction de maximum et le minimum
  return H;
}


vector<double> InitRDF(int M, const double L, double& a, double& h)
{
  vector<double> H(M,0); //on crée le vecteur H de M élements et on les initialise tous à 0  
  double b = (L/2)*BohrToAng;
  h= (b-a)/M;                //on déclare le pas en fonction de maximum et le minimum
  return H;
}


void AddHistogramme(const vector< vector<double> >& V, vector<double> &H, double amin, double h)
{
    for(size_t i=0; i<V[0].size(); i++)
    {
		for(size_t k=0;k<V.size();k++)
            	{
                	int p= (V[k][i]-amin)/h;
			// H.size() est de type size_t => entier >=0
			//  int entier qui peut être <0 ou >=0
			// int(H.size()) : conversion de size_t à int
			// on compare ensuite p de type int et int(H.size()) de type int
                	if(p>=0 and p<int(H.size())) H[p]++; 
		}
    }
}

void AddHistoVsquare(const vector< vector<double> >& V, vector<double> &H, double amin, double h)
{
    for(size_t i=0; i<V[0].size(); i++)
    {
		double Vs=sqrt(V[0][i]*V[0][i]+V[1][i]*V[1][i]+V[2][i]*V[2][i]);
    int p= (Vs-amin)/h;
			// H.size() est de type size_t => entier >=0
			//  int entier qui peut être <0 ou >=0
			// int(H.size()) : conversion de size_t à int
			// on compare ensuite p de type int et int(H.size()) de type int
    if(p>=0 and p<int(H.size())) H[p]++; 
    }
}

void AddRDF(const vector< vector<double> >& R, vector<double> &H_rdf, double h_rdf, const double L, double amin_rdf)
{
   double xij,yij,zij;
   double rijs,rij;
   int N=R[0].size();

   for(int i=0; i<N; i++)
	{
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

			rij=sqrt(rijs); //racine de rijs 

      int p= (abs(-rij*BohrToAng+amin_rdf))/h_rdf;
			if(p>=0 and p<int(H_rdf.size())) H_rdf[p]++;
    }
  }
}


void EnregistreHistogramme(const vector<double>& H,double amin, double h, string nomFichier)
{
    ofstream fh(nomFichier.c_str());
    
    for (size_t i=0; i< H.size();i++){
        	fh<<amin+i*h<<"\t"<<H[i]<<endl; 
        }     
       //h est calculé dans initHisto et représente le pas pour tracer l'histogramme
   
    fh.close();
}

void EnregistreHistoTH(const vector<double>& H,double amin, double h, string nomFichier, const double m, const double T, const int S)
{
    ofstream fth(nomFichier.c_str());
    
    vector<double> H_th (H.size(),0);
    double sommeh=0.0;
    double factor;
    if(S==1)
      factor=3;
    if(S==2)
      factor=1;

    for (size_t i=0; i< H.size();i++)
        {
        	double x=amin+i*h;
          H_th[i]=DistributionTh_V(m,T,x,S);
          sommeh+=H_th[i];
        }    

    for (size_t i=0; i< H.size();i++)
        {
        	double x=amin+i*h;
          fth<<x<<"\t"<<factor*H_th[i]/sommeh<<endl; 
        }     
   
    fth.close();
}

void NormaliseHistogramme(vector<double> &H, int nIter, int nParts)
{
    for(size_t i=0; i<H.size(); i++)
		H[i] /= (nIter*nParts); //on divise par le nombre d'itérations et de particules pour normaliser
}

void NormaliseRDF(vector<double> &H, int nIter, int nParts, const double h_rdf, double rho, double amin_rdf)
{
    for(size_t i=0; i<H.size(); i++)
    {
      double r=i*h_rdf+amin_rdf;
		  H[i] /= (nIter*nParts)*(4*M_PI*r*r*h_rdf*(rho/(BohrToAng*BohrToAng*BohrToAng))); //on divise par le nombre d'itérations et de particules pour normaliser
    }
}


double calculDif(const vector< vector<double> >& R, int i, int j, double dr[], double L)
{
        double ll= L/2;
        for(int c=0;c<3;c++)
        {                //nous avons décrémenté et incrémenté de façon que les particules soient dans le cube
                dr[c]=R[c][i]-R[c][j];
                if(dr[c]>ll) dr[c] -= L;
                if(dr[c]<-ll) dr[c] += L;
        }
        return sqrt( dr[0]*dr[0] +dr[1]*dr[1] +dr[2]*dr[2]); // retourne la norme de la distance entre les particules
}

