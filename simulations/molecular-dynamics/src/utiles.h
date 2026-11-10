#include<vector>
using namespace std;

double getMinDiff(double v, const double L);
void applyPeriodicite(vector<vector<double>>& R, double L);
vector<double> InitHistogramme(const vector <vector<double> >& V, int M, double& a, double& h);
void AddHistogramme(const vector< vector<double> >& V, vector<double> &H, double amin, double h);
void EnregistreHistogramme(const vector<double>& H,double amin, double h, string nomFichier);
void NormaliseHistogramme(vector<double> &H, int nIter, int nParts);
double calculDif(const vector< vector<double> >& R, int i, int j, double dr[], double L);
vector<double> InitRDF(int M, const double L, double& a, double& h);
void AddRDF(const vector< vector<double> >& R, vector<double> &H_rdf, double h_rdf, const double L, double amin_rdf);
void NormaliseRDF(vector<double> &H, int nIter, int nParts, const double h_rdf, double rho, double amin_rdf);
void AddHistoVsquare(const vector< vector<double> >& V, vector<double> &H, double amin, double h);
vector<double> InitHistoVsquare(const vector <vector<double> >& V, int M, double& a, double& h);
void EnregistreHistoTH(const vector<double>& H,double amin, double h, string nomFichier, const double m, const double T, const int S);


