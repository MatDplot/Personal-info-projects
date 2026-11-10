#include<vector>
using namespace std;

double calculEnergiePotentielle(vector<vector<double>>& R, const double L, const double epsilon, const double sigma);
double calculEnergieCinetique(vector<vector<double>>& V, const double m);
double calculPression(const double rho, const int N, const double m, vector<vector<double>>& V, vector<vector<double>>& R, const double sigma, const double epsilon, const double  L);
