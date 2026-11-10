#include<vector>
#include<string>
using namespace std;
void dynamiqueNVET(vector<vector<double>>& R, vector<vector<double>>& V, double L, double m, double sigma, double epsilon, double tpas, double tfin, int method, int savenSteps, double tthermal, string conservation, double T, const double rho, int M, int M_rdf, double constcouplage, vector<double>& NV_param, vector<ofstream>& ofstream_vector);
void dynamiqueNV(vector<vector<double>>& R, vector<vector<double>>& V, double L, double m, double sigma, double epsilon, double tpas, double tfin, int method, int savenSteps, double tthermal, string conservation, double T, const double rho, double constcouplage, vector<double>& NV_param, vector<ofstream>& ofstream_vector);



