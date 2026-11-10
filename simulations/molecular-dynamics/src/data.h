#include<vector>
using namespace std;

void LectureCI(double& epsilon, double& sigma, double& m, int& N, double& T, double& rho, 
	double& tpas, double& tfin, int &method, int& savenSteps, double& tthermal, string& conservation, string fname, int& M, int& M_rdf, double& constcouplage, vector<double>& NV_param);
void AffichageCI(double epsilon, double sigma, double m, int N, int N_original, double T, double rho, 
	double tpas, double tfin, int method, int savenSteps, double tthermal, string conservation, double constcouplage, vector<double>& NV_param);
vector<ofstream> init_result(int type);
vector<ofstream> init_resultNV(int type);
