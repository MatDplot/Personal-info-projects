#include<vector>
#include<string>
using namespace std;
void CalculAccelerationV2(vector<vector<double>>& A, vector<vector<double>>& R, const double L, const double m, const double epsilon, const double sigma);
void CalculAcceleration(vector<vector<double>>& A, vector<vector<double>>& R, const double L, const double m, const double epsilon, const double sigma);
void Euler(vector<vector<double>>& A, vector<vector<double>>& R, vector<vector<double>>& V, const double tpas, const double L, const double m, const double sigma, const double epsilon);
void Verlet(vector<vector<double>>& A, vector<vector<double>>& R, vector<vector<double>>& V, const double tpas, const double L, const double m, const double sigma, const double epsilon);
void Leapfrog(vector<vector<double>>& A, vector<vector<double>>& R, vector<vector<double>>& V, const double tpas, const double L, const double m, const double sigma, const double epsilon);
void RK4(vector<vector<double>>& A, vector<vector<double>>& R, vector<vector<double>>& V, const double tpas, const double L, const double m, const double sigma, const double epsilon);


