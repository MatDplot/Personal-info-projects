#include<vector>
#include<string>
using namespace std;


double gauss ();
double calculeTemperatureTnstantanee(vector<vector<double>>& V, int N, const double m);
void rescaleVitesses(vector<vector<double>>& V, int N, double T, const double m);
void initaliseVitesses(vector<vector<double>>& V, int N, double T, const double m);
void saveVitesses(vector<vector<double>>& V);
void LectureVitessesI(vector<vector<double>>& A, const int N, string fname);
void ApplyBerendsen(vector<vector<double>>& V, int N, double T, double tau, double m);
