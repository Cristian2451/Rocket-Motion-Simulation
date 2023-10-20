#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
using namespace std;

// Define Class containing values of environmental parameters
class EParameters {
public:
    double rho;      // Density of air (kg/m^2)
    double g;        // Gravitational acceleration (m/s^2)
    double Cd;       // Drag coefficient of air
    double v0;       // Initial velocity (m/s)
    double h0;       // Initial height (m)
};


// Define Class containing values of rocket parameters
// Parameters stored using vectors, as if rocket is multistage, each stage will have its own parameters
class RParameters {
public:
    vector<double> m0;       // Initial fuel mass (kg)
    vector<double> mr;       // Rocket mass (kg)
    vector<double> A;        // Area (m^2)
    vector<double> mdot_f;   // Fuel mass flow rate (kg/s)
    vector<double> ue;       // Exhaust velocity (m/s)
    vector<double> dt;       // Time step (s)
    
    // Function that sets size to all the vector instance variables of class (where n = n° of stages)
    void SetSize(int n){
        m0 = vector<double>(n);     
        mr = vector<double>(n);       
        A = vector<double>(n);        
        mdot_f = vector<double>(n);   
        ue = vector<double>(n);       
        dt = vector<double>(n);
    }
};


// Define Class containing instance variables (t,h,v,m)
class Data {
public:
    double t;      // Time (s)
    double h;      // Height (m)
    double v;      // Velocity (m/s)
    double m;      // Mass (kg)
    
    // Data class simple constructor
    Data(double tc, double hc, double vc, double mc){
        t = tc;
        h = hc;
        v = vc;
        m = mc;
    }
    //Data class deep-copy constructor
    Data(const Data& d){
        t = d.t;
        h = d.h;
        v = d.v;
        m = d.m;
    }
    // Stream insertion operator overload for Data class
    friend ostream& operator<<(ostream& out, const Data& d){
        out << setw(9) << d.t << "  " << setw(9) << d.h << "  " << setw(9) << d.v << "  " << setw(7) << d.m;
        return out;
    }
};


// Function that finds n° of stages by counting lines of txt file
int CountStages() {
    int n = 0;
    ifstream Parameters("parameters.txt");
    string line;
    while (getline(Parameters, line)) {
        n++;
    }
    return n-1; // remove environmental parameters line from count
}


// Define function that reads parameters.txt file and assigns values to EParameters and RParameters classes
void readVariables(EParameters &EV, RParameters &RV, int n) {
    vector<double>data;
    ifstream Parameters("parameters.txt");
    
    if (Parameters.good()) {   // Check file opens ok
        Parameters >> EV.rho >> EV.g >> EV.Cd >> EV.v0 >> EV.h0; // Asigns values of first line to the environmental parameters
        string line;
        // Read the next lines one at a time into the variable 'line'
        while(getline(Parameters, line)) {
            stringstream lineStream(line);
            double value;
            vector<double>lineData;
            // Read a value at a time from 'line' and store in a vector
            while(lineStream >> value) {
                lineData.push_back(value);
            }
            data.insert(data.end(),lineData.begin(), lineData.end()); //concatenate each vector to a sigle vector 'data'
        }
        // Assign values in vector 'data' to the equivalent rocket parameter (and stage)
        for(int i = 0; i < n; i++){
            RV.m0[i] = data[i*6];       
            RV.mr[i] = data[1+i*6];    
            RV.A[i] = data[2+i*6];        
            RV.mdot_f[i] = data[3+i*6]; 
            RV.ue[i] = data[4+i*6];      
            RV.dt[i] = data[5+i*6]; 
        }
    }
    else {
        cout << "Failed to open file" << endl << endl;
    }
}

// Function that ouputs vector containing runge kutta constant for h, v and m respectively
vector<double> f(double h, double v, double m, EParameters &EV, RParameters &RV, int stage){
    vector<double> k(3);
    k[0] = RV.dt[stage-1]*v; 
    k[1] = RV.dt[stage-1]*(-EV.g - EV.rho*abs(v)*v*EV.Cd*RV.A[stage-1]/(2*m) + (RV.mdot_f[stage-1]*RV.ue[stage-1]/(m))*((v >= 0) - (v < 0)));   // (v >= 0) - (v < 0) used to find sign of velocity
    k[2] = RV.dt[stage-1]*(-RV.mdot_f[stage-1]);
    return k;
}


// Function performs 4th order Runge Kutta for a single rocket stage and returns final Data (t,h,v,m)
Data RungeKutta(Data var, EParameters &EV, RParameters &RV, double fuel_mass, int stage, int N){
    
    vector<double> k1(3), k2(3), k3(3), k4(3), k(3); // Runge Kutta constants
    
    // When code initialized (t=0), delete output.txt file if it already exists
    if (var.t == 0) {
        ofstream ClearFile("output.txt", ios::trunc);
    }
    ofstream vOut("output.txt", ios::out | ios::app); // Open file for output and append data to end of file if it already exists
    if (vOut.good()) {   // Check file opens ok
        vOut.precision(7);
        vOut << Data(var.t, var.h, var.v, var.m) << '\n'; // Output Data to output.txt file 
        
        int n = 3;   // Number of variables undergoing Runge Kutta (h,v,m)
        
        while(var.h >= 0){    // While height > 0, iterate Runge Kutta
            k1 = f(var.h, var.v, var.m, EV, RV, stage);
            k2 = f(var.h + k1[0] / 2, var.v + k1[1] / 2, var.m + k1[2] / 2, EV, RV, stage);
            k3 = f(var.h + k2[0] / 2, var.v + k2[1] / 2, var.m + k2[2] / 2, EV, RV, stage);
            k4 = f(var.h + k3[0], var.v + k3[1], var.m + k3[2], EV, RV, stage);
            for (int j = 0; j < n; j++){
                k[j] = (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
            }
            // If new mass greater than or equal to rocket mass, update values
            if (var.m + k[2] >= RV.m0[stage-1] - fuel_mass) {
                var = Data(var.t + RV.dt[stage-1], var.h + k[0], var.v + k[1], var.m + k[2]);
                vOut << Data(var.t, var.h, var.v, var.m) << '\n';
            }
            // If new mass less than rocket mass and not in final stage, don't update values and stop iterating Runge Kutta
            else if (var.m + k[2] < (RV.m0[stage-1] - fuel_mass) && stage < N) {
                break;
            }
            // If new mass less than rocket mass and in final stage, don't update values, set mass equal to rocket mass
            else if (var.m + k[2] < (RV.m0[stage-1] - fuel_mass) && stage == N){
                var = Data(var.t, var.h, var.v, RV.m0[stage-1] - fuel_mass);
                vOut << Data(var.t, var.h, var.v, var.m) << '\n';
                RV.mdot_f[stage-1] = 0; // set fuel mass flow rate to zero
            }
            // (for last stage) if mass equal to rocket mass update values (except mass)
            else if (var.m == RV.mr[N-1]){
                var = Data(var.t + RV.dt[stage-1], var.h + k[0], var.v + k[1], var.m);
                vOut << Data(var.t, var.h, var.v, var.m) << '\n';
                RV.mdot_f[stage-1] = 0; // set again just in case it wasn't called before
            }
        }
    }
    else {
        cout << "Failed to open file" << endl;
    }
    return var; // returns final values (t,h,v,m)
}

// Function call RungeKutta function for each stage, and modifies required parameters on each stage sepparation
void RocketSimulator(EParameters &EV, RParameters &RV, int N){
    int t0 = 0;  // Initial time (s)
    int stage = 1; // Current stage
    Data var (t0, EV.h0, EV.v0, RV.m0[0]); // Specify initial values of t,h,v,m in the Data class called 'var' (if multistage, initial mass not yet accurate)
    // Iterate through 'N' stages
    for(int i = 0; i < N; i++){
        double fuel_mass = RV.m0[i] - RV.mr[i]; // Calculate fuel mass for stage number (i+1)
        for(int j = i+1; j < N; j++){
            RV.m0[i] += RV.m0[j];  // Calculate total mass of rocket with the stages left
        }
        var.m = RV.m0[i]; // Update initial mass, or mass after stage sepparation (updates only for multistage rocket)
        var = RungeKutta(var, EV, RV, fuel_mass, stage, N); // Performs runge kutta for current stage
        stage++; 
    }
    cout << "Results obtained, check output.txt file" << endl;
}


int main()
{
    int N = CountStages(); // Total stages
    EParameters EV; // Environmental parameters class
    RParameters RV; // Rocket parameters class
    RV.SetSize(N); // Sets size to vector objects in RV class
    readVariables(EV, RV, N); // Assigns values of "parameters.txt" to its equivalent class instance variable 
    RocketSimulator(EV, RV, N);
}