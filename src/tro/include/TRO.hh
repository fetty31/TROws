#ifndef TRO_HH
#define TRO_HH

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <eigen3/Eigen/Dense>
#include "kdtree.h"


using namespace std;
using namespace Eigen;

// This class is for the KDTree
class Point : public std::array<float,2> {
    public :
        static const int DIM = 2;
};

// Optimized trajectory struct
struct trajectory {

    int N; // Dimension

    MatrixXd pointsTraj, coefsSplinesTrajX, coefsSplinesTrajY;// Optimized trajectory points and its coefficients of the splines
    MatrixXd middlePoints; // Central trajectory points
    MatrixXd Xopt; // Optimized stages (9xN)

    VectorXd Pheading; // Heading of central trajectory
    VectorXd splinesLengths; // Length of every spline
    VectorXd radiCurv; // Curvature of the trajectory
    VectorXd freeL, freeR; // Space from the pointsTraj to the limit of the track, left and right
    VectorXd Vx; // Vx profile
    VectorXd Vy; // Vy profile
    VectorXd w; // Yaw rate profile
    kdt::KDTree<Point> trajTree; // KDTree for the planner 

};

// TRO class
class TRO{

    private:

        // Internal variables of TRO

        const std::string x_opt = "/home/fetty/Escritorio/LaptimeSimulator/TROpy/x_opt.csv"; // Output (optimized stages) of TRO.py
        const std::string problem = "/home/fetty/Escritorio/LaptimeSimulator/TROpy/problem.csv"; // Problem characteristics of TRO.py
        const std::string Pmiddle = "/home/fetty/Escritorio/LaptimeSimulator/TROpy/filtered_points.csv"; // Middle trajectory of TRO.py

        int n_states = 7; // Number of states variables
        int n_controls = 2; // Number of controls variables
        int N;
        int Npar;

        double T;

        void get_data();
        void heading();
        void get_trajectory();
        void radi_curv();
        void create_KDTree();
        MatrixXd coefs_splines(VectorXd x);
        VectorXd polyval(Vector4d coeffs, VectorXd t);
        double polyval2(Vector4d coeffs, double t);
        double integral_length(Vector4d coefsX, Vector4d coefsY);
        MatrixXd read_csv(const std::string filename, bool firstout = true); 

    public:

        TRO(); // Constructor
        void init(); // Initialization function

        // Optimized trajectory data
        trajectory traj = trajectory();
};


#endif