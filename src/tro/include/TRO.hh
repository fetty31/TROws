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

    MatrixXd Pleft, Pright; // Gates
    int N; // Dimension
    MatrixXd pointsSol, coefsSplinesTrajX, coefsSplinesTrajY;// Points of the minimum curvature trajectory solution and its coefficients of the splines
    VectorXd splinesLengths; // Length of every spline
    VectorXd radiCurv; // Curvature of the trajectory
    VectorXd freeL, freeR; // Space from the pointTraj to the limit of the track, left and right
    MatrixXd pointsTraj; 
    VectorXd velocity; // Velocity profile
    kdt::KDTree<Point> trajTree; // KDTree for the planner 

};

// TRO class
class TRO{

    private:

        // Internal variables of TRO
        MatrixXd filecontent; // Eigen Matrix to store file content
        MatrixXd Xopt; // Matrix to store MPC stages
        MatrixXd middlePoints; // Matrix to store middle path trajectory

        VectorXd midX; // Middle trajectory x coord
        VectorXd midY; // Middle trajectory y coord
        VectorXd Pheading; // Heading of points from middle trajectory

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