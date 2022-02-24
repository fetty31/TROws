#include "../include/TRO.hh"

using namespace std;
using namespace Eigen;

// Constructor
TRO::TRO(){}

// Initialization of TRO
void TRO::init(){

    // auto start_time = std::chrono::system_clock::now();

    // Run principal functions
    get_data();
    heading();
    get_trajectory();
    // radi_curv();
    // create_KDTree();

    // auto end_time = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed = end_time - start_time;
    // cout << elapsed << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//------------------------Principal functions--------------------------------------------------

// get_data: read TRO.py output and transform the data to desired format
void TRO::get_data(){

    MatrixXd filecontent = read_csv(x_opt); 

    MatrixXd middlePoints = read_csv(Pmiddle);
    traj.middlePoints = middlePoints;

    MatrixXd Problem = read_csv(problem);

    N = Problem(0,0);
    Npar = Problem(2,0);
    T = Problem(1,0);

    // Index variables
    int idx = 0;
    int idu = 0;

    // Concatenate stages in Xopt (9 x N+1)
    Map<MatrixXd> stages(filecontent.topRows(n_states*(N+1)).data(),n_states,N+1);
    Map<MatrixXd> control(filecontent.bottomRows(n_controls*(N+1)).data(),n_controls,N+1);
    MatrixXd Xopt(stages.rows()+control.rows(), stages.cols());
    Xopt << control, 
            stages;
    traj.Xopt = Xopt;

}

// heading: get heading of the middle trajectory 
void TRO::heading(){
    
    double xv, yv;
    
    VectorXd midX = traj.middlePoints.col(0);
    VectorXd midY = traj.middlePoints.col(1);

    VectorXd Pheading(midX.size()-1);

    for(int i = 0; i < midX.size()-1; i++){

        xv = midX(i+1) - midX(i);
        yv = midY(i+1) - midY(i);
        Pheading(i) = atan2(yv,xv);

    }
    traj.Pheading = Pheading;

}

// get_trajectory: get x,y coordinates from calculated MPC stages (xopt) and calculate splines
void TRO::get_trajectory(){

    MatrixXd finalTraj(N,2);
    Map<VectorXd,0,InnerStride<5> > middle_x5(traj.middlePoints.col(0).data(), N);
    Map<VectorXd,0,InnerStride<5> > middle_y5(traj.middlePoints.col(1).data(), N);
    Map<VectorXd,0,InnerStride<5> > Pheading5(traj.Pheading.data(), N);

    // finalTraj.col(0) = traj.middlePoints.col(0) - traj.Xopt.row(4).transpose()*traj.Pheading.array().sin();
    // finalTraj.col(1) = traj.middlePoints.col(1) - traj.Xopt.row(4).transpose()*traj.Pheading.array().cos();

    for(int i = 0; i < N; i++){

        finalTraj(i,0) = middle_x5(i) - traj.Xopt(4,i)*sin(Pheading5(i));
        finalTraj(i,1) = middle_y5(i) + traj.Xopt(4,i)*cos(Pheading5(i));
        
    }
    traj.pointsTraj = finalTraj;
}

// radi_curv: get curvature of trajectory
void TRO::radi_curv(){

}

// Create KDTree for the planner
void TRO::create_KDTree(){

    // Initialize vector
    vector<Point> tr;

    // Fill vector
    for(int i=0; i < traj.pointsTraj.rows(); i++){
        Point p;
        p[0] = traj.pointsTraj(i,0);
        p[1] = traj.pointsTraj(i,1);

        tr.push_back(p);
    }

    // Initialize KDTree
    traj.trajTree.build(tr);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//------------------------Auxiliar functions--------------------------------------------------

// coefs_splines: returns the coefficients of the splines that interpolates all N points
/* Explanation:
        coefs = [ d1 c1 b1 a1
                  d2 c2 b2 a2
                    .....
                  dn cn bn an ]
        where pi(t) = di*t^3 + ci*t^2 + bi*t + ai;
        with 0<t<1 */
MatrixXd TRO::coefs_splines(VectorXd x){

    // Getting the lenght
    int N = x.size();

    // Initializing matrixs
    MatrixXd Deqs = MatrixXd::Zero(4*N,4*N);
    MatrixXd matrixEqs(4,8);
    matrixEqs << 0,  0,  0, 1, 0, 0, 0,  0,
                -1, -2, -3, 0, 1, 0, 0,  0,
                    0, -1, -3, 0, 0, 1, 0,  0,
                    0,  0,  0, 1, 1, 1, 1, -1; // Matrix 4x8

    // Building Deqs matrix

        // Building core:
        for(int i=2; i<=N-1; i++) {
            Deqs.block(4*(i-1),4*(i-1)-3, 4, 8) = matrixEqs;
        }
        // Building corners:
        Deqs.topLeftCorner(4,5) = matrixEqs.rightCols(5);
        Deqs.topRightCorner(4,3) = matrixEqs.leftCols(3);
        Deqs.bottomLeftCorner(4,1) = matrixEqs.rightCols(1);
        Deqs.bottomRightCorner(4,7) = matrixEqs.leftCols(7);
    
        // Creating column of x with zeros
        VectorXd xp = VectorXd::Zero(4*N,1);
        for(int i=0; i<=N-1; i++){
            xp(4*i) = x(i);
        }

        // Obtaining coefficients
        VectorXd coefsVec(4*N,1);
        // MatrixXd Dinv = Deqs.inverse();
        // coefsVec = Dinv*xp;
        coefsVec = Deqs.partialPivLu().solve(xp);

        // Reshape and reorganizing
        Map<MatrixXd> coefsMat1(coefsVec.data(), 4,N); 
        MatrixXd coefsMat2 = coefsMat1.transpose();
        MatrixXd coefsMat = coefsMat2.rowwise().reverse();

    return coefsMat;

}

// Auxiliar function: evaluates a 3rd order polynomial given its parameters and a vector x
// (same as Matlab)
VectorXd TRO::polyval(Vector4d coeffs, VectorXd t) {

    VectorXd results(t.size());
    double result;
    int aux;
    for (int i=0; i<t.size(); i++){
        result = 0;
        aux = 0;
        for (int deg = coeffs.size() - 1; deg >= 0; deg--){
           result += coeffs(aux) * pow(t(i), deg);
           aux++;
        }
        results(i) = result;
    }
    return results;
}

// Same but evaluates polynomial to only one point, so t is a double not a vector
double TRO::polyval2(Vector4d coeffs, double t) {
    return coeffs(0)*pow(t,3) + coeffs(1)*pow(t,2) + coeffs(2)*t + coeffs(3);
}

// integral_length: calculate the length of a 3rd order spline given its coefficients. Rectangular integration.
double TRO::integral_length(Vector4d coefsX, Vector4d coefsY){
    
    // This is the precision of integration
    double diffx = 0.0005;
    double dx, dy, cx, cy, bx, by, t, n, area, function;
    
    dx = coefsX(0);
    cx = coefsX(1);
    bx = coefsX(2);
    dy = coefsY(0);
    cy = coefsY(1);
    by = coefsY(2);
    n = 1/diffx;
    area = 0;

    // Numerical integration with the area of rectangles
    for(int j=0; j<n; j++){
        t = j*diffx;
        function = sqrt(pow(3*dx*pow(t, 2) + 2*cx*t + bx, 2) + pow(3*dy*pow(t, 2) + 2*cy*t + by, 2));
        area += function*diffx;
    }
    return area;
}

MatrixXd TRO::read_csv(const std::string filename, bool firstout){ 

    MatrixXd FileContent; // Eigen Matrix to store csv data
    std::vector<std::vector<double>> result; // Vector of vectors to store each line content

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line;
    int ncols = 0;

    // Read data, line by line
    while(std::getline(myFile, line)){

        // Create a stringstream of the current line
        std::vector<double> linevec;
        std::stringstream ss(line);
        std::string val;

        while(std::getline(ss,val,',')){
            linevec.push_back(std::stod(val));
        }
        result.push_back(linevec);
        ncols = linevec.size();
    }

    int first = 0;
    if(firstout){
        FileContent.resize(result.size(),ncols-1); // Eigen Matrix to store all csv content
        first++;
    }else{
        FileContent.resize(result.size(),ncols); 
    }

    // Store csv data into Eigen Matrix
    long int cols, rows;
    rows = 0;
    vector<vector<double>>::iterator row;
    vector<double>::iterator col;
    for(row = result.begin(); row != result.end(); ++row){
        cols = 0;
        for(col = row->begin() + first; col != row->end(); ++col){
            FileContent(rows,cols) = *col; 
            cols++;
        }
        rows++;
    }

    // Close file
    myFile.close();
    return FileContent;
}

// void printFrequency(string st)
// {
//     // Each word it mapped to
//     // it's frequency
//     map<string, int>FW;
   
//     // Used for breaking words
//     stringstream ss(st);
   
//     // To store individual words
//     string Word;
 
//     while (ss >> Word)
//         FW[Word]++;
 
//     map<string, int>::iterator m;
//     for (m = FW.begin(); m != FW.end(); m++)
//         cout << m->first << "-> "
//              << m->second << "\n";
// }

// template<typename M>
// M read_csv(const std::string & path) {

//     std::ifstream indata;
//     indata.open(path);

//     std::string line;
//     std::vector<double> values;

//     uint rows = 0;
//     while (std::getline(indata, line)) {

//         std::stringstream lineStream(line);
//         std::string cell;

//         while (std::getline(lineStream, cell, ',')) {
//             values.push_back(std::stod(cell));
//         }
//         ++rows;
//     }

//     indata.close();

//     return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>> (values.data(), rows, values.size()/rows);
// }

///////////////////////////////////////////////////////////////////////////////////////////////
//------------------------Planner functions--------------------------------------------------

