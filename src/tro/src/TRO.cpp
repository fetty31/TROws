#include "../include/TRO.hh"

using namespace std;
using namespace Eigen;

// Constructor
TRO::TRO(){}

// Initialization of TRO
void TRO::init(){

    auto start_time = std::chrono::system_clock::now();

    // Run principal functions
    get_data();
    heading();
    get_trajectory();
    radi_curv();
    create_KDTree();

    // Running flag set to true
    isRun = true;

    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    ROS_WARN( "TRO ELAPSED TIME: %f ms", elapsed.count()*1000);
    // cout << elapsed.count() << endl;

}

// Return flag, useful for planner
bool TRO::isRunning(){
    return this->isRun;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//------------------------Principal functions--------------------------------------------------

// get_data: read TRO.py output and transform the data to desired format
void TRO::get_data(){

    MatrixXd filecontent = read_csv(x_opt); 

    MatrixXd middlePoints = read_csv(Pmiddle);
    traj.middlePoints = middlePoints;

    MatrixXd Problem = read_csv(problem);

    MatrixXd FreeSpace = read_csv(freeSpace);
    traj.centralFree = FreeSpace;

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

// get_trajectory: get x,y coordinates from calculated MPC stages (xopt) and interpolate splines
void TRO::get_trajectory(){

    traj.dimension = N/5;

    MatrixXd finalTraj(N,2);
    VectorXd FreeR(N), FreeL(N);
    Map<VectorXd,0,InnerStride<5> > middle_x5(traj.middlePoints.col(0).data(), N);
    Map<VectorXd,0,InnerStride<5> > middle_y5(traj.middlePoints.col(1).data(), N);
    Map<VectorXd,0,InnerStride<5> > freeR(traj.centralFree.col(0).data(), N);
    Map<VectorXd,0,InnerStride<5> > freeL(traj.centralFree.col(1).data(), N);
    Map<VectorXd,0,InnerStride<5> > Pheading5(traj.Pheading.data(), N);
    
    // finalTraj.col(0) = traj.middlePoints.col(0) - traj.Xopt.row(4).transpose()*traj.Pheading.array().sin();
    // finalTraj.col(1) = traj.middlePoints.col(1) - traj.Xopt.row(4).transpose()*traj.Pheading.array().cos();

    for(int i = 0; i < N; i++){

        finalTraj(i,0) = middle_x5(i) - traj.Xopt(4,i)*sin(Pheading5(i));
        finalTraj(i,1) = middle_y5(i) + traj.Xopt(4,i)*cos(Pheading5(i));
        FreeR(i) = traj.Xopt(4,i) + freeR(i);
        FreeL(i) = -traj.Xopt(4,i) + freeL(i);
        
    }
    traj.pointsSol = finalTraj; // Optimized trajectory with N points 
    traj.freeL = FreeL;
    traj.freeR = FreeR;
    Map<MatrixXd,0,InnerStride<5> > pointsSol5(traj.pointsSol.data(),traj.dimension,2); // Reduce number of points in order to compute splines faster

    // Getting the splines coefficients of the optimized trajectory
    traj.coefsSplinesTrajX = coefs_splines(pointsSol5.col(0));
    traj.coefsSplinesTrajY = coefs_splines(pointsSol5.col(1));

    // Getting the length of every spline of the trajectory
    traj.splinesLengths = VectorXd(traj.dimension);
    for (int i=0; i<traj.dimension; i++){
        traj.splinesLengths(i) = integral_length(traj.coefsSplinesTrajX.row(i), traj.coefsSplinesTrajY.row(i));
    }

}

// radi_curv: get curvature of trajectory
void TRO::radi_curv(){

    // Declare variables 
    int numberPointsPerSpline, totalPoints=0;
    double curvature, dx,cx,bx,dy,cy,by;

    // This is the aproximated dimension that will have the vectors
    int aproxL = (int) (traj.splinesLengths.sum()/(this->spacing)) + traj.dimension; 

    // Initializing points, radius, fL, fR, Vx, Vy, w
    MatrixXd points(aproxL, 2);
    VectorXd radius(aproxL);
    VectorXd fL(aproxL), fR(aproxL), VX(aproxL), VY(aproxL), W(aproxL);

    // Setting freeR, freeL to have size = traj.dimension
    Map<VectorXd,0,InnerStride<5> > freeL5(traj.freeL.data(), traj.dimension);
    Map<VectorXd,0,InnerStride<5> > freeR5(traj.freeR.data(), traj.dimension);

    // Setting optimized stages for MPC to have size = traj.dimension
    VectorXd vx = traj.Xopt.row(6);
    VectorXd vy = traj.Xopt.row(7);
    VectorXd ww = traj.Xopt.row(8);
    Map<VectorXd,0,InnerStride<5> > Vx(vx.data(), traj.dimension);
    Map<VectorXd,0,InnerStride<5> > Vy(vy.data(), traj.dimension);
    Map<VectorXd,0,InnerStride<5> > w(ww.data(), traj.dimension);

    for(int i=0; i < traj.dimension; i++){

        // Getting the coefficients for the formula
        dx = traj.coefsSplinesTrajX(i,0);
        cx = traj.coefsSplinesTrajX(i,1);
        bx = traj.coefsSplinesTrajX(i,2);
        dy = traj.coefsSplinesTrajY(i,0);
        cy = traj.coefsSplinesTrajY(i,1);
        by = traj.coefsSplinesTrajY(i,2);

        // Number of points per spline
        numberPointsPerSpline = (int) (traj.splinesLengths(i) / (this->spacing));

        // Parametrization of the spline, numberPointsPerSpline between 0 and almost 1:
        VectorXd t = VectorXd::LinSpaced(numberPointsPerSpline, 0, 1-1/(double)numberPointsPerSpline);

        // Specific coordinates for particular points (the third is the spline index)
        points.block(totalPoints, 0, numberPointsPerSpline, 1) = polyval(traj.coefsSplinesTrajX.row(i), t);
        points.block(totalPoints, 1, numberPointsPerSpline, 1) = polyval(traj.coefsSplinesTrajY.row(i), t);

        // Curvature radius formula:
        for(int j=0; j<numberPointsPerSpline; j++){

            radius(totalPoints + j) = pow(sqrt(pow(3*dx*pow(t(j),2)+2*cx*t(j)+bx, 2) + pow(3*dy*pow(t(j),2)+2*cy*t(j)+by, 2)), 3) / ((3*dx*pow(t(j),2)+2*cx*t(j)+bx)*(6*dy*t(j)+2*cy) - (3*dy*pow(t(j),2)+2*cy*t(j)+by)*(6*dx*t(j)+2*cx));

            // We suppose freeR, freeL stay constant over points in the same spline (polynomial)
            fL(totalPoints + j) = freeL5(i);
            fR(totalPoints + j) = freeR5(i);

            // We suppose Vx,Vy,w stay constant over points in the same spline (polynomial)
            VX(totalPoints + j) = Vx(i);
            VY(totalPoints + j) = Vy(i);
            W(totalPoints + j) = w(i);
        }

        totalPoints += numberPointsPerSpline;
    }

    // Now that it is known the exact value of the size of the vector, they are resized (last empty places are removed)
    traj.radiCurv = radius.head(totalPoints);
    traj.pointsTraj = points.topRows(totalPoints);
    traj.freeL = fL.head(totalPoints); traj.freeR = fR.head(totalPoints);
    traj.Vx = VX.head(totalPoints); traj.Vy = VY.head(totalPoints); traj.w = W.head(totalPoints);

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

// read_csv files (firstout argument is used to avoid first column of indexes)
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

// Main function of the planner: creates the message and returns it.
dv_msgs::ObjectiveArrayCurv TRO::plannerTRO(const dv_msgs::CarState::ConstPtr &data){

    // Create message
    dv_msgs::ObjectiveArrayCurv msg;
    msg.header = data->header;

    // Planning
    MatrixXd plan = planning(data);

    // Fill the message
    for(int i=0; i < steps; i++){
        dv_msgs::ObjectiveCurv objective = dv_msgs::ObjectiveCurv();

        objective.x = plan(i,0);
        objective.y = plan(i,1);
        objective.s = plan(i,2);
        objective.k = plan(i,3);
        objective.vx = plan(i,4);
        objective.vy = plan(i,5);
        objective.w = plan(i,6);
        objective.L = plan(i,7);
        objective.R = plan(i,8);
        
        msg.objectives.push_back(objective);
    }

    msg.smax = traj.pointsTraj.rows() * spacing;

    return msg;
}


// Auxiliar function of the planner: creates the plan matrix of the MPC-Curv dimensions.
MatrixXd TRO::planning(const dv_msgs::CarState::ConstPtr &data){
    
    // Actual position of the car
    Point state;
    state[0] = data->ekf.position.x;
    state[1] = data->ekf.position.y;

    // Search the nearest point of the trajectory
    int nnid = traj.trajTree.nnSearch(state);

    // Initialize planning matrix
    MatrixXd plan(steps, 9); // [x, y, s, k, vx, vy, w, L, R]

    // All steps
    for(int i=0; i < steps; i++){
        if(nnid >= traj.pointsTraj.rows()) nnid -= traj.pointsTraj.rows(); // restart loop 

        plan(i,0) = traj.pointsTraj(nnid,0);
        plan(i,1) = traj.pointsTraj(nnid,1);
        plan(i,2) = nnid*spacing;
        plan(i,3) = 1/traj.radiCurv(nnid);
        plan(i,4) = traj.Vx(nnid);
        plan(i,5) = traj.Vy(nnid);
        plan(i,6) = traj.w(nnid);
        plan(i,7) = traj.freeL(nnid);
        plan(i,8) = traj.freeR(nnid);

        nnid++;
    }

    return plan;
}

// Visualization function that creates the message to plot the trajectory in RVIZ
nav_msgs::Path TRO::get_path(){

    // Init message:
    nav_msgs::Path pathMsg;
	geometry_msgs::PoseStamped pose;

    pathMsg.header.stamp = ros::Time::now();
    pathMsg.header.frame_id = "global";

    for(unsigned int i = 0; i < traj.pointsTraj.rows(); i++ ){

        pose.pose.position.x = traj.pointsTraj(i, 0);
        pose.pose.position.y = traj.pointsTraj(i, 1);

        pathMsg.poses.push_back(pose);
        
    }

    return pathMsg;
}