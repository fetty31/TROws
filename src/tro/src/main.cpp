#include <ros/ros.h>
#include "../include/TRO.hh"
// #include <dynamic_reconfigure/server.h>

// Publishers are initialized here
ros::Publisher troPath;
ros::Publisher troPlanner;

// TRO is initialized here
TRO TROobject;

// This is the inital callback
void callback_init(const dv_msgs::ConeArrayOrdered::ConstPtr &data){

    if(not TROobject.isRunning()){

        // Initialization
        TROobject.init();
    
        // // Visualization of the path
        troPath.publish(TROobject.get_path());

}
}

// This is the planner callback
void callback_planner(const dv_msgs::CarState::ConstPtr &data){

    // Don't do anything until the TRO is initialized
    if (TROobject.isRunning()){

        dv_msgs::ObjectiveArrayCurv msgCurv = TROobject.plannerTRO(data);
        troPlanner.publish(msgCurv);

    }
}


int main(int argc, char **argv){

    // Init Node:
    ros::init(argc, argv, "tro");

    // Handle Connections:
    ros::NodeHandle nh;

    // Publisher & Subscriber:
    ros::Subscriber subMap = nh.subscribe("/cones/loop", 1, callback_init);
    ros::Subscriber subState = nh.subscribe("/state/car", 1, callback_planner);
    troPath = nh.advertise<nav_msgs::Path>("/gro/path", 1); // Visualization purposes only
    troPlanner = nh.advertise<dv_msgs::ObjectiveArrayCurv>("/gro/planner/curv", 1);

    TROobject.init();
    cout << "tot ok" <<endl;

    ros::spin();

}

