#include <ros/ros.h>
#include "../include/TRO.hh"

TRO TROobject;

// TRO is initialized here
int main(int argc, char **argv){

    ros::init(argc, argv, "tro");

    ros::NodeHandle nh;

    TROobject.init();
    cout << "tot guai" <<endl;

    // ros::spin();

    return 0;
};
