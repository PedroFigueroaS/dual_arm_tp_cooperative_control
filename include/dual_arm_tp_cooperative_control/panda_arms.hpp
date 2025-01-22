#pragma once

#include <string>

#include <Eigen/Eigen>
#include <controller_interface/controller_interface.hpp>
#include "franka_semantic_components/franka_robot_model.hpp"
#include <rclcpp/rclcpp.hpp>
#include <Eigen/Dense>
#define VERSOR_LEMMA_TH 0.000001


class panda_arm
      {
    public:
      Eigen::Matrix4d wTe, wTb, bTe,wTo,eTt,wTt,wTg;
      Eigen::Matrix<double,7,1>q,qdot;
      Eigen::Matrix<double,6,7>bJe,bJt,wJt,wJe;
      Eigen::Matrix<double,7,1>jlmin,jlmax;
      Eigen::Matrix<double,6,1>xdot_tool;
      Eigen::Matrix<double,6,6>Ste,wWb;
      Eigen::Matrix<double,3,1>pos_error,rot_error;
      panda_arm(Eigen::Matrix4d init_pos,Eigen::Matrix4d offset)
      {
        bTe=init_pos;
        wTb=offset;
        wTe=offset*init_pos;
      }
      Eigen::MatrixXd init_Jacobian();
      Eigen::MatrixXd init_const_tf(panda_arm obj);
    private:
    /* data */
      };
      //END PANDA ARM CLASS

double DecreasingBellShapedFunction(double xmin, double xmax, double ymin, double ymax, double x);
Eigen::MatrixXd RegPseudoInverse(Eigen::MatrixXd W,double lambda, double threshold);
Eigen::MatrixXd iCAT_pseudoInverse(Eigen::MatrixXd J,Eigen::MatrixXd A,Eigen::MatrixXd Q,double lambda,double threshold,double weight);
std::tuple<Eigen::MatrixXd,Eigen::MatrixXd> iCAT_task(Eigen::MatrixXd A,Eigen::MatrixXd J,Eigen::MatrixXd Qold,Eigen::VectorXd rhoold,Eigen::VectorXd xdot,double lambda,double threshold,double weight);
Eigen::VectorXd Saturate(Eigen::VectorXd x,double xmax);
Eigen::Vector3d VersorLemma(Eigen::Matrix3d R1, Eigen::Matrix3d R2);
Eigen::Matrix3d rotation(double thx,double thy,double thz);
Eigen::Matrix3d skew(Eigen::Matrix<double,3,1> v);
std::tuple<Eigen::Matrix<double,3,1>,Eigen::Matrix<double,3,1>> cart_error(Eigen::MatrixXd wTgoal,Eigen::MatrixXd wTcp);
std::tuple<panda_arm,panda_arm> Update_transform(panda_arm left_arm,panda_arm right_arm,std::unique_ptr<franka_semantic_components::FrankaRobotModel> &left_sensor,std::unique_ptr<franka_semantic_components::FrankaRobotModel> &right_sensor);
std::tuple<panda_arm,panda_arm> Compute_Jacobian(panda_arm left_arm,panda_arm right_arm,std::unique_ptr<franka_semantic_components::FrankaRobotModel> &left_sensor,std::unique_ptr<franka_semantic_components::FrankaRobotModel> &right_sensor);
std::tuple<panda_arm,panda_arm> ComputeTaskReferences(panda_arm left_arm,panda_arm right_arm);