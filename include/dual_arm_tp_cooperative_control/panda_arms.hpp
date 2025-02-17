#pragma once

#include <string>

#include <Eigen/Eigen>
#include <controller_interface/controller_interface.hpp>
#include "franka_semantic_components/franka_robot_model.hpp"
#include <rclcpp/rclcpp.hpp>
#include <rclcpp/logger.hpp>
#include <Eigen/Dense>
#define VERSOR_LEMMA_TH 0.000001

namespace dual_arm_tp{
class panda_arm
      {
    public:
      Eigen::Matrix4d wTe, wTb, bTe,wTo,eTt,wTt,wTg,tTo;
      Eigen::Matrix<double,7,1>q,qdot;
      Eigen::Matrix<double,6,7>bJe,bJt,wJt,wJe,wJo;
      Eigen::Matrix<double,7,1>jlmin,jlmax;
      Eigen::Matrix<double,6,1>xdot_tool;
      Eigen::Matrix<double,7,1>xdot_jl;
      Eigen::Matrix<double,6,1>xdot_alt;
      Eigen::Matrix<double,6,6>Ste,wWb,Sto;
      Eigen::Matrix<double,3,1>pos_error,rot_error,pos_error_goal,rot_error_goal;
      Eigen::Matrix<double,6,6>A_ma;
      double alt;
      panda_arm()=default;
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

class mission
{
  public:
  int phase;
  double phase_time;
  std::string prev_action, current_action;
  std::vector<std::string> go_to_tasks;
  std::vector<std::string> coop_manip_tasks;
  std::vector<std::string> end_motion_tasks;

  mission()=default;
  private:

};

class act_functions
{

public:
Eigen::Matrix<double,6,6> A_tool,A_rc;
Eigen::Matrix<double,14,14>A_jl;
act_functions()=default;
private:
  /* data */
};

class task_jacobians
{
  public:
  Eigen::Matrix<double,6,14> RC_jacobian, l_ma_jacobian,r_ma_jacobian,left_tool_jacobian,right_tool_jacobian;
  task_jacobians()=default;

  private:

};

//DUAL ARM COOP TASK AND GOAL DEFINITION
Eigen::Matrix<double,6,1> coop_rc_task;
Eigen::Matrix<double,4,4> wTog;
Eigen::Matrix<double,14,14>J_jl;
double DecreasingBellShapedFunction(double xmin, double xmax, double ymin, double ymax, double x);
double IncreasingBellShapedFunction(double xmin, double xmax, double ymin, double ymax, double x);
Eigen::MatrixXd RegPseudoInverse(Eigen::MatrixXd W,double lambda, double threshold);
Eigen::MatrixXd iCAT_pseudoInverse(Eigen::MatrixXd J,Eigen::MatrixXd A,Eigen::MatrixXd Q,double lambda,double threshold,double weight);
std::tuple<Eigen::MatrixXd,Eigen::MatrixXd> iCAT_task(Eigen::MatrixXd A,Eigen::MatrixXd J,Eigen::MatrixXd Qold,Eigen::VectorXd rhoold,Eigen::VectorXd xdot,double lambda,double threshold,double weight);
Eigen::VectorXd Saturate(Eigen::VectorXd x,double xmax);
Eigen::Vector3d VersorLemma(Eigen::Matrix3d R1, Eigen::Matrix3d R2);
Eigen::Matrix3d rotation(double thx,double thy,double thz);
Eigen::Matrix3d skew(Eigen::Matrix<double,3,1> v);
std::tuple<Eigen::Matrix<double,3,1>,Eigen::Matrix<double,3,1>> cart_error(Eigen::MatrixXd wTgoal,Eigen::MatrixXd wTcp);
std::tuple<panda_arm,panda_arm> Update_transform(panda_arm left_arm,panda_arm right_arm,Eigen::MatrixXd left_bTec, Eigen::MatrixXd right_bTec,mission state);

std::tuple<panda_arm,panda_arm> Compute_Jacobian(panda_arm left_arm,panda_arm right_arm,Eigen::MatrixXd left_bJec,Eigen::MatrixXd right_bJec, mission state);
std::tuple<panda_arm,panda_arm> ComputeTaskReferences(panda_arm left_arm,panda_arm right_arm,mission state);

std::tuple<panda_arm,panda_arm,mission> update_mission_phase(panda_arm left_arm,panda_arm right_arm,mission state);
double ActionTransition(std::string task,std::vector<std::string> prev_tasks, std::vector<std::string> curr_tasks, double time);

std::tuple<panda_arm,panda_arm,act_functions> ComputeActivationFunctions(panda_arm left_arm,panda_arm right_arm,mission state,act_functions A);
// Compute task jacobians
std::tuple<task_jacobians> generate_jacobians(panda_arm left_arm,panda_arm right_arm,mission state,task_jacobians dual_arm_jacob);



}