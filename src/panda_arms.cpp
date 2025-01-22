#include "dual_arm_tp_cooperative_control/panda_arms.hpp"
#include <cstring>
#include <iostream>

Eigen::MatrixXd panda_arm::init_Jacobian()
{
    Eigen::Matrix<double,6,7>Jac0=Eigen::MatrixXd::Zero(6,7);
    return Jac0;
}
Eigen::MatrixXd panda_arm::init_const_tf(panda_arm arm)
{
    arm.wWb<<Eigen::MatrixXd::Identity(6,6);
    arm.wWb.block(0,0,3,3)<<arm.wTb.block(0,0,3,3);
    arm.wWb.block(3,3,3,3)<<arm.wTb.block(0,0,3,3);

    return arm.wWb;
}

double DecreasingBellShapedFunction(double xmin, double xmax, double ymin, double ymax, double x)
{
    double y, cosarg;
    if (x<=xmin)
    {
        y=ymax;
    }
    else if(x>=xmax)
    {
        y=ymin;
    }
    else
    {
        cosarg=(x-xmin)*M_PI/(xmax-xmin);
        y=(ymax-ymin)*(0.5*cos(cosarg)+0.5)+ymin;
    }
    return y;
}

Eigen::MatrixXd RegPseudoInverse(Eigen::MatrixXd W,double lambda, double threshold)
{
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(W, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::MatrixXd SigmaPlus(W.cols(),W.rows());
	Eigen::VectorXd eigen = svd.singularValues();
	int r = eigen.size();
	SigmaPlus.setZero();	
	for (int i = 0; i < r; i++)
	{ 
		double p=DecreasingBellShapedFunction(0,threshold,0,lambda,eigen(i));
        SigmaPlus(i,i) = eigen(i)/(eigen(i)*eigen(i) + p);	
	}
	Eigen::MatrixXd V = svd.matrixV();
	Eigen::MatrixXd UT = svd.matrixU().transpose().eval();
	return V*SigmaPlus*UT;
}

Eigen::MatrixXd iCAT_pseudoInverse(Eigen::MatrixXd J,Eigen::MatrixXd A,Eigen::MatrixXd Q,double lambda,double threshold,double weight)
{
    int m=J.rows();
    int n=J.cols();
    Eigen::MatrixXd Rtor=J.transpose()*(Eigen::MatrixXd::Identity(m,m)-A)*A*J;
    Eigen::MatrixXd Rctrl=weight*(Eigen::MatrixXd::Identity(n,n)-Q).transpose()*(Eigen::MatrixXd::Identity(n,n)-Q);
    Eigen::MatrixXd JJpinv=RegPseudoInverse(J.transpose()*A.transpose()*A*J+Rtor+Rctrl,lambda,threshold);
    JJpinv=JJpinv*J.transpose()*A.transpose()*A;
    return JJpinv;
}
std::tuple<Eigen::MatrixXd,Eigen::MatrixXd> iCAT_task(Eigen::MatrixXd A,Eigen::MatrixXd J,Eigen::MatrixXd Qold,Eigen::VectorXd rhoold,Eigen::VectorXd xdot,double lambda,double threshold,double weight)
{
    int n=J.cols();
    if(A.cols()!=J.rows())
    {
        std::cout<<"Error: A and J have different size"<<std::endl;
        //return;
    }
    Eigen::MatrixXd JQpinv1=iCAT_pseudoInverse(J*Qold,A,Qold,lambda,threshold,weight);
    Eigen::MatrixXd JQpinv2=iCAT_pseudoInverse(J*Qold,A,Eigen::MatrixXd::Identity(n,n),lambda,threshold,weight);
    Eigen::MatrixXd Q_new=Qold*(Eigen::MatrixXd::Identity(n,n)-JQpinv2*J*Qold);
    Eigen::MatrixXd W=J*Qold*JQpinv1;
    Eigen::MatrixXd T=Eigen::MatrixXd::Identity(n,n)-Qold*JQpinv2*W*J;
    Eigen::MatrixXd rho_new=T*rhoold+Qold*JQpinv2*W*xdot;
    return std::make_tuple(Q_new,rho_new);
}
Eigen::VectorXd Saturate(Eigen::VectorXd x,double xmax)
 {
    Eigen::VectorXd out;
    int nv=x.size();
    //std::cout <<nv<< std::endl;
    double max=0;
    for (int i=0;i<nv;i++)
    {
        if(std::abs(x(i))>max)
    {
      max=std::abs(x(i));
    }
    }
    if(max>xmax)
    {
        out=(x/max)*xmax;
    }
    else{
        out=x;        
    }
   return out;
 }
 Eigen::Vector3d VersorLemma(Eigen::Matrix3d R1, Eigen::Matrix3d R2)
{
    Eigen::Vector3d i1,j1,k1,i2,k2,j2,out;
    out<<0,0,0;
    i1<<R1.col(0);
    j1<<R1.col(1);
    k1<<R1.col(2);
    i2<<R2.col(0);
    j2<<R2.col(1);
    k2<<R2.col(2);
    Eigen::Vector3d sigma=i1.cross(i2)+j1.cross(j2)+k1.cross(k2);
    Eigen::Vector3d rosintheta=sigma*(0.5);
    double sintheta=rosintheta.norm();
    double costheta=((i1.dot(i2)+j1.dot(j2)+k1.dot(k2))-1)/2;
    //std::cout <<sintheta<< std::endl;
    Eigen::Vector3d maxnormvector;
    if (sintheta>VERSOR_LEMMA_TH)
    {
        double theta=atan2(sintheta,costheta);
        out<<(rosintheta*(theta/sintheta));
    }
    else
    {
        if(costheta>0)
        {
          out<<0,0,0;
        }
        else
        {
          auto h=R1+R2;
        Eigen::Vector3d temp0=h.col(0);
        Eigen::Vector3d temp1=h.col(1);
        Eigen::Vector3d temp2=h.col(2);
        double t0=temp0.norm();
        double t1=temp1.norm();
        double t2=temp2.norm();
        if(t1>t0)
        {
            if(t1>t2){
                maxnormvector=temp1;
            }
            else
            {
                maxnormvector=temp2;
            }
        }
        else
        {
             if(t2>t0){
                maxnormvector=temp2;
            }
            else
            {
                maxnormvector=temp0;
            }
        }

        if(maxnormvector.norm()!=0)
        {
            out<<(maxnormvector*(M_PI/maxnormvector.norm()));
        }
        else
        {
            out<<0,0,0;
        }

        }
    }
    return out;
}
Eigen::Matrix3d rotation(double thx,double thy,double thz)
{
    Eigen::Matrix3d Rx,Ry,Rz,rot;
    Rx<<1,0,0,
        0,cos(thx),-sin(thx),
        0,sin(thx),cos(thx);
    Ry<<cos(thy),0,sin(thy),
                0,1,0,
        -sin(thy),0,cos(thy);
    Rz<<cos(thz),-sin(thz),0,
        sin(thz),cos(thz),0,
        0,0,1;
    rot<<Rz*Ry*Rx;
    return rot;
}
Eigen::Matrix3d skew(Eigen::Matrix<double,3,1> v)
{
    Eigen::Matrix3d skew_v=Eigen::MatrixXd::Zero(3,3);

    skew_v(0,1) = -v(2);
    skew_v(0,2) =  v(1);
    skew_v(1,2) = -v(0);

    skew_v(1,0) =  v(2);
    skew_v(2,0) = -v(1);
    skew_v(2,1) =  v(0);

    return skew_v;
}
std::tuple<Eigen::Matrix<double,3,1>,Eigen::Matrix<double,3,1>> cart_error(Eigen::MatrixXd wTgoal,Eigen::MatrixXd wTcp)
{
    Eigen::Matrix<double,3,1> error_lin,error_ang;
    error_lin<<wTgoal.block(0,3,3,1)-wTcp.block(0,3,3,1);
    error_ang<<(VersorLemma(wTgoal.block(0,0,3,3),wTcp.block(0,0,3,3)))*(-1);
    return std::make_tuple(error_lin,error_ang);
}

std::tuple<panda_arm,panda_arm> Update_transform(panda_arm left_arm,panda_arm right_arm,std::unique_ptr<franka_semantic_components::FrankaRobotModel> &left_sensor,std::unique_ptr<franka_semantic_components::FrankaRobotModel> &right_sensor)
{
    //UPDATE LEFT EE WRT B FRAME
    Eigen::Map<const Eigen::Matrix<double, 4, 4>> left_current_tf(left_sensor->getPoseMatrix(franka::Frame::kEndEffector).data());
    Eigen::Map<const Eigen::Matrix<double, 4, 4>> right_current_tf(right_sensor->getPoseMatrix(franka::Frame::kEndEffector).data());
    left_arm.bTe<<left_current_tf;
    right_arm.bTe<<right_current_tf;
    //CONVERT FROM B FRAME TO W FRAME
    left_arm.wTe<<left_arm.wTb*left_arm.bTe;
    right_arm.wTe<<right_arm.wTb*right_arm.bTe;
    //CONVERT EE FRAME TO TOOL FRAME
    left_arm.wTt<<left_arm.wTe*left_arm.eTt;
    right_arm.wTt<<right_arm.wTe*right_arm.eTt;
    return std::make_tuple(left_arm,right_arm);
}
std::tuple<panda_arm,panda_arm> Compute_Jacobian(panda_arm left_arm,panda_arm right_arm,std::unique_ptr<franka_semantic_components::FrankaRobotModel> &left_sensor,std::unique_ptr<franka_semantic_components::FrankaRobotModel> &right_sensor)
{
  //COMPUTE LEFT TOOL JACOBIAN MATRIX
    Eigen::Matrix<double, 6, 7> left_jacobian(left_sensor->getZeroJacobian(franka::Frame::kEndEffector).data());
    left_arm.bJe<<left_jacobian;
    left_arm.Ste=Eigen::MatrixXd::Zero(6,6);
    left_arm.Ste.block(0,0,3,3)<<Eigen::MatrixXd::Identity(3,3);
    left_arm.Ste.block(0,3,3,3)<<-skew(left_arm.wTe.block(0,0,3,3)*left_arm.eTt.block(0,3,3,1));
    left_arm.Ste.block(3,3,3,3)<<Eigen::MatrixXd::Identity(3,3);
    //std::cout <<wWbl<< std::endl;
    left_arm.wJt<<left_arm.Ste*left_arm.wWb*left_arm.bJe;
    //COMPUTE RIGHT TOOL JACOBIAN MATRIX
    Eigen::Matrix<double, 6, 7> right_jacobian(right_sensor->getZeroJacobian(franka::Frame::kEndEffector).data());
    right_arm.bJe<<right_jacobian;
    right_arm.Ste=Eigen::MatrixXd::Zero(6,6);
    right_arm.Ste.block(0,0,3,3)<<Eigen::MatrixXd::Identity(3,3);
    right_arm.Ste.block(0,3,3,3)<<-skew(right_arm.wTe.block(0,0,3,3)*right_arm.eTt.block(0,3,3,1));
    right_arm.Ste.block(3,3,3,3)<<Eigen::MatrixXd::Identity(3,3);
    //wWbr.block(0,0,3,3)<<right_arm.wTb.block(0,0,3,3);
    //wWbr.block(3,3,3,3)<<right_arm.wTb.block(0,0,3,3);
    right_arm.wJt<<right_arm.Ste*right_arm.wWb*right_arm.bJe;
  return std::make_tuple(left_arm,right_arm);
}
std::tuple<panda_arm,panda_arm> ComputeTaskReferences(panda_arm left_arm,panda_arm right_arm)
{
  //LEFT ARM
    Eigen::Matrix<double,3,1> lin_error_left,ang_error_left;    
    std::tie(lin_error_left,ang_error_left)=cart_error(left_arm.wTg,left_arm.wTt);
    left_arm.pos_error<<lin_error_left;
    left_arm.rot_error<<ang_error_left;
    left_arm.xdot_tool<<1.0*lin_error_left,1.0*ang_error_left;
    left_arm.xdot_tool.block(0,0,3,1)<<Saturate(left_arm.xdot_tool.block(0,0,3,1),0.1);
    //left_arm.xdot_tool.block(3,0,3,1)<<Saturate(left_arm.xdot_tool.block(3,0,3,1),0.1);
    //RIGHT ARM
    Eigen::Matrix<double,3,1> lin_error_right,ang_error_right;    
    std::tie(lin_error_right,ang_error_right)=cart_error(right_arm.wTg,right_arm.wTt);
    right_arm.pos_error<<lin_error_right;
    right_arm.rot_error<<ang_error_right;
    right_arm.xdot_tool<<1.0*lin_error_right,1.0*ang_error_right;
    right_arm.xdot_tool.block(0,0,3,1)<<Saturate(right_arm.xdot_tool.block(0,0,3,1),0.1);
    //right_arm.xdot_tool.block(3,0,3,1)<<Saturate(right_arm.xdot_tool.block(3,0,3,1),0.1);
    return std::make_tuple(left_arm,right_arm);
}