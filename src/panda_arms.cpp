#include "dual_arm_tp_cooperative_control/panda_arms.hpp"
#include <cstring>
#include <iostream>




namespace dual_arm_tp{

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

double IncreasingBellShapedFunction(double xmin, double xmax, double ymin, double ymax, double x)
{
    double y, cosarg;
    if (x<=xmin)
    {
        y=ymin;
    }
    else if(x>=xmax)
    {
        y=ymax;
    }
    else
    {
        cosarg=((x-xmin)*M_PI/(xmax-xmin))+M_PI;
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
    double costheta=((i1.dot(i2)+j1.dot(j2)+k1.dot(k2))-1)*0.5;
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

std::tuple<panda_arm,panda_arm> Update_transform(panda_arm left_arm,panda_arm right_arm,Eigen::MatrixXd left_bTec, Eigen::MatrixXd right_bTec,mission state)
{
    
    left_arm.bTe<<Eigen::MatrixXd::Identity(4,4);
    left_arm.bTe.block(0,0,3,3)<<left_bTec.block(0,0,3,3);
    left_arm.bTe.block(0,3,3,1)<<left_bTec.block(0,3,3,1);
    right_arm.bTe<<Eigen::MatrixXd::Identity(4,4);
    right_arm.bTe.block(0,0,3,3)<<right_bTec.block(0,0,3,3);
    right_arm.bTe.block(0,3,3,1)<<right_bTec.block(0,3,3,1);
    //right_arm.bTe=right_current_tf;
    //CONVERT FROM B FRAME TO W FRAME
    left_arm.wTe<<left_arm.wTb*left_arm.bTe;
    right_arm.wTe<<right_arm.wTb*right_arm.bTe;
    //CONVERT EE FRAME TO TOOL FRAME
    left_arm.wTt<<left_arm.wTe*left_arm.eTt;
    right_arm.wTt<<right_arm.wTe*right_arm.eTt;

    //DUAL ARM MANIP MISSION PHASE
    if (state.phase==2)
    {
        left_arm.wTo=left_arm.wTt*left_arm.tTo;
        right_arm.wTo=right_arm.wTt*right_arm.tTo;
    }

    
    //left_arm.bTe<<Eigen::MatrixXd::Identity(4,4);
    return std::make_tuple(left_arm,right_arm);
}
std::tuple<panda_arm,panda_arm> Compute_Jacobian(panda_arm left_arm,panda_arm right_arm,Eigen::MatrixXd left_bJec,Eigen::MatrixXd right_bJec, mission state)
{
  //COMPUTE LEFT TOOL JACOBIAN MATRIX
    //Eigen::Matrix<double, 6, 7> left_jacobian(left_sensor->getZeroJacobian(franka::Frame::kEndEffector).data());
    left_arm.bJe<<left_bJec;
    left_arm.Ste=Eigen::MatrixXd::Zero(6,6);
    left_arm.Ste.block(0,0,3,3)<<Eigen::MatrixXd::Identity(3,3);
    left_arm.Ste.block(0,3,3,3)<<-skew(left_arm.wTe.block(0,0,3,3)*left_arm.eTt.block(0,3,3,1));
    left_arm.Ste.block(3,3,3,3)<<Eigen::MatrixXd::Identity(3,3);
    //std::cout <<wWbl<< std::endl;
    left_arm.wJt<<left_arm.Ste*left_arm.wWb*left_arm.bJe;
    //COMPUTE RIGHT TOOL JACOBIAN MATRIX
    //Eigen::Matrix<double, 6, 7> right_jacobian(right_sensor->getZeroJacobian(franka::Frame::kEndEffector).data());
    right_arm.bJe<<right_bJec;
    right_arm.Ste=Eigen::MatrixXd::Zero(6,6);
    right_arm.Ste.block(0,0,3,3)<<Eigen::MatrixXd::Identity(3,3);
    right_arm.Ste.block(0,3,3,3)<<-skew(right_arm.wTe.block(0,0,3,3)*right_arm.eTt.block(0,3,3,1));
    right_arm.Ste.block(3,3,3,3)<<Eigen::MatrixXd::Identity(3,3);
    //wWbr.block(0,0,3,3)<<right_arm.wTb.block(0,0,3,3);
    //wWbr.block(3,3,3,3)<<right_arm.wTb.block(0,0,3,3);
    right_arm.wJt<<right_arm.Ste*right_arm.wWb*right_arm.bJe;

    left_arm.Sto=Eigen::MatrixXd::Zero(6,6);
    right_arm.Sto=Eigen::MatrixXd::Zero(6,6);
    if (state.phase==2)
    {
        
        left_arm.Sto.block(0,0,3,3)<<Eigen::MatrixXd::Identity(3,3);
        left_arm.Sto.block(0,3,3,3)<<-skew(left_arm.wTt.block(0,0,3,3)*left_arm.tTo.block(0,3,3,1));
        left_arm.Sto.block(3,3,3,3)<<Eigen::MatrixXd::Identity(3,3);
        left_arm.wJo=left_arm.Sto*left_arm.wJt;
        
        right_arm.Sto.block(0,0,3,3)<<Eigen::MatrixXd::Identity(3,3);
        right_arm.Sto.block(0,3,3,3)<<-skew(right_arm.wTt.block(0,0,3,3)*right_arm.tTo.block(0,3,3,1));
        right_arm.Sto.block(3,3,3,3)<<Eigen::MatrixXd::Identity(3,3);
        right_arm.wJo=right_arm.Sto*right_arm.wJt;
        
        //left_arm.wTo=left_arm.wTt*left_arm.tTo;
        //right_arm.wTo=right_arm.wTt*right_arm.tTo;
    }

 J_jl=Eigen::MatrixXd::Identity(14,14);



  return std::make_tuple(left_arm,right_arm);
}

std::tuple<panda_arm,panda_arm> ComputeTaskReferences(panda_arm left_arm,panda_arm right_arm,mission state)
{
    //JOINT LIMIT TASK ERROR
    Eigen::MatrixXd jl_median_l=(left_arm.jlmin+left_arm.jlmax)*0.5;
    left_arm.xdot_jl=1*(jl_median_l-left_arm.q);
    Eigen::MatrixXd jl_median_r=(right_arm.jlmin+right_arm.jlmax)*0.5;
    right_arm.xdot_jl=1*(jl_median_r-right_arm.q);

    //MINIMUM ALTITUD TASK

    Eigen::Vector3d kw,tool_left_position,tool_right_position;
    left_arm.xdot_alt=Eigen::MatrixXd::Zero(6,1);
    right_arm.xdot_alt=Eigen::MatrixXd::Zero(6,1);
    kw<<0,0,1;
    tool_left_position<<left_arm.wTt.block(0,3,3,1);
    tool_right_position<<right_arm.wTt.block(0,3,3,1);
    left_arm.alt=std::abs(kw.dot(tool_left_position));
    right_arm.alt=std::abs(kw.dot(tool_right_position));
    left_arm.xdot_alt(2)=0.2*(0.2-left_arm.alt);
    right_arm.xdot_alt(2)=0.2*(0.2-right_arm.alt);
    
    //TOOL POSITION TASK

    switch (state.phase)
    {
    case 1:
    {
        Eigen::Matrix<double,3,1> lin_error_left,ang_error_left;
        std::tie(lin_error_left,ang_error_left)=cart_error(left_arm.wTg,left_arm.wTt);
        left_arm.pos_error<<lin_error_left;
        left_arm.rot_error<<ang_error_left;
    //lin_error_left=left_arm.wTg.block(0,3,3,1)-left_arm.wTt.block(0,3,3,1);
    //ang_error_left=rot_err(left_arm.wTg.block(0,0,3,3),left_arm.wTt.block(0,0,3,3));
        left_arm.xdot_tool<<1.0*lin_error_left,1.0*ang_error_left;
        // left_arm.xdot_tool.block(0,0,3,1)<<Saturate(left_arm.xdot_tool.block(0,0,3,1),0.3);
        // left_arm.xdot_tool.block(3,0,3,1)<<Saturate(left_arm.xdot_tool.block(3,0,3,1),0.3);
        Eigen::Matrix<double,3,1> lin_error_right,ang_error_right; 
        std::tie(lin_error_right,ang_error_right)=cart_error(right_arm.wTg,right_arm.wTt);
    //lin_error_left=left_arm.wTg.block(0,3,3,1)-left_arm.wTt.block(0,3,3,1);   
    //left_arm.rot_error<<ang_error_left;
    //left_arm.bTe<<Eigen::MatrixXd::Identity(4,4);
        right_arm.pos_error<<lin_error_right;
        right_arm.rot_error<<ang_error_right;
        right_arm.xdot_tool<<1.0*lin_error_right,1.0*ang_error_right;
        // right_arm.xdot_tool.block(0,0,3,1)<<Saturate(right_arm.xdot_tool.block(0,0,3,1),0.3);
        // right_arm.xdot_tool.block(3,0,3,1)<<Saturate(right_arm.xdot_tool.block(3,0,3,1),0.3);
        break;
    }
    case 2:
    {
        coop_rc_task=Eigen::MatrixXd::Zero(6,1);
        Eigen::Matrix<double,3,1> lin_error_left,ang_error_left;
        std::tie(lin_error_left,ang_error_left)=cart_error(wTog,left_arm.wTo);
        left_arm.pos_error_goal<<lin_error_left;
        left_arm.rot_error_goal<<ang_error_left;    
        left_arm.xdot_tool<<1.0*lin_error_left,1.0*ang_error_left;
        // left_arm.xdot_tool.block(0,0,3,1)<<Saturate(left_arm.xdot_tool.block(0,0,3,1),0.05);
        // left_arm.xdot_tool.block(3,0,3,1)<<Saturate(left_arm.xdot_tool.block(3,0,3,1),0.05);
        Eigen::Matrix<double,3,1> lin_error_right,ang_error_right; 
        std::tie(lin_error_right,ang_error_right)=cart_error(wTog,right_arm.wTo);
        right_arm.pos_error_goal<<lin_error_right;
        right_arm.rot_error_goal<<ang_error_right;    
        right_arm.xdot_tool<<1.0*lin_error_right,1.0*ang_error_right;
        // right_arm.xdot_tool.block(0,0,3,1)<<Saturate(right_arm.xdot_tool.block(0,0,3,1),0.05);
        // right_arm.xdot_tool.block(3,0,3,1)<<Saturate(right_arm.xdot_tool.block(3,0,3,1),0.05);

        break;
    }

    case 3:
    {
        left_arm.xdot_tool<<Eigen::MatrixXd::Zero(6,1);
        right_arm.xdot_tool<<Eigen::MatrixXd::Zero(6,1);

    }

    }

    

    return std::make_tuple(left_arm,right_arm);

}

std::tuple<panda_arm,panda_arm,mission> update_mission_phase(panda_arm left_arm,panda_arm right_arm,mission state)
{
    Eigen::Matrix<double,3,3> tRol, tRor;
    Eigen::Matrix<double,3,1> tPol, tPor;
    Eigen::Matrix<double,4,4> tTol, tTor;
    tRol=Eigen::MatrixXd::Identity(3,3);
    tRor=Eigen::MatrixXd::Identity(3,3);
    tPol<<0,0,0;
    tPor<<0,0,0;
    tTol=Eigen::MatrixXd::Identity(4,4);
    tTor=Eigen::MatrixXd::Identity(4,4);

    switch (state.phase)
    {
    case 1:


    if(left_arm.pos_error.norm()<0.03 && right_arm.pos_error.norm()<0.03 && left_arm.rot_error.norm()<0.2 && right_arm.rot_error.norm()<0.2)
            {
            state.phase=2;
            state.prev_action=state.current_action;
            state.current_action="coop_manip";
            // ADD OBJ FRAME RIGID BODY TO MANIP CHAIN
            tRol=left_arm.wTt.block(0,0,3,3).transpose()*left_arm.wTo.block(0,0,3,3);
            tPol<<-left_arm.wTt.block(0,3,3,1)+left_arm.wTo.block(0,3,3,1);
            left_arm.tTo<<tRol,tPol,
              Eigen::MatrixXd::Zero(1,3),1;
            tRor=right_arm.wTt.block(0,0,3,3).transpose()*right_arm.wTo.block(0,0,3,3);
            tPor<<-right_arm.wTt.block(0,3,3,1)+right_arm.wTo.block(0,3,3,1);
            right_arm.tTo<<tRor,tPor,
              Eigen::MatrixXd::Zero(1,3),1;

            }
        /* code */
        break;
    
    case 2:
        if(left_arm.pos_error_goal.norm()<0.02 && right_arm.pos_error_goal.norm()<0.3 && left_arm.rot_error_goal.norm()<0.0698 && right_arm.rot_error_goal.norm()<0.0698)
            {
            state.phase=3;
            state.prev_action=state.current_action;
            state.current_action="end_motion";
            state.phase_time=0;
            }
        break;
    case 3:
        break;

    }





    return std::make_tuple(left_arm,right_arm,state);
}

double ActionTransition(std::string task,std::vector<std::string> prev_tasks, std::vector<std::string> curr_tasks, double time)
{
    auto prev = std::find(prev_tasks.begin(),prev_tasks.end(),task);
    auto curr = std::find(curr_tasks.begin(),curr_tasks.end(),task);
    double A=0;
    if(prev!=prev_tasks.end() && curr!=curr_tasks.end() )
    {
        A=1;
        //std::cout <<"A:"<<A<<std::endl;
    }
    else if(prev==prev_tasks.end() && curr!=curr_tasks.end() )
    {
        A=IncreasingBellShapedFunction(0,1,0,1,time);
        //std::cout <<"Increase bell function:"<<std::endl;
    }
    else if(prev!=prev_tasks.end() && curr==curr_tasks.end() )
    {
        //std::cout <<"Decrease bell function:"<<std::endl;
        A=DecreasingBellShapedFunction(0,1,0,1,time);
    }
    else
    {
        A=0;
        //std::cout <<"A:"<<A<<std::endl;
    }
    return A;
}

std::tuple<panda_arm,panda_arm,act_functions> ComputeActivationFunctions(panda_arm left_arm,panda_arm right_arm,mission state,act_functions A)
{
    std::vector<std::string> prev_tasks;
    std::vector<std::string> current_tasks;
    if (state.prev_action=="go_to")
    {
        prev_tasks=state.go_to_tasks;
    }
    else if (state.prev_action=="coop_manip")
    {
        prev_tasks=state.coop_manip_tasks;
    }
    
    if (state.current_action=="go_to")
    {
        current_tasks=state.go_to_tasks;
    }
    else if (state.current_action=="coop_manip")
    {
        current_tasks=state.coop_manip_tasks;
    }

    switch (state.phase)
    {
    case 1:
        A.A_tool=Eigen::MatrixXd::Identity(6,6)*ActionTransition("T",prev_tasks,current_tasks,state.phase_time);
        
        break;
    case 2:
        A.A_rc=Eigen::MatrixXd::Identity(6,6);
        A.A_tool=Eigen::MatrixXd::Identity(6,6)*ActionTransition("T",prev_tasks,current_tasks,state.phase_time);
        
        break;
    case 3:
        A.A_tool=Eigen::MatrixXd::Zero(6,6);
        A.A_rc=Eigen::MatrixXd::Zero(6,6);
        break;
    }
    //Minimum altitude task
    left_arm.A_ma=Eigen::MatrixXd::Identity(6,6)*DecreasingBellShapedFunction(0.15, 0.20,0,1,left_arm.alt)*ActionTransition("MA", prev_tasks, current_tasks,state.phase_time);
    right_arm.A_ma=Eigen::MatrixXd::Identity(6,6)*DecreasingBellShapedFunction(0.15, 0.20,0,1,right_arm.alt)*ActionTransition("MA", prev_tasks, current_tasks,state.phase_time);
    //Joint limit task
    Eigen::MatrixXd jl_delta_l=0.1*(left_arm.jlmax-left_arm.jlmin)*0.5;
    Eigen::MatrixXd jl_delta_r=0.1*(right_arm.jlmax-right_arm.jlmin)*0.5;
    // //LEFT ARM
     Eigen::Matrix<double,7,1> left_A_jl_min,left_A_jl_max,right_A_jl_min,right_A_jl_max;
     left_A_jl_min=Eigen::MatrixXd::Zero(7,1);
     left_A_jl_max=Eigen::MatrixXd::Zero(7,1);
     right_A_jl_min=Eigen::MatrixXd::Zero(7,1);
     right_A_jl_max=Eigen::MatrixXd::Zero(7,1);
    for (int i = 0; i < 7; i++)
    {
        left_A_jl_min(i)=DecreasingBellShapedFunction(left_arm.jlmin(i),left_arm.jlmin(i)+jl_delta_l(i),0,1,left_arm.q(i));
        left_A_jl_max(i)=IncreasingBellShapedFunction(left_arm.jlmax(i),left_arm.jlmax(i)+jl_delta_l(i),0,1,left_arm.q(i));
        right_A_jl_min(i)=DecreasingBellShapedFunction(right_arm.jlmin(i),right_arm.jlmin(i)+jl_delta_r(i),0,1,right_arm.q(i));
        right_A_jl_max(i)=IncreasingBellShapedFunction(right_arm.jlmax(i),right_arm.jlmax(i)+jl_delta_r(i),0,1,right_arm.q(i));
    }
     Eigen::Matrix<double,14,1> Ajl=Eigen::MatrixXd::Zero(14,1);
    Ajl.block(0,0,7,1)<<(left_A_jl_min+left_A_jl_max)*ActionTransition("JL",prev_tasks,current_tasks,state.phase_time);
    Ajl.block(7,0,7,1)<<(right_A_jl_min+right_A_jl_max)*ActionTransition("JL",prev_tasks,current_tasks,state.phase_time);
    A.A_jl=Ajl.array().matrix().asDiagonal();

    return std::make_tuple(left_arm,right_arm,A);
}

std::tuple<task_jacobians> generate_jacobians(panda_arm left_arm,panda_arm right_arm,mission state,task_jacobians dual_arm_jacob)
{
  Eigen::Matrix<double,6,7> left_tool_jac,right_tool_jac;

  if(state.phase==1)
  {
    left_tool_jac<<left_arm.wJt;
    right_tool_jac<<right_arm.wJt;

  }
  else if(state.phase==2)
  {
    left_tool_jac<<left_arm.wJo;
    right_tool_jac<<right_arm.wJo;
  }
    //RIGID CONSTRAINT TASK
    dual_arm_jacob.RC_jacobian.block(0,0,6,7)=left_tool_jac;
    dual_arm_jacob.RC_jacobian.block(0,7,6,7)=-right_tool_jac;
    //MINIMUM ALTITUDE TASK
    dual_arm_jacob.l_ma_jacobian.block(0,0,6,7)=left_tool_jac;
    dual_arm_jacob.l_ma_jacobian.block(0,7,6,7)=Eigen::MatrixXd::Zero(6,7);
    dual_arm_jacob.r_ma_jacobian.block(0,0,6,7)=Eigen::MatrixXd::Zero(6,7);
    dual_arm_jacob.r_ma_jacobian.block(0,7,6,7)=right_tool_jac;
    //TOOL JACOBIAN
    dual_arm_jacob.left_tool_jacobian.block(0,0,6,7)=left_tool_jac;
    dual_arm_jacob.left_tool_jacobian.block(0,7,6,7)=Eigen::MatrixXd::Zero(6,7);
    dual_arm_jacob.right_tool_jacobian.block(0,0,6,7)=Eigen::MatrixXd::Zero(6,7);
    dual_arm_jacob.right_tool_jacobian.block(0,7,6,7)=right_tool_jac;

  return std::make_tuple(dual_arm_jacob);
}


}
