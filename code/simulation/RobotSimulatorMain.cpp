#ifdef B3_USE_ROBOTSIM_GUI
#include "b3RobotSimulatorClientAPI.h"
#else
#include "b3RobotSimulatorClientAPI_NoGUI.h"
#endif
#include <chrono>
using namespace std::chrono;

#include "../Utils/b3Clock.h"

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <ctime>
#include <assert.h>
#define ASSERT_EQ(a, b) assert((a) == (b));
#include "MinitaurSetup.h"
#include <random>


#include "../CommonInterfaces/CommonExampleInterface.h"
#include "../CommonInterfaces/CommonGUIHelperInterface.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h"

#include "btBulletDynamicsCommon.h"

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "../CommonInterfaces/CommonRigidBodyBase.h"
#include "LinearMath/btTransform.h"
#include "LinearMath/btHashMap.h"
using namespace std;
//planar 2D
#include "BulletCollision/CollisionShapes/btBox2dShape.h"
#include "BulletCollision/CollisionDispatch/btEmptyCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btBox2dBox2dCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btConvex2dConvex2dAlgorithm.h"

#include "BulletCollision/CollisionShapes/btBox2dShape.h"
#include "BulletCollision/CollisionShapes/btConvex2dShape.h"
#include "BulletCollision/NarrowPhaseCollision/btMinkowskiPenetrationDepthSolver.h"
#define _USE_MATH_DEFINES

#include <cmath>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int main(int argc, char* argv[])
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
	srand(time(0));

	auto start_run = high_resolution_clock::now();
#ifdef B3_USE_ROBOTSIM_GUI
	b3RobotSimulatorClientAPI* sim = new b3RobotSimulatorClientAPI();
	bool isConnected = sim->connect(eCONNECT_GUI);
#else
	b3RobotSimulatorClientAPI_NoGUI* sim = new b3RobotSimulatorClientAPI_NoGUI();
	bool isConnected = sim->connect(eCONNECT_DIRECT);
#endif
	if (!isConnected)
	{
		printf("Cannot connect\n");
		return -1;
	}
	//Can also use eCONNECT_DIRECT,eCONNECT_SHARED_MEMORY,eCONNECT_UDP,eCONNECT_TCP, for example:
	//sim->connect(eCONNECT_UDP, "localhost", 1234);
	sim->configureDebugVisualizer(COV_ENABLE_GUI, 0);
	//	sim->configureDebugVisualizer( COV_ENABLE_SHADOWS, 0);//COV_ENABLE_WIREFRAME
	sim->setTimeOut(10);
	//syncBodies is only needed when connecting to an existing physics server that has already some bodies
	sim->syncBodies();
	btScalar fixedTimeStep = 1. / 240.;

	sim->setTimeStep(fixedTimeStep);

	btQuaternion q = sim->getQuaternionFromEuler(btVector3(0.1, 0.2, 0.3));
	btVector3 rpy;
	rpy = sim->getEulerFromQuaternion(q);


	// sim->setGravity(btVector3(0, 0, -10));
	sim->setGravity(btVector3(0, 0, -9.8));


	//int blockId = sim->loadURDF("cube.urdf");
	//b3BodyInfo bodyInfo;
	//sim->getBodyInfo(blockId,&bodyInfo);

	int plane = sim->loadURDF("./blockDesign/plane.urdf");
	int ang = 90;
	// MinitaurSetup minitaur;
	// int minitaurUid = minitaur.setupMinitaur(sim, btVector3(0, 0, .3));

	std::normal_distribution<double> pos_distribution(0.0, 0.1);
	std::normal_distribution<double> orn_distribution(5*M_PI/180, 12.5*M_PI/180);


mylabel2:
// // ################### 3D block ###################
//
// 	stringstream file5;
//
//   file5 << "./blockDesign/peg_3d_2.urdf";
//   string block_file3 = file5.str();
//
//   stringstream file6;
//
//   file6 << "./blockDesign/socket_3d2.urdf";
//   string base_file3 = file6.str();
//
//   b3RobotSimulatorLoadUrdfFileArgs base_args9;
//
//   base_args9.m_startPosition.setValue(1.5,0,0.0);
//   base_args9.m_startOrientation.setEulerZYX(M_PI, 0, 0);
//
//   int base_underwater2 = sim->loadURDF(base_file3,base_args9);
//
//   b3RobotSimulatorLoadUrdfFileArgs block_arg2;
//
//   block_arg2.m_startPosition.setValue(0,0.6,0.13);
//   block_arg2.m_startOrientation.setEulerZYX(0.1,0,0);
//
//
//   int *block = new int;
//   *block = sim->loadURDF(block_file3,block_arg2);
//
//   int *block2 = new int;
//   int *block3 = new int;



// ################2D block ###############
  // stringstream file1;
  // 
  // file1 << "./blockDesign/peg_1.urdf";
	// string block_file = file1.str();
  // 
  // stringstream file3;
  // 
  // file3 << "./blockDesign/peg_2.urdf";
  // string block_file2 = file3.str();
  // 
	// stringstream file5;
  // 
  // file5 << "./blockDesign/peg_3.urdf";
  // string block_file3 = file5.str();
  // 
  // stringstream file2;
  // 
  // file2 << "./blockDesign/socket_1.urdf";
	// string base_file = file2.str();
  // 
  // stringstream file4;
  // 
  // file4 << "./blockDesign/socket_2.urdf";
  // string base_file2 = file4.str();
  // 
  // stringstream file6;
  // 
  // file6 << "./blockDesign/socket_3.urdf";
  // string base_file3 = file6.str();

  stringstream file;
  // file << "./blockDesign/block_design_anchor.urdf";
  // file << "./blockDesign/socket_1.urdf";
  file << "./blockDesign/block_anchor.urdf";
  string block_file = file.str();
  
  stringstream file_others;
  // file_others << "./blockDesign/block_design.urdf";
  file_others << "./blockDesign/block.urdf";
  string block_file_other = file_others.str();
  
  //
	// b3RobotSimulatorLoadUrdfFileArgs base_args1;
	// base_args1.m_startPosition.setValue(0,1.0,0.0);
	// base_args1.m_startOrientation.setEulerZYX(0, 0, 1*M_PI/4);
  //
	// int base1 = sim->loadURDF(base_file,base_args1);
  //
	// b3RobotSimulatorLoadUrdfFileArgs base_args2;
	// base_args2.m_startPosition.setValue(0.1284,1.0,0.0);
	// base_args2.m_startOrientation.setEulerZYX(0, 0, 1*M_PI/4);
  //
	// int base2 = sim->loadURDF(base_file,base_args2);
  //
	// b3RobotSimulatorLoadUrdfFileArgs base_args3;
	// base_args3.m_startPosition.setValue(0.2568,1.0,0.0);
	// base_args3.m_startOrientation.setEulerZYX(0, 0, 1*M_PI/4);
  //
	// int base3 = sim->loadURDF(base_file,base_args3);
  //
	// b3RobotSimulatorLoadUrdfFileArgs base_args4;
	// base_args4.m_startPosition.setValue(0.3852,1.0,0.0);
	// base_args4.m_startOrientation.setEulerZYX(0, 0, 1*M_PI/4);
  //
	// int base4 = sim->loadURDF(base_file,base_args4);
  //
	// b3RobotSimulatorLoadUrdfFileArgs base_args5;
	// base_args5.m_startPosition.setValue(-0.1284,1.0,0.0);
	// base_args5.m_startOrientation.setEulerZYX(0, 0, 1*M_PI/4);
  //
	// int base5 = sim->loadURDF(base_file,base_args5);
  //
	// b3RobotSimulatorLoadUrdfFileArgs base_args6;
	// base_args6.m_startPosition.setValue(-0.2568,1.0,0.0);
	// base_args6.m_startOrientation.setEulerZYX(0, 0, 1*M_PI/4);
  //
	// int base6 = sim->loadURDF(base_file,base_args6);
  //
	// b3RobotSimulatorLoadUrdfFileArgs base_args7;
	// base_args7.m_startPosition.setValue(-0.3852,1.0,0.0);
	// base_args7.m_startOrientation.setEulerZYX(0, 0, 1*M_PI/4);
  //
	// int base7 = sim->loadURDF(base_file,base_args7);
	// 
	// 
	// 
	// 
	// HERE Loading URDF

  // b3RobotSimulatorLoadUrdfFileArgs base_args1;
  // base_args1.m_startPosition.setValue(-1.25,0,0.0);
  // base_args1.m_startOrientation.setEulerZYX(0, 0, 0);
  // //
  // int base_underwater3 = sim->loadURDF(base_file2,base_args1);
  // 
  // //
  // b3RobotSimulatorLoadUrdfFileArgs base_args9;
  // base_args9.m_startPosition.setValue(0,0,0.0);
  // base_args9.m_startOrientation.setEulerZYX(0, 0, 0);
  // // base_args9.m_startPosition.setValue(1.5,0,0.0);
  // // base_args9.m_startOrientation.setEulerZYX(M_PI, 0, 0);
  // 
  // int base_underwater2 = sim->loadURDF(base_file3,base_args9);
  // 
  // b3RobotSimulatorLoadUrdfFileArgs base_args8;
  // base_args8.m_startPosition.setValue(1.25,0,0.0);
  // base_args8.m_startOrientation.setEulerZYX(0, 0, 0);
  // 
  // int base_underwater = sim->loadURDF(base_file,base_args8);

  b3RobotSimulatorLoadUrdfFileArgs block_arg1;
  block_arg1.m_startPosition.setValue(0,0,2);
  block_arg1.m_startOrientation.setEulerZYX(0,0,0);
  
  int *block = new int;
  *block = sim->loadURDF(block_file, block_arg1);


  // b3RobotSimulatorLoadUrdfFileArgs block_arg1;
  // block_arg1.m_startPosition.setValue(0.9,0.8,0.01);
  // block_arg1.m_startOrientation.setEulerZYX(0.34,0,0);
  // 
  // 
  // int *block2 = new int;
  // *block2 = sim->loadURDF(block_file,block_arg1);
  // 
  // b3RobotSimulatorLoadUrdfFileArgs block_arg2;
  // block_arg2.m_startPosition.setValue(-0.35,0.8,0.01);
  // block_arg2.m_startOrientation.setEulerZYX(0.32,0,0);

  // block_arg2.m_startPosition.setValue(0,0.6,0.13);
  // block_arg2.m_startOrientation.setEulerZYX(0.1,0,0);
  //
  // block_arg2.m_startPosition.setValue(-0.05,0.6,0.01);
  // block_arg2.m_startOrientation.setEulerZYX(0.2,0,0);

  // int *block = new int;
  // *block = sim->loadURDF(block_file3,block_arg2);
  // 
  // b3RobotSimulatorLoadUrdfFileArgs block_arg3;
  // block_arg3.m_startPosition.setValue(-1.6,0.8,0.01);
  // block_arg3.m_startOrientation.setEulerZYX(0.34,0,0);
  // 
  // 
  // int *block3 = new int;
  // *block3 = sim->loadURDF(block_file2,block_arg3);








	btVector3 basePos{0, 0, 60};

	// sim->resetDebugVisualizerCamera(3, -89,20, basePos);
  sim->resetDebugVisualizerCamera(3, -60, 0, basePos);

	int offset_sign = 1;
	int n = 13;
	int k = 1;
	int suc = 0;
	int loop = 1;
	float m = -30;
	float offset_x = 0.0642+-3*0.1284;

	btQuaternion orn;
	btQuaternion quat;

	const btScalar RADIANS_PER_DEGREE = M_PI / btScalar(180.0);

	btTransform tr;
	tr.setIdentity();


	int count2 = 0;


	mylabel:
	count2 ++;
	start_run = high_resolution_clock::now();
	double pos_noise = pos_distribution(generator);
	double orn_noise = orn_distribution(generator);
	orn.setEulerZYX(0, 0, 0);

	int sign_orn = sgn((rand()%101 - 50));

	quat.setEulerZYX(m*RADIANS_PER_DEGREE, 0, 0);
	

	b3Clock clock;
	double startTime = clock.getTimeInSeconds();
	double simWallClockSeconds = 20.;
#if 0
	while (clock.getTimeInSeconds()-startTime < simWallClockSeconds)
	{
		sim->stepSimulation();
	}
#endif
	sim->setRealTimeSimulation(false);
	int vidLogId = -1;
	int ballLogId = -1;
	int rotateCamera = 0;


	stringstream filename1;
	filename1 << "data.dat";
	string file_data = filename1.str();

	stringstream filename2;
	filename2 << "pos.dat";
	string file_pos = filename2.str();

	stringstream filename3;
	filename3 << "orn.dat";
	string file_orn = filename3.str();

	stringstream filename4;
	filename4 << "vel.dat";
	string file_vel = filename4.str();

  stringstream filename5;
  filename5 << "time.dat";
  string file_time = filename5.str();

	ofstream myfile1 (file_data, ios_base::app);
	ofstream myfile2 (file_pos, ios_base::app);
	ofstream myfile3 (file_orn, ios_base::app);
	ofstream myfile4 (file_vel, ios_base::app);
  ofstream myfile5 (file_time, ios_base::app);

	if (k != 1)
	{
		myfile1 <<"\n";
		myfile2 <<"\n";
		myfile3 <<"\n";
		myfile4 <<"\n";
    myfile5 <<"\n";
	}



	float x;
	float y;
	float z;

	float frame_x;
	float frame_y;
	float frame_z;

	btVector3 basePos_ori;
	btQuaternion baseOrn_ori;
	sim->getBasePositionAndOrientation(*block, basePos_ori, baseOrn_ori);

	btVector3 baseLinearVelocity;
	btVector3 baseAngularVelocity;

	const b3RobotSimulatorClientAPI_NoGUI* sim1 = sim;

	sim1->getBaseVelocity(*block, baseLinearVelocity, baseAngularVelocity);



	btVector3 rpy1;
	rpy1 = sim->getEulerFromQuaternion(baseOrn_ori);




	x = basePos_ori.x();
	y = basePos_ori.y();
	z = basePos_ori.z();
	frame_x = basePos_ori.x();
	frame_y = basePos_ori.y();
	frame_z = basePos_ori.z();
	printf("Pos[%d] = [%f,%f,%f], Orn[%d] = [%f,%f,%f]\n", k, basePos_ori.x(), basePos_ori.y(), basePos_ori.z(), k, baseOrn_ori.x()
	,baseOrn_ori.y(), baseOrn_ori.z());

	printf("LinearVelocity[%d] = [%f,%f,%f], AngularVelocity[%d] = [%f,%f,%f]\n", k, baseLinearVelocity[0], baseLinearVelocity[1], baseLinearVelocity[2], k, baseAngularVelocity[0]
	,baseAngularVelocity[1], baseAngularVelocity[2]);

  myfile1 <<"Position = [" <<basePos_ori.x() <<", " <<basePos_ori.y() << ", "<<basePos_ori.z() <<"], ";
  myfile1 <<"Orientation = [" <<baseOrn_ori.x() <<", " <<baseOrn_ori.y() << ", "<<baseOrn_ori.z() <<", "<< baseOrn_ori.w()<<"]. "<<endl;
  myfile1 <<"LinearVelocity = [" <<baseLinearVelocity[0] <<", " <<baseLinearVelocity[1] << ", "<<baseLinearVelocity[2] <<"], ";
  myfile1 <<"AngularVelocity = [" <<baseAngularVelocity[0] <<", " <<baseAngularVelocity[1] << ", "<<baseAngularVelocity[2] <<"]. "<<endl;
  myfile2 <<basePos_ori.x() <<", " <<basePos_ori.y() << ", "<<basePos_ori.z() <<endl;
  myfile3 <<baseOrn_ori.x() <<", " <<baseOrn_ori.y() << ", "<<baseOrn_ori.z() <<", "<<baseOrn_ori.w()<<endl;
  myfile4 <<baseLinearVelocity[0] <<", " <<baseLinearVelocity[1] << ", "<<baseLinearVelocity[2]<<endl;
  myfile5 <<0.0 <<endl;

	std::clock_t t_start;
  std::clock_t t_start_cons;
	double span;
  double t_record;
	t_start = std::clock();
  t_start_cons = std::clock();


  double volume;
  double d_water;

  volume = 5.746E-05;
  d_water = 1200;

  btVector3 Pos;
  btQuaternion Orn;
  btVector3 Fb = btVector3(0, 0, 0);
  btVector3 C_Fb = btVector3(-30.0, 0, 0);
  btVector3 C_Fd = btVector3(0, -30, 0);
  // btVector3 C_Fb = btVector3(-2.069, -1.922, 2.0);
  btVector3 F_left = btVector3(-2000, 0, 0);
  // btVector3 F_left = btVector3(0, 0, 0);
  btVector3 C_Fz = btVector3(0, 0, 0);
  btVector3 F_z = btVector3(0, 0, -3);
  btVector3 F_down = btVector3(0, -2000, 0);
  btVector3 F_im = btVector3(0, 0, -10000);
  
  b3RobotSimulatorChangeDynamicsArgs D_args;
  // D_args.m_linearDamping = 1;
  // D_args.m_angularDamping = 1;
  
  std::clock_t tc;
  
  int *new_block = new int;
  b3RobotSimulatorLoadUrdfFileArgs new_block_arg;
  
  std::vector<bool> status(2);
  status[0] = false;
  status[1] = false;
  
  std::vector<int*> new_blocks(3);
  new_blocks[0] = new int;
  new_blocks[1] = new int;
  new_blocks[2] = new int;

  int count = 0;
  

	while (sim->canSubmitCommand())
	{



    // sim ->applyExternalForce(*block, -1, Fb, C_Fb, 1);

    // cout<<volume*d_water*9.8<<endl;

    tc = std::clock();
    if (tc > 1000000 && tc < 1500000 && status[0] == false) {
      new_block_arg.m_startPosition.setValue(30,0,1);
      new_block_arg.m_startOrientation.setEulerZYX(0,0,0);
      
      // *new_block = sim->loadURDF(block_file_other, new_block_arg);
      *new_blocks[0] = sim->loadURDF(block_file_other, new_block_arg);

      status[0] = true;
    } 
    else if (tc > 2000000 && status[1] == false) {
      new_block_arg.m_startPosition.setValue(0,30,1);
      new_block_arg.m_startOrientation.setEulerZYX(0,0,0);
      
      // *new_block = sim->loadURDF(block_file_other, new_block_arg);
      *new_blocks[1] = sim->loadURDF(block_file_other, new_block_arg);
      
      // sim->changeDynamics(*new_blocks[1], -1, mass=10000);
      D_args.m_mass = 1000000;
      sim->changeDynamics(*new_blocks[1], *new_blocks[1], D_args);
      

      new_block_arg.m_startPosition.setValue(30,30,1);
      new_block_arg.m_startOrientation.setEulerZYX(0,0,0);
      *new_blocks[2] = sim->loadURDF(block_file_other, new_block_arg);
      

      // cout << "here again" << endl;
      status[1] = true;
    }
    
    // cout << "count: " << count << endl;
    if (tc > 3000000 && count == 0) {
      // sim ->applyExternalForce(*new_blocks[0], -1, F_z, C_Fz, EF_LINK_FRAME);
    	// F_left = btVector3(-20, 0, 0);
    	sim ->applyExternalForce(*new_blocks[0], -1, F_left, C_Fb, EF_LINK_FRAME);
    	sim ->applyExternalForce(*new_blocks[2], -1, F_left, C_Fb, EF_LINK_FRAME);
    	sim ->applyExternalForce(*new_blocks[1], -1, F_im, C_Fz, EF_LINK_FRAME);
    }

    if (count == 1) {
		  sim ->applyExternalForce(*new_blocks[1], -1, F_down, C_Fd, EF_LINK_FRAME);
    	sim ->applyExternalForce(*new_blocks[2], -1, F_down, C_Fd, EF_LINK_FRAME);
    
    }


		sim->stepSimulation();
    // sim->getCameraImage(640, 640, ca1, ca2);

	


		b3KeyboardEventsData keyEvents;
		sim->getKeyboardEvents(&keyEvents);
		if (keyEvents.m_numKeyboardEvents) {
			
			for (int i = 0; i < keyEvents.m_numKeyboardEvents; i++)
			{
				b3KeyboardEvent& e = keyEvents.m_keyboardEvents[i];

				if (e.m_keyCode == '0')
				{
					if (e.m_keyState & eButtonTriggered)
					{
						if (vidLogId < 0)
						{
							vidLogId = sim->startStateLogging(STATE_LOGGING_VIDEO_MP4, "video.mp4");
						}
						else
						{
							sim->stopStateLogging(vidLogId);
							vidLogId = -1;
						}
					}
				}

				if (e.m_keyCode == 'm')
				{
					if (ballLogId < 0 && e.m_keyState & eButtonTriggered)
					{
						ballLogId = sim->startStateLogging(STATE_LOGGING_MINITAUR, "simlog.bin");
					}
					if (ballLogId >= 0 && e.m_keyState & eButtonReleased)
					{
						sim->stopStateLogging(ballLogId);
						ballLogId = -1;
					}
				}

				if (e.m_keyCode == 'r' && e.m_keyState & eButtonTriggered)
				{
					rotateCamera = 1 - rotateCamera;
				}

				if (e.m_keyCode == 'a'&& e.m_keyState & eButtonTriggered )
				{
					b3RobotSimulatorLoadUrdfFileArgs args;

					args.m_startPosition.setValue(1+pos_noise,0.5+pos_noise,3);
					args.m_startOrientation.setEulerZYX(orn_noise, orn_noise, orn_noise);

					int small_ball = sim->loadURDF("ball_0.1.urdf",args);


					btVector3 Pos;
					btQuaternion Orn;
					sim->getBasePositionAndOrientation(small_ball, Pos, Orn);


					std::cout << Orn.x() << '\n';
					printf("Pos = [%f,%f,%f], Orn = [%f,%f,%f]\n", Pos.x(), Pos.y(), Pos.z(), Orn.x()
				,Orn.y(), Orn.z());


				}

				if (e.m_keyCode == 's'&& e.m_keyState & eButtonTriggered )
				{
					b3RobotSimulatorLoadUrdfFileArgs args;

					args.m_startPosition.setValue(1+pos_noise,0.5+pos_noise,3);
					args.m_startOrientation.setEulerZYX(orn_noise, orn_noise, orn_noise);
					int small_cube = sim->loadURDF("test_cube.urdf",args);
					btVector3 Pos;
					btQuaternion Orn;
					sim->getBasePositionAndOrientation(small_cube, Pos, Orn);


					std::cout << Orn.x() << '\n';
					printf("Pos = [%f,%f,%f], Orn = [%f,%f,%f]\n", Pos.x(), Pos.y(), Pos.z(), Orn.x()
				,Orn.y(), Orn.z());

				}

				if (e.m_keyCode == 'n'&& e.m_keyState & eButtonTriggered ) {
          D_args.m_mass = 100;
          sim->changeDynamics(*new_blocks[1], *new_blocks[1], D_args);
          D_args.m_mass = 1000000;
          sim->changeDynamics(*new_blocks[0], *new_blocks[0], D_args);
		        count += 1;
				}

				// if (e.m_keyCode == 'd'&& e.m_keyState & eButtonTriggered )
				// {
		  //       Fb = btVector3(18,7, 0);
				// }

    //     if (e.m_keyCode == 'z'&& e.m_keyState & eButtonTriggered )
    //     {
    //         Fb = btVector3(0,0, 0);
    //     }

    //     if (e.m_keyCode == 'x'&& e.m_keyState & eButtonTriggered )
    //     {
    //         Fb = btVector3(0,-18, 0);
    //     }

				
			}
		}



		// b3RobotSimulatorLoadUrdfFileArgs args;
		// double number = distribution(generator);
		// args.m_startPosition.setValue(number,number,3);
		// args.m_startOrientation.setEulerZYX(number, number, number);
		// int small_cube = sim->loadURDF("test_cube.urdf",args);


		if (rotateCamera) {
			static double yaw = 0;
			double distance = 1;
			yaw += 0.1;
			btVector3 basePos;
			btQuaternion baseOrn;
			sim->getBasePositionAndOrientation(plane, basePos, baseOrn);
			sim->resetDebugVisualizerCamera(distance, -20, yaw, basePos);
		}
		b3Clock::usleep(1000. * 1000. * fixedTimeStep);


    frame_x = Pos.x();
    frame_y = Pos.y();
    frame_z = Pos.z();


	}



	printf("sim->disconnect\n");

	sim->disconnect();

	printf("delete sim\n");
	delete sim;

	printf("exit\n");

}
