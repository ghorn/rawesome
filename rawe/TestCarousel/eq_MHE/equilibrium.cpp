/*
 *    This file is part of ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2009 by Boris Houska and Hans Joachim Ferreau, K.U.Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC) under
 *    supervision of Moritz Diehl. All rights reserved.
 *
 *    ACADO Toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


 /**
 *    SINGLE POWER KITE START-UP WITH PROPELLER
 *    CARTESIAN COORDINATES (ODE FORMULATION)
 *    JULY 2011 SEBASTIEN GROS, HIGHWIND, OPTEC
 *    SEBASTIEN GROS
 */


#include <acado_optimal_control.hpp>
#include <include/acado_gnuplot/gnuplot_window.hpp>



int main( int argc, char * const argv[] ){


	//======================================================================
	// PARSE THE INPUT ARGUMENTS:
	// ----------------------------------
	double IC[4];
	
	/* arguments are passed to the main function by string
	 *  there are 'int argc' arguments
	 *		the first one is the name of the program
	 *		the next ones are arguments passed in the call like
	 *			program number1 number2
	 *  in this stupid simple parser, we directly read doubles from the arguments
	 */
	 
	int i=1;
	while (i < argc) {
		// load the i-th string argument passed to the main function
		char* input = argv[i];
		// parse the string to a double value
		IC[i-1] = atof(input);
		i++;
	}
	
	cout << "Initial Conditions" << endl;
	cout << "------------" << endl;
	for (i = 0; i < argc-1; ++i) {
		cout << i+1 << ":\t" << IC[i] << endl;
	}
	//======================================================================	
	
	
	 
	 
   USING_NAMESPACE_ACADO


// DIFFERENTIAL STATES :
// -------------------------
  
	Parameter      x;      // Position
	Parameter	   y;      //  
	Parameter      z;      //  
// -------------------------      //  -------------------------------------------
  
	Parameter     dx;      //  Speed
	Parameter     dy;      //  
	Parameter     dz;      //  
	
	Parameter    e11;
	Parameter    e12;
	Parameter    e13;
	Parameter    e21;
	Parameter    e22;
	Parameter    e23;
	Parameter    e31;
	Parameter    e32;
	Parameter    e33;
	
	
	Parameter    w1;
	Parameter    w2;
	Parameter    w3;	
// -------------------------      //  -------------------------------------------
// 	Parameter      r;      //  Kite distance
// 	Parameter     dr;      //  Kite distance / dt
	double r = IC[1];

//-------------------------      //  -------------------------------------------
 
	Parameter  delta;      //  Carousel
	Parameter ddelta;      //  
	Parameter     ur;
	Parameter     up;      //  Ailerons  
	
	// Collect all the states in a vector
	const int n_XD = 22; // Number of states
	
	IntermediateState XD[n_XD];    // Vector collecting all states
	XD[0]  = x;
	XD[1]  = y;
	XD[2]  = z;
	XD[3]  = dx;
	XD[4]  = dy;
	XD[5]  = dz;
	XD[6]  = e11;
	XD[7]  = e12;
	XD[8]  = e13;
	XD[9]  = e21;
	XD[10] = e22;
	XD[11] = e23;
	XD[12] = e31;
	XD[13] = e32;
	XD[14] = e33;
	XD[15] = w1;
	XD[16] = w2;
	XD[17] = w3;
	XD[18] = delta;
	XD[19] = ddelta;
	XD[20] = ur;
	XD[21] = up;
	
// CONTROL :
// -------------------------
	Parameter             dddelta;  //  Carousel acceleration
	Parameter             ddr;
	Parameter             dur;
	Parameter             dup;      //  Ailerons  
	
	// Collect all the controls in a vector
	const int n_U = 4; // Number of controls
	
	IntermediateState U[n_U];    // Vector collecting all states
	U[0] = dddelta;
	U[1] = ddr;
	U[2] = dur;
	U[3] = dup;
	
// PARAMETERS
// -----------------------
// 	Parameter SlackCLmax;
// 	Parameter SlackLambda;

	
	// DEFINITION OF PI :
	// ------------------------
	
	double PI = 3.1415926535897932;
	
	
	  
	
// CONSTANTS :
// ------------------------
		//  PARAMETERS OF THE KITE :
		//  -----------------------------
		double mk =  0.626;      //  mass of the kite               //  [ kg    ]
		//  A =  0.2;      //  effective area                 //  [ m^2   ]


		//   PHYSICAL CONSTANTS :
		//  -----------------------------
		double   g =    9.81;      //  gravitational constant         //  [ m /s^2]
		double rho =    1.23;      //  density of the air             //  [ kg/m^3]

		//  PARAMETERS OF THE CABLE :
		//  -----------------------------
		double rhoc = 1450.00;      //  density of the cable           //  [ kg/m^3]
// 		double cc =   1.00;      //  frictional constant            //  [       ]
		double dc = 1e-3;      //  diameter                       //  [ m     ]



		double AQ      =  PI*dc*dc/4.0;


		//CAROUSEL ARM LENGTH
		double rA = 1.085; //(dixit Kurt)

		double zT = -0.02;

// 		double ZT = 0;
		//             YT = 0.005;

		//INERTIA MATRIX (Kurt's direct measurements)
		double I1 = 0.0163;
		double I31 = 0.0006;
		double I2 = 0.0078;
		double I3 = 0.0229;

		//IMU POSITION & ANGLE
// 		double XIMU = [0.0246;-0.0116;-0.0315];
// 		double alphaIMU = 0*pi/180;//4
// 		double betaIMU = 0*pi/180;
// 		double deltaIMU = 0*pi/180;

// 		double alpha0 = -0*PI/180; 

		//TAIL LENGTH
		double LT = 0.4;


		//ROLL DAMPING
		double RD = 1e-2; 
		double PD = 0*5e-3;
		double YD = 0*5e-3;
		//WIND-TUNNEL PARAMETERS

		//Lift (report p. 67)
		double CLA = 5.064;

		double CLe = -1.924;//e-5;//0.318;//

		double CL0 = 0.239;

		//Drag (report p. 70)
		double CDA = -0.195;
		double CDA2 = 4.268;
		double CDB2 = 5;//0;//
// 		double CDe = 0.044;
// 		double CDr = 0.111;
		double CD0 = 0.026;

		//Roll (report p. 72)
		double CRB = -0.062;
		double CRAB = -0.271; 
		double CRr = -5.637e-1;//e-6;//-0.244;//

		//Pitch (report p. 74)
		double CPA = 0.293;
		double CPe = -4.9766e-1;//e-6;//-0.821;//

		double CP0 = 0.03;

		//Yaw (report p. 76)
		double CYB = 0.05;
		double CYAB = 0.229;

		double SPAN = 0.96;
		double CHORD = 0.1;



// OTHER VARIABLES :
// ------------------------
	
	IntermediateState     mc;      //  mass of the cable
	IntermediateState     m ;      //  effective inertial mass
	IntermediateState  mgrav;      //  gravific mass
	//  IntermediateState     dmc;      //  first  derivative of m     with respect to t

// ORIENTATION AND FORCES :
// ------------------------
	
	IntermediateState wind               ;      //  the wind at altitude 	
	
	IntermediateState Cf              ;      //  cable drag
	IntermediateState CD              ;      //  the aerodynamic drag coefficient
	IntermediateState CL              ;      //  the aerodynamic lift coefficient
	IntermediateState CR              ;      //  the aerodynamic roll coefficient
	IntermediateState CP              ;      //  the aerodynamic pitch coefficient
	IntermediateState CY              ;      //  the aerodynamic yaw coefficient
	
	IntermediateState F           [ 3];      //  aero forces + gravity
	IntermediateState FL          [ 3];      //  the lift force
	IntermediateState FD          [ 3];      //  the drag force
// 	IntermediateState Ff          [ 3];      //  the frictional force
// 	IntermediateState Fcable          ;      //  force in the cable
	
// 	IntermediateState er          [ 3];      // X normed to 1
	IntermediateState eTe         [ 3];      //unrotated transversal vector (psi = 0)
	IntermediateState eLe         [ 3];      //unrotated lift vector (psi = 0)
	IntermediateState we          [ 3];      //  effective wind vector
	IntermediateState wE          [ 3];      //  effective wind vector
	IntermediateState wp              ;		
// 	IntermediateState wep         [ 3];		// effective wind projected in the plane orthogonal to X
	
	IntermediateState VKite           ;     //Kite (relative) speed
	IntermediateState VKite2          ;     //Squared (relative) kite speed
		
	IntermediateState Vp;
	IntermediateState VT			[3];
	IntermediateState alpha;
	IntermediateState beta;
	IntermediateState betaTail;
	IntermediateState alphaTail;
	IntermediateState T			    [3];
	IntermediateState TB		    [3];
	
	
	// TERMS ON RIGHT-HAND-SIDE
	// OF THE DIFFERENTIAL
	// EQUATIONS              :
	// ------------------------
	
	IntermediateState de11;
	IntermediateState de12;
	IntermediateState de13;
	IntermediateState de21;
	IntermediateState de22;
	IntermediateState de23;
	IntermediateState de31;
	IntermediateState de32;
	IntermediateState de33;
	
	
//                        MODEL EQUATIONS :
// ===============================================================

	// CROSS AREA OF THE CABLE :
	// ---------------------------------------------------------------
	
// 	AQ      =  PI*dc*dc/4.0                                       ;
	
	// THE EFECTIVE MASS' :
	// ---------------------------------------------------------------
	
	mc      =  0*rhoc*AQ*r  ;   // mass of the cable
	m       =  mk + mc/3.0;   // effective inertial mass
	mgrav   =  mk + mc/2.0;   // effective inertial mass
	
	// -----------------------------   // ----------------------------
	//   dm      =  (rhoc*AQ/ 3.0)*dr;   // time derivative of the mass
	
	
	// WIND SHEAR MODEL :
	// ---------------------------------------------------------------
	
	wind       =  0.;
	
	
	// EFFECTIVE WIND IN THE KITE`S SYSTEM :
	// ---------------------------------------------------------------
	
	we[0]   = -wind + dx - ddelta*y;
	we[1]   =		  dy + ddelta*rA + ddelta*x;
	we[2]   =		  dz;
	
	VKite2 = (we[0]*we[0] + we[1]*we[1] + we[2]*we[2]); 
	VKite = sqrt(VKite2); 
	
	// CALCULATION OF THE FORCES :
	// ---------------------------------------------------------------
	
// 	// er
//     er[0] = x/r;
// 	er[1] = y/r;
// 	er[2] = z/r;
// 	
// 	//Velocity accross X (cable drag)
// 	wp = er[0]*we[0] + er[1]*we[1] + er[2]*we[2];
// 	wep[0] = we[0] - wp*er[0];
// 	wep[1] = we[1] - wp*er[1];
// 	wep[2] = we[2] - wp*er[2];
	
	//Aero coeff.
	
	
	// LIFT DIRECTION VECTOR
	// -------------------------

	//Relative wind speed in Airfoil's referential 'E'
// 	wE[0] = e11*we[0]  + e12*we[1]  + e13*we[2];
// 	wE[1] = e21*we[0]  + e22*we[1]  + e23*we[2];
// 	wE[2] = e31*we[0]  + e32*we[1]  + e33*we[2];
	wE[0] = e11*we[0]  + e21*we[1]  + e31*we[2];
	wE[1] = e12*we[0]  + e22*we[1]  + e32*we[2];
	wE[2] = e13*we[0]  + e23*we[1]  + e33*we[2];


	//Airfoil's transversal axis in fixed referential 'e'
// 	eTe[0] = e21;
// 	eTe[1] = e22;
// 	eTe[2] = e23;
	eTe[0] = e12;
	eTe[1] = e22;
	eTe[2] = e32;


	// Lift axis ** Normed to we !! **
	eLe[0] = - eTe[1]*we[2] + eTe[2]*we[1];
	eLe[1] = - eTe[2]*we[0] + eTe[0]*we[2];
	eLe[2] = - eTe[0]*we[1] + eTe[1]*we[0];

	// AERODYNAMIC COEEFICIENTS
	// ----------------------------------
	//VT = cross([w1;w2;w3],[-LT;0;0]) + wE;

	VT[0] =          wE[0];
	VT[1] = -LT*w3 + wE[1];
	VT[2] =  LT*w2 + wE[2];

	alpha = -wE[2]/wE[0];

	//Note: beta & alphaTail are compensated for the tail motion induced by omega
// 	beta = VT[1]/sqrt(VT[0]*VT[0] + VT[2]*VT[2]);
	betaTail = VT[1]/VT[0];
//     betaTail = wE(2)/sqrt(wE(1)*wE(1) + wE(3)*wE(3));
    beta = wE[1]/wE[0];
	alphaTail = -VT[2]/VT[0];

// 	CL = CLA*alpha + CLe*up     + CLr*ur + CL0;
// 	CD = CDA*alpha + CDA2*alpha*alpha + CDB2*beta*beta + CDe*up + CDr*ur + CD0;
// 	CR = -RD*w1 + CRB*betaTail + CRAB*alphaTail*betaTail + CRr*ur;
// 	CP = CPA*alphaTail + CPe*up + CPr*ur + CP0;
// 	CY = CYB*betaTail + CYAB*alphaTail*betaTail;
	CL = CLA*alpha + CLe*up + CL0;
	CD = CDA*alpha + CDA2*alpha*alpha + CDB2*beta*beta + CD0;
	CR = -RD*w1 + CRB*betaTail + CRr*ur + CRAB*alphaTail*betaTail;
	CP = -PD*w2 + CPA*alphaTail + CPe*up + CP0;
	CY = -YD*w3 + CYB*betaTail + CYAB*alphaTail*betaTail;
	
	
	
// 	Cf = rho*dc*r*VKite/8.0;

	// THE FRICTION OF THE CABLE :
	// ---------------------------------------------------------------

// 	Ff[0] = -rho*dc*r*VKite*cc*wep[0]/8.0;
// 	Ff[1] = -rho*dc*r*VKite*cc*wep[1]/8.0;
// 	Ff[2] = -rho*dc*r*VKite*cc*wep[2]/8.0;

	// LIFT :
	// ---------------------------------------------------------------
	CL = 0.2*CL;
	FL[0] =  rho*CL*eLe[0]*VKite/2.0;
	FL[1] =  rho*CL*eLe[1]*VKite/2.0;
	FL[2] =  rho*CL*eLe[2]*VKite/2.0;

	// DRAG :
	// ---------------------------------------------------------------

	FD[0] = -rho*VKite*CD*we[0]/2.0;
	FD[1] = -rho*VKite*CD*we[1]/2.0; 
	FD[2] = -rho*VKite*CD*we[2]/2.0; 


	// FORCES (AERO)
	// ---------------------------------------------------------------

	F[0] = FL[0] + FD[0];
	F[1] = FL[1] + FD[1];
	F[2] = FL[2] + FD[2];


	// ATTITUDE DYNAMICS
	// -----------------------------------------------------------

// 	de11 = e21*w3 - e31*w2 + ddelta*e12; 
// 	de21 = e31*w1 - e11*w3 + ddelta*e22; 
// 	de31 = e11*w2 - e21*w1 + ddelta*e32; 
// 	de12 = e22*w3 - e32*w2 - ddelta*e11; 
// 	de22 = e32*w1 - e12*w3 - ddelta*e21; 
// 	de32 = e12*w2 - e22*w1 - ddelta*e31; 
// 	de13 = e23*w3 - e33*w2; 
// 	de23 = e33*w1 - e13*w3; 
// 	de33 = e13*w2 - e23*w1;
    de11 = e12*w3 - e13*w2 + ddelta*e21;
    de12 = e13*w1 - e11*w3 + ddelta*e22;
    de13 = e11*w2 - e12*w1 + ddelta*e23;
    de21 = e22*w3 - e23*w2 - ddelta*e11;
    de22 = e23*w1 - e21*w3 - ddelta*e12; 
    de23 = e21*w2 - e22*w1 - ddelta*e13;
    de31 = e32*w3 - e33*w2;
    de32 = e33*w1 - e31*w3;
    de33 = e31*w2 - e32*w1; 



	//////////////////////////////////////////////////////////////////////// 
	//                                                                    // 
	//  AUTO-GENERATED EQUATIONS (S. Gros, HIGHWIND, OPTEC, KU Leuven)    // 
	//                                                                    // 
	//////////////////////////////////////////////////////////////////////// 

	// Inverse of the Generalized inertia matrix 
	IntermediateState IMA(3,3); 
	IMA(0,0) = 1/m; 
	IMA(0,1) = 0; 
	IMA(0,2) = 0; 
	IMA(1,0) = 0; 
	IMA(1,1) = 1/m; 
	IMA(1,2) = 0; 
	IMA(2,0) = 0; 
	IMA(2,1) = 0; 
	IMA(2,2) = 1/m; 

	// G1 (right hand-side up) 
	IntermediateState G1(3,1); 
// 	G1(0,0) = F[0] + ddelta*dy*m + dddelta*m*y; 
// 	G1(1,0) = F[1] - (dddelta*m*(2*rA + 2*x))/2 - ddelta*dx*m; 
// 	G1(2,0) = F[2] - g*mgrav; 
	G1(0,0) = F[0] + (m*(2*ddelta*dy + 2*ddelta*ddelta*rA + 2*ddelta*ddelta*x))/2 + ddelta*dy*m + dddelta*m*y; 
	G1(1,0) = F[1] - (m*(2*ddelta*dx - 2*ddelta*ddelta*y))/2 - (dddelta*m*(2*rA + 2*x))/2 - ddelta*dx*m; 
	G1(2,0) = F[2] - g*mgrav; 

	// G2 (right hand-side down) 
	IntermediateState G2; 
	G2 = - dx*dx - dy*dy - dz*dz; 

	// NabG 
	IntermediateState NabG(3,1); 
	NabG(0,0) = x; 
	NabG(1,0) = y; 
	NabG(2,0) = z; 

	// NabGT 
	IntermediateState NabGT(1,3); 
	NabGT(0,0) = x; 
	NabGT(0,1) = y; 
	NabGT(0,2) = z; 

	// LambdaFac * lambda = lambdaright 
	IntermediateState LambdaFac; 
	LambdaFac = NabGT*IMA*NabG; 

	// lambdaright 
	IntermediateState lambdaright; 
	lambdaright = NabGT*IMA*G1 - G2; 

	// lambda 
	IntermediateState lambda; 
	lambda = lambdaright/LambdaFac; 

	// ddq (accelerations) 
	IntermediateState ddq(3,1); 
	ddq = IMA*(G1-NabG*lambda); 

	// Consistency Conditions 
	IntermediateState Const, dConst; 
	Const = - r*r/2 + x*x/2 + y*y/2 + z*z/2; 
	dConst = dx*x + dy*y + dz*z; 

	///////////////////////////// END OF AUTO-GENERATED CODE ////////////////////////////////////////////////////// 
	
	// TORQUES (BRIDLE)
	// ---------------------------------------------------------------

	TB[0] =  zT*(e12*lambda*x + e22*lambda*y + e32*lambda*z);
	TB[1] = -zT*(e11*lambda*x + e21*lambda*y + e31*lambda*z);
	TB[2] =  0;
	
	// TORQUES (AERO)
	// ---------------------------------------------------------------

	T[0] =  0.5*rho*VKite2*SPAN*CR   + TB[0];
	T[1] =  0.5*rho*VKite2*CHORD*CP  + TB[1];
	T[2] =  0.5*rho*VKite2*SPAN*CY   + TB[2];


	IntermediateState dw1, dw2, dw3;
// 	dw1 = (I31*(T[2] + w2*(I1*w1 + I31*w3) - I2*w1*w2))/(I31*I31 - I1*I3) - (I3*(T[0] - w2*(I31*w1 + I3*w3) + I2*w2*w3))/(I31*I31 - I1*I3); 
// 	dw2 = (T[1] + w1*(I31*w1 + I3*w3) - w3*(I1*w1 + I31*w3))/I2; 
// 	dw3 = (I31*(T[0] - w2*(I31*w1 + I3*w3) + I2*w2*w3))/(I31*I31 - I1*I3) - (I1*(T[2] + w2*(I1*w1 + I31*w3) - I2*w1*w2))/(I31*I31 - I1*I3); 
	dw1 = (I31*T[2])/(I31*I31 - I1*I3) - (I3*T[0])/(I31*I31 - I1*I3); 
	dw2 = T[1]/I2; 
	dw3 = (I31*T[0])/(I31*I31 - I1*I3) - (I1*T[2])/(I31*I31 - I1*I3); 
	
	
// 	Fcable = lambda*r;
	
// 	IntermediateState ddx, ddy, ddz;
// 	ddx = ddq(0,0);
// 	ddy = ddq(1,0);
// 	ddz = ddq(2,0);
	
	
	// ===============================================================
	//                        END OF MODEL EQUATIONS
	// ===============================================================
	
	
	
	IntermediateState Cost = 0;
	for ( i=0; i < n_U; i++ ) {
		Cost += U[i]*U[i];
	}

    // DEFINE THE NLP:
    // ----------------------------------
	NLP nlp;
	nlp.minimize( Cost );
// 	nlp.minimize( ur*ur + up*up );
	nlp.minimize( dz*dz );
	
    
    // CONSTRAINTS:
    // ---------------------------------
	
	// State bounds
	nlp.subjectTo(  -1 <= up  <= 1 );
	nlp.subjectTo(  -3.2767 <= ur  <= 3.2767 );
// 	nlp.subjectTo(  0.0 <= x  );
// 	nlp.subjectTo(   rA <= y  );
	
/*	nlp.subjectTo( -1.0 <= e11 <= 1.0  );
	nlp.subjectTo( -1.0 <= e12 <= 1.0  );
	nlp.subjectTo( -1.0 <= e13 <= 1.0  );
	nlp.subjectTo( -1.0 <= e21 <= 1.0  );
	nlp.subjectTo( -1.0 <= e22 <= 1.0  );
	nlp.subjectTo( -1.0 <= e23 <= 1.0  );
	nlp.subjectTo( -1.0 <= e31 <= 1.0  );
	nlp.subjectTo( -1.0 <= e32 <= 1.0  );
	nlp.subjectTo( -1.0 <= e33 <= 1.0  );*/
	
	nlp.subjectTo(  delta == IC[2] );
	nlp.subjectTo( ddelta == IC[3] );
	
	// Flight envelope
	nlp.subjectTo( 0 <= CL <= 1  );
// 	nlp.subjectTo( CL == 0.5  );
	nlp.subjectTo( 1 <= lambda/10  );
// 	nlp.subjectTo( -10*PI/180 <= beta <= 10*PI/180 );
// 	nlp.subjectTo( -10*PI/180 <= alpha <= 10*PI/180 );
// 	nlp.subjectTo( -15*PI/180 <= alphaTail <= 15*PI/180 );
	nlp.subjectTo( 0 <= wE[0]  );
	
	// Invariants
	nlp.subjectTo( Const == 0 );
	nlp.subjectTo( dConst == 0 );
	
	nlp.subjectTo( e11*e11 + e12*e12 + e13*e13 - 1 == 0 );
	nlp.subjectTo( e11*e21 + e12*e22 + e13*e23     == 0 );
	nlp.subjectTo( e11*e31 + e12*e32 + e13*e33     == 0 );
	nlp.subjectTo( e21*e21 + e22*e22 + e23*e23 - 1 == 0 );
	nlp.subjectTo( e21*e31 + e22*e32 + e23*e33     == 0 );
	nlp.subjectTo( e31*e31 + e32*e32 + e33*e33 - 1 == 0 );
	
	// Equilibrium
// ROTATE BACK THE SYSTEM DYNAMICS:
// ---------------------------------------------------------------

// 	nlp.subjectTo(      z == IC[0] );

// 	nlp.subjectTo( dx == 0 );
	nlp.subjectTo( dy == 0 );
// 	nlp.subjectTo( dz == 0 ); // needs to be fixed

// 	nlp.subjectTo( 1*ddq(0,0) == 0 );
	nlp.subjectTo( 1*ddq(1,0) == 0 );
	nlp.subjectTo( 1*ddq(2,0) == 0 );

	nlp.subjectTo( w1 - ddelta*e31 == 0 );
	nlp.subjectTo( w2 - ddelta*e32 == 0 );
	nlp.subjectTo( w3 - ddelta*e33 == 0 );
	
	nlp.subjectTo( dw1 == 0 );
	nlp.subjectTo( dw2 == 0 );
	nlp.subjectTo( dw3 == 0 );
	
	
	nlp.subjectTo( dddelta == 0 );
// 	nlp.subjectTo( ddr == 0 );
	nlp.subjectTo( dur == 0 );
	nlp.subjectTo( dup == 0 );
	
	nlp.subjectTo( ur == 0 );
	nlp.subjectTo( up == 0 );
	
	nlp.subjectTo( w2 >= 0 );
	
	// DEFINE AN OPTIMIZATION ALGORITHM AND SOLVE THE OCP:
    // ---------------------------------------------------
	OptimizationAlgorithm algorithm(nlp);	

	algorithm.set                         ( MAX_NUM_ITERATIONS, 100 );
    algorithm.set                         ( KKT_TOLERANCE    , 1e-6 );
	

	algorithm.initializeParameters("EQ_init.txt"      );
	
	algorithm.set( MIN_LINESEARCH_PARAMETER, 1e-4 );
	algorithm.set( HESSIAN_APPROXIMATION, EXACT_HESSIAN );
	algorithm.set( KKT_TOLERANCE, 1e-10 );
	
	// 	algorithm.set( PRINT_SCP_METHOD_PROFILE, YES );
	
	algorithm.solve();
	
	algorithm.getParameters("EQ_params.txt"    );
	

    return 0;
}








// 	IntermediateState EC1 = (-x)*ddelta       + dy;
// 	IntermediateState EC2 = (-y)*ddelta       - dx;
// 	IntermediateState EC3 = dz;
// 	IntermediateState EC4 = ( ( -dddelta )*x + ( -ddelta*ddelta )*y ) + 2*(-dx)*ddelta  + ddX(1,0);
// 	IntermediateState EC5 = ( (  ddelta*ddelta )*x + ( -dddelta )*y ) + 2*(-dy)*ddelta  - ddX(0,0);
// 	IntermediateState EC6 = 1*ddX(2,0);
// 	IntermediateState EC7 = dr;
// 	IntermediateState EC8 = z - IC[0];
// 	IntermediateState EC9 = r - IC[1];
	
// 	Parameter Z11;
// 	Parameter Z11, Z12, Z13, Z14, Z15, Z16, Z17, Z18, Z19;
// 	Parameter Z21, Z22, Z23, Z24, Z25, Z26, Z27, Z28, Z29;
// 	Parameter Z31, Z32, Z33, Z34, Z35, Z36, Z37, Z38, Z39;
// 	Parameter Z41, Z42, Z43, Z44, Z45, Z46, Z47, Z48, Z49;
// 	Parameter Z51, Z52, Z53, Z54, Z55, Z56, Z57, Z58, Z59;
// 	Parameter Z61, Z62, Z63, Z64, Z65, Z66, Z67, Z68, Z69;
// 	Parameter Z71, Z72, Z73, Z74, Z75, Z76, Z77, Z78, Z79;
// 	IntermediateState Z(7,9);
// 	Z(0,0) = Z11; Z(0,1) = Z12; Z(0,2) = Z13; Z(0,3) = Z14; Z(0,4) = Z15; Z(0,5) = Z16; Z(0,6) = Z17; Z(0,7) = Z18; Z(0,8) = Z19;
// 	Z(1,0) = Z21; Z(1,1) = Z22; Z(1,2) = Z23; Z(1,3) = Z24; Z(1,4) = Z25; Z(1,5) = Z26; Z(1,6) = Z27; Z(1,7) = Z28; Z(1,8) = Z29;
// 	Z(2,0) = Z31; Z(2,1) = Z32; Z(2,2) = Z33; Z(2,3) = Z34; Z(2,4) = Z35; Z(2,5) = Z36; Z(2,6) = Z37; Z(2,7) = Z38; Z(2,8) = Z39;
// 	Z(3,0) = Z41; Z(3,1) = Z42; Z(3,2) = Z43; Z(3,3) = Z44; Z(3,4) = Z45; Z(3,5) = Z46; Z(3,6) = Z47; Z(3,7) = Z48; Z(3,8) = Z49;
// 	Z(4,0) = Z51; Z(4,1) = Z52; Z(4,2) = Z53; Z(4,3) = Z54; Z(4,4) = Z55; Z(4,5) = Z56; Z(4,6) = Z57; Z(4,7) = Z58; Z(4,8) = Z59;
// 	Z(5,0) = Z61; Z(5,1) = Z62; Z(5,2) = Z63; Z(5,3) = Z64; Z(5,4) = Z65; Z(5,5) = Z66; Z(5,6) = Z67; Z(5,7) = Z68; Z(5,8) = Z69;
// 	Z(6,0) = Z71; Z(6,1) = Z72; Z(6,2) = Z73; Z(6,3) = Z74; Z(6,4) = Z75; Z(6,5) = Z76; Z(6,6) = Z77; Z(6,7) = Z78; Z(6,8) = Z79;
// 	
// 	
// 	IntermediateState JEC(2,9);
// 	JEC(0,0) =  x-xA;
// 	JEC(0,1) =  y-yA;
// 	JEC(0,2) =  z;
// 	JEC(0,3) =  0;
// 	JEC(0,4) =  0;
// 	JEC(0,5) =  0;
// 	JEC(0,6) = -0.5;
// 	JEC(0,7) =  0;
// 	JEC(1,0) =  dx-dxA;
// 	JEC(1,1) =  dy-dyA;
// 	JEC(1,2) =  dz;
// 	JEC(1,3) =  x-xA;
// 	JEC(1,4) =  y-yA;
// 	JEC(1,5) =  z;
// 	JEC(1,6) = -dr;
// 	JEC(1,7) = -r;
// 	
// 	IntermediateState NS(2,9);
// 	
// 	int k;
// 	int j;
// 	for ( i=0; i<2; i++ ){
// 		for ( k=0; k<9; k++ ){
// 			NS(i,k) = 0;
// 			for ( j=0; j<7; j++ ){
// // 				cout << i << k << j << endl;
// 				NS(i,k) +=  JEC(i,j)*Z(j,k);
// 			}
// 		}
// 	}

