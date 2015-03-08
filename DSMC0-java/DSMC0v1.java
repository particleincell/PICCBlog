/*
 * Bird's DSMC0.FOR translated to JAVA

* PROGRAM DSMC0
*
//test of collision procedures in a uniform gas
*
//SI units are used throughout
*/
package demos;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Random;

public class DSMC0v1 
{
    public static void main(String[] args) 
    {
        new DSMC().run();
    }
}
class DSMC
{
	final int MNM=1000;		//MNM  is the maximum number of molecules
	final int MNC=50;		//MNC  is the maximum number of sub 
	final int MNSC=400;		//MNSC is the maximum number of sub-cells 
	final int MNSP=5;		//MNSP is the maximum number of molecular species
	final int MNSG=1;		//MNSG is the number of species groups for collision sampling

	double COL[][] = new double[MNSP+1][MNSP+1];	//COL(M,N) is the number of collisions between species N-M molecules
	double MOVT;    //MOVT the total number of molecular moves
	long NCOL;    //NCOL is the total number of collisions
	double SELT;	//SELT the total number of pair selections
	double SEPT;	//SEPT the sum of collision pair separations
	double CS[][][] = new double[7+1][MNC+1][MNSP+1];	//CS(N,M,L) sampled information on species L in cell M
												//---- N=1 number sum
												//----N=2,3,4 sum of u,v,w
												//----N=5,6,7 sum of u*u,v*v,w*w
	/*MOLS block*/
	int NM;                          //NM is the number of molecules
	double PP[] = new double[MNM+1];      //PP(M)  is the x coordinate molecule M
	double PV[][] = new double[3+1][MNM+1]; //PV(1 to 3,M)  u,v,w velocity components of molecule M
	int IPL[] = new int[MNM+1];     //IPL(M) sub-cell number for molecule M
	int IPS[] = new int[MNM+1];     //IPS(M) species code number
	int IR[] = new int[MNM+1];      //IR(M)  cross-reference array (molecule numbers in order of sub-cells)

	/*cells block*/
	double CC[] = new double[MNC+1];      //CC(M) is the cell volume
	double CG[][] = new double[3+1][MNC+1];   //--CG(N,M) is the geometry related information on cell M
										  //----N=1 the minimum x coordinate
										  //----N=2 the maximum x coordinate
										  //----N=3 the cell width
	int IC[][][] = new int[2+1][MNC+1][MNSG+1]; //--IC(N,M,L) information on the molecules of species group L in cell M
												//----N=1 (start address -1) of the molecule numbers in the array IR
												//----N=2 the number of molecules in the cell
	int ISC[] = new int[MNSC+1];    //ISC(M) the cell in which the sub-cell lies
	double CCG[][][][] = new double [2+1][MNC+1][MNSG+1][MNSG+1]; //--CCG(N,M,L,K) is for collisions between species groups L-K in cell M
														  //----N=1 is the maximum value of (relative speed)*(coll. cross-section)
														  //----N=2 is the remainder when the selection number is rounded
	int ISCG[][][] = new int[2+1][MNSC+1][MNSG+1]; //--ISCG(N,M,L) is the information on species group L in sub-cell M
												   //----N=1 (start address -1) of the molecule numbers in the array IR
												   //----N=2 the number of molecules in the sub-cell
	int IG[][] = new int[2+1][MNSG+1];    //--IG(2,M) information on group L molecules
											//----N=1 (start address -1) of the molecule numbers in the array IR
											//----N=2 the number of molecules in the cell

	/*gas block*/
	double SP[][] = new double[5+1][MNSP+1];   //--SP(N,M) information on species M
										   //----N=1 the reference cross-section (diameter in the data)
										   //----N=2 the reference temperature
										   //----N=3 the viscosity-temperature power law
										   //----N=4 the reciprocal of the VSS scattering parameter
										   //----N=5 the molecular mass
	double SPM[][][] = new double[6+1][MNSP+1][MNSP+1]; ////SPM(N,M,L) information on the interaction between L-M molecules
												////--N=1  the reference cross-section (diameter in the data)
												////--N=2  the reference temperature
												////--N=3  the viscosity-temperature power law
												////--N=4  the reciprocal of the VSS scattering parameter
												////--N=5  the reduced mass
												////--N=6  the Gamma function of (5/2 - viscosity-temperature power law)
	int ISP[] = new int[MNSP+1];    //ISP(M) the colision sampling group in which species M lies


	/*SAMP block*/
	double TIME;            //TIME time
	int NPR;         ////NPR the number of output/restart file update cycles
	int NSMP;        ////NSMP the total number of samples
	double ND;         ////FND the stream number density
	double TMP;        ////FTMP the stream temperature
	double FSP[] = new double[MNSP+1];    ////FSP(M) the fraction of species M in the stream
	int ISPD;		////ISPD relates to the setting of data for colls. between unlike mols.
					////--set to 0 if data is set automatically to the mean values
					////--set to 1 if the values are set explicitly in the data

	/*comp block*/
	double FNUM;    ////FNUM  is the number of real molecules represented by a simulated mol.
	double DTM;     ////DTM is the time step
	int NIS;     ////NIS is the number of time steps between samples
	int NSP;     ////NSP is the number of samples between restart and output file updates
	int NPS;     ////NPS is the estimated number of samples to steady flow
	int NPT;     ////NPT is the number of file updates to STOP

	/*geom block*/
	double CW;      ////CW is the cell width
	int NSC;     ////NSC is the number of sub-cells per cell
	double XF;      ////XF is the minimum x coordinate
	double XR;      ////XR is the maximum x coordinate

	/*const block*/
	double PI;      ////PI is pi and SPI is the square root of pi
	double SPI;
	double BOLTZ;   ////BOLTZ is the Boltzmann constant
	
	/*elast block*/  
	double VRC[]=new double[3+1];
	double VRR;
	double VR;
	int LS,MS;
	int L,M;
	int MSC;
	double CVR;
	//int K,L,M,N;
	int MM,MN,NN;
	
	public void run()
	{
		//System.out.printf("RF %.10g %.10g %.10g %.10g %.10g\n",RF(0),RF(0),RF(0),RF(0),RF(0));
		System.out.println(" INPUT 0,1 FOR CONTINUING,NEW CALCULATION:- ");
		int NQL= 1; /*instead of READ*/

		if (NQL == 1)
		{
			INIT0();
		}
		else
		{
	   /*     WRITE (*,*) ' READ THE RESTART FILE'
			OPEN (4,FILE='DSMC0.RES',STATUS='OLD',FORM='UNFORMATTED')
			READ (4) BOLTZ,CC,CCG,CG,COL,CS,CW,DTM,FNUM,FTMP,IC,IPL,IPS,IR,
			 &           ISC,ISCG,ISP,MOVT,NCOL,NIS,NM,NSC,NSMP,NPR,NPT,NSP,PI,
			 &           PP,PV,SELT,SEPT,SP,SPI,SPM,TIME,XF,XR
			CLOSE (4)
 */       }

		if (NQL == 1) SAMPI0();

		while (NPR<NPT)
		{
			NPR=NPR+1;
			for (int JJJ=1;JJJ<=NSP;JJJ++)
			{
				for (int III=1;III<=NIS;III++)
				{
					TIME=TIME+DTM;
	
					System.out.printf(" DSMC0:- Move %d %d   of  %d %d\t%d Collisions\n",III,JJJ,NIS,NSP,NCOL);
			
					MOVE0();
					INDEXM();
					COLLM();
					}

				SAMPLE0();
			}	/*for*/

			/*WRITE (*,*) ' WRITING RESTART AND OUTPUT FILES',NPR,'  OF ',NPT
			OPEN (4,FILE='DSMC0.RES',FORM='UNFORMATTED')
			WRITE (4) BOLTZ,CC,CCG,CG,COL,CS,CW,DTM,FNUM,FTMP,IC,IPL,IPS,IR,
			&          ISC,ISCG,ISP,MOVT,NCOL,NIS,NM,NSC,NSMP,NPR,NPT,NSP,PI,
			&          PP,PV,SELT,SEPT,SP,SPI,SPM,TIME,XF,XR
			CLOSE (4)
			**/
			System.out.println(NCOL+" Collisions");
	
			OUT0();
		} /*while*/ 
}

	/**
	* INIT0
	*/
	public void INIT0()
	{
		PI=3.141592654;
		SPI=Math.sqrt(PI);
		BOLTZ=1.3806E-23;

		DATA0();

		//--set information on the cross-species collisions
		if (MNSP==1) ISPD=0;
		for (int N=1;N<=MNSP;N++)
		{
			for (M=1;M<=MNSP;M++)
			{
				if ((ISPD==0)||(N==M)) 
				{
					//--the collision cross section is assumed to be given by eqn (1.35)
					SPM[1][N][M]=0.25*PI*(SP[1][N]+SP[1][M])*(SP[1][N]+SP[1][M]);
				
					//--mean values are used for ISPD=0
					SPM[2][N][M]=0.5*(SP[2][N]+SP[2][M]);
					SPM[3][N][M]=0.5*(SP[3][N]+SP[3][M]);
					SPM[4][N][M]=0.5*(SP[4][N]+SP[4][M]);
				}
				else
				{
					//--the cross-collision diameter is converted to the cross-section
					SPM[1][N][M]=PI*SPM[1][N][M]*SPM[1][N][M];
				}	/*end if*/
				
				//*--the reduced mass is defined in eqn (2.7)
				SPM[5][N][M]=(SP[5][N]/(SP[5][N]+SP[5][M]))*SP[5][M];
				SPM[6][N][M]=GAM(2.5-SPM[3][N][M]);
			}	/*for*/
		}	/*for*/

		//*--initialise variables
		TIME=0.;
		NM=0;
	  
		CG[1][1]=XF;
		CW=(XR-XF)/MNC;
		for (M=1;M<=MNC;M++)
		{
			if (M>1) CG[1][M]=CG[2][M-1];
			CG[2][M]=CG[1][M]+CW;
			CG[3][M]=CW;
			CC[M]=CW;
			for (L=1;L<=MNSG;L++)
				for (int K=1;K<=MNSG;K++)
				{
					//*--the maximum value of the (rel. speed)*(cross-section) is set to a
					//*--reasonable, but low, initial value and will be increased as necessary
					CCG[2][M][L][K]=RF(0);
					CCG[1][M][L][K]=SPM[1][1][1]*300.*Math.sqrt(TMP/300.);
				}
		}	/*for*/

		//--set sub-cells
		for (int N=1;N<=MNC;N++)
			for (M=1;M<=NSC;M++)
			{
				L=(N-1)*NSC+M;
				ISC[L]=N;
			}

		//--generate initial gas in equilibrium at temperature FTMP
		for (L=1;L<=MNSP;L++)
		{
			double REM=0;
			//--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
			double VMP=Math.sqrt(2.*BOLTZ*TMP/SP[5][L]);	
			for (int N=1;N<=MNC;N++)
			{
				//--A is the number of simulated molecules of species L in cell N to
				//--simulate the required concentrations at a total number density of FND
				double A=ND*CG[3][N]*FSP[L]/FNUM+REM;
				if (N<MNC) 
				{
					MM=(int)A;		//! added (int)
					REM=(A-MM);		//--the remainder REM is carried forward to the next cell
				} else
				{
					MM=(int)(A+0.5);	/*this was NINT, round to nearest whole number*/
				}
				
				for (M=1;M<=MM;M++)
				{
					if (NM<MNM) /*this was <=*/
					{
						//--round-off error could have taken NM to MNM+1
						NM=NM+1;
						IPS[NM]=L;
						PP[NM]=CG[1][N]+RF(0)*(CG[2][N]-CG[1][N]);
						/*!added (int)*/
						IPL[NM]=(int)((PP[NM]-CG[1][N])*(NSC-0.001)/CG[3][N]+1+NSC*(N-1));
						//--species, position, and sub-cell number have been set
					  
						//--velocity components have been set
						for (int K=1;K<=3;K++)
							PV[K][NM] = RVELC(A,VMP);
					}	/*end if*/
				} /*for*/
			}	/*for*/
		}
		
		System.out.printf ("%d MOLECULES\n",NM);
	}

/*   SAMPI0.FOR
*
**/
	public void SAMPI0()
	{
		NPR=0;
		NCOL=0;
		NSMP=0;
		MOVT=0.;
		SELT=0.;
		SEPT=0.;
		
		for (L=1;L<=MNSP;L++)
			for (int N=1;N<=MNC;N++)
			{
				CS[1][N][L]=1.E-6;
				for (M=2;M<=7;M++)
					CS[M][N][L]=0.0;
			}
	
		for (M=1;M<=MNSP;M++)
			for (int N=1;N<=MNSP;N++)
				COL[M][N]=0.;
	}
	
	public void MOVE0()
	{
		for (int N=1;N<=NM;N++)
		{
			MOVT=MOVT+1;
			MSC=IPL[N];
			int MC=ISC[MSC];			//MC is the initial cell number
			double XI=PP[N];
			double DX=PV[1][N]*DTM;
			double X=XI+DX;
			//molecule N at XI is moved by DX to X
			if (X<XF)
			{
				//specular reflection from the minimum x boundary at x=XF (eqn (11.7))
				X=2.*XF-X;
				PV[1][N]=-PV[1][N];
			}
			
			if (X>XR)
			{
				//specular reflection from the maximum x boundary at x=XR (eqn (11.7))
				X=2.*XR-X;
				PV[1][N]=-PV[1][N];
			}
			
			if (X<CG[1][MC]||X>CG[2][MC])
			{
				//the molecule has moved from the initial cell
				MC = (int)((X-XF)/CW+0.99999);
				if (MC==0) MC=1;
				//MC is the new cell number (note avoidance of round-off error)
			}
			
			MSC=(int)(((X-CG[1][MC])/CG[3][MC])*(NSC-.001)+1+NSC*(MC-1));
			//MSC is the new sub-cell number
			IPL[N]=MSC;
			PP[N]=X;
		}
	}
  
	//the NM molecule numbers are arranged in order of the molecule groups
	//and, within the groups, in order of the cells and, within the cells,
	//in order of the sub-cells
    public void INDEXM()
	{
		for (MM=1;MM<=MNSG;MM++)
		{
			IG[2][MM]=0;
			for (NN=1;NN<=MNC;NN++)
				IC[2][NN][MM]=0;
      
			for (NN=1;NN<=MNSC;NN++)
				ISCG[2][NN][MM]=0;
		}
	  
		for (int N=1;N<=NM;N++)
		{
			LS=IPS[N];
			int MG=ISP[LS];
			IG[2][MG]=IG[2][MG]+1;
			MSC=IPL[N];
			ISCG[2][MSC][MG]=ISCG[2][MSC][MG]+1;
			int MC=ISC[MSC];
			IC[2][MC][MG]=IC[2][MC][MG]+1;
		}
		//--number in molecule groups in the cells and sub-cells have been counte
		
		M=0;
		for (L=1;L<=MNSG;L++)
		{
			IG[1][L]=M;	//--the (start address -1) has been set for the groups
			M=M+IG[2][L];
		}
		
		for (L=1;L<=MNSG;L++)
		{
			//--the (start address -1) has been set for the cells
			M=IG[1][L];
			for (int N=1;N<=MNC;N++)
			{
				IC[1][N][L]=M;
				M=M+IC[2][N][L];
			}
					
			//--the (start address -1) has been set for the sub-cells
			M=IG[1][L];
			for (int N=1;N<=MNSC;N++)
			{
				ISCG[1][N][L]=M;
				M=M+ISCG[2][N][L];
				ISCG[2][N][L]=0;
			}
		}
		
		for (int i=1;i<MNM+1;i++)
			IR[i]=-1;
					
		for (int N=1;N<=NM;N++)
		{
			//*--the molecule number N has been set in the cross-reference array
			LS=IPS[N];
			int MG=ISP[LS];
			MSC=IPL[N];
			ISCG[2][MSC][MG]=ISCG[2][MSC][MG]+1;
			int K=ISCG[1][MSC][MG]+ISCG[2][MSC][MG];
			IR[K]=N;
		}
	}
	
	//calculates collisions appropriate to DTM in a monatomic gas mixture
	//VRC(3) are the pre-collision components of the relative velocity
	public void COLLM()
	{
		for (int N=1;N<=MNC;N++)
		{
			//--consider collisions in cell N
			for (NN=1;NN<=MNSG;NN++)
			{
					for (MM=1;MM<=MNSG;MM++)
					{
						double SN=0.;
						for (int K=1;K<=MNSP;K++)
						  if (ISP[K]==MM) SN=SN+CS[1][N][K];
						
						//*--AVN is the average number of group MM molecules in the cell
						double AVN;
						if (SN>1.0)
						  AVN=SN/(double)(NSMP);
						else
						  AVN=IC[2][N][MM];
						//*--ASEL is the number of pairs to be selected, see eqn (11.5)
						double ASEL=0.5*IC[2][N][NN]*AVN*FNUM*CCG[1][N][NN][MM]*DTM/CC[N]+CCG[2][N][NN][MM];
						
						int NSEL=(int)ASEL;
						CCG[2][N][NN][MM]=ASEL-NSEL;
						
						if (NSEL>0)
						{
						//*--if there are insufficient molecules to calculate collisions,
						//*--the number NSEL is added to the remainer CCG(2,N,NN,MM)
							if (((NN!=MM)&&(IC[2][N][NN]<1 || IC[2][N][MM]<1))
								||((NN==MM)&&(IC[2][N][NN]<2))) 
							{
								CCG[2][N][NN][MM]=CCG[2][N][NN][MM]+NSEL;
							}
							else
							{
								double CVM=CCG[1][N][NN][MM];
								SELT=SELT+NSEL;
							
								for (int ISEL=1;ISEL<=NSEL;ISEL++)
								{
									/*note, this sets L and M as well as LS and MS*/
									SELECT(N);
									
									//*--if necessary, the maximum product in CVM is upgraded
									if (CVR>CVM) CVM=CVR;
							  
									//*--the collision is accepted with the probability of eqn (11.6)
									if (RF(0)<CVR/CCG[1][N][NN][MM]) 
									{
										NCOL=NCOL+1;
										SEPT=SEPT+Math.abs(PP[L]-PP[M]);
										COL[LS][MS]=COL[LS][MS]+1.00;
										COL[MS][LS]=COL[MS][LS]+1.00;
			
										ELASTIC();
									}/*if*/
								}	/*for*/
								CCG[1][N][NN][MM]=CVM;
							} /*if*/
						} /*if*/
				}
			}
		}	/*for*/
	} 

	public void SELECT(int N)
	{
		//selects a potential collision pair and calculates the product of the
		//collision cross-section and relative speed
		int K=(int)(RF(0)*(IC[2][N][NN]-0.0001))+IC[1][N][NN]+1;
		L=IR[K];
		M = L;
	
		//the first molecule L has been chosen at random from group NN in cell
		do {
			MSC=IPL[L];
			if ((NN==MM && ISCG[2][MSC][MM]==1) || (NN!=MM&&ISCG[2][MSC][MM]==0))
			{
				//if MSC has no type MM molecule find the nearest sub-cell with one
				int NST=1;
				int NSG=1;
				do 
				{
					int INC=NSG*NST;
					NSG=-NSG;
					NST=NST+1;
					MSC=MSC+INC;
				}	while (MSC<1 || MSC>MNSC ||
						   ISC[MSC] != N || ISCG[2][MSC][MM]<1);
			}
			
			//the second molecule M is now chosen at random from the group MM
			//molecules that are in the sub-cell MSC
			K=(int)(RF(0)*(ISCG[2][MSC][MM]-0.0001))+ISCG[1][MSC][MM]+1;
			M=IR[K];
			//choose a new second molecule if the first is again chosen
		} while (L==M);		
			
		for (K=1;K<=3;K++)
		{
			VRC[K] = PV[K][L]-PV[K][M];
			//VRC(1 to 3) are the components of the relative velocity
		}

		VRR=VRC[1]*VRC[1]+VRC[2]*VRC[2]+VRC[3]*VRC[3];
		VR=Math.sqrt(VRR);//VR is the relative speed
		LS=IPS[L];
		MS=IPS[M];
		CVR=VR*SPM[1][LS][MS]*Math.pow(2.*BOLTZ*SPM[2][LS][MS]/(SPM[5][LS][MS]*VRR),SPM[3][LS][MS]-0.5)/SPM[6][LS][MS];
		//the collision cross-section is based on eqn (4.63)
	}
  
	//generates the post-collision velocity components
	public void ELASTIC()
	{
		double VRCP[] = new double[3+1];		//VRCP(3) are the post-collision components of the relative velocity
		double VCCM[] = new double[3+1];		//VCCM(3) are the components of the centre of mass velocity
		double RML=SPM[5][LS][MS]/SP[5][MS];
		double RMM=SPM[5][LS][MS]/SP[5][LS];
		double A,B,C,D;
		double OC,SC;
		
		for (int K=1;K<=3;K++)
		{
			VCCM[K]=RML*PV[K][L]+RMM*PV[K][M];
		}
		
		//VCCM defines the components of the centre of mass velocity (eqn 2.1)
		if (Math.abs(SPM[4][LS][MS]-1.)<1.E-3)
		{
			//use the VHS logic
			B=2.*RF(0)-1.0;	//B is the cosine of a random elevation angle
			A=Math.sqrt(1.-B*B);
			VRCP[1]=B*VR;
			C=2.*PI*RF(0); //C is a random azimuth angle
			VRCP[2]=A*Math.cos(C)*VR;
			VRCP[3]=A*Math.sin(C)*VR;		
		}
		else
		{
			//use the VSS logic
			B=2.*(Math.pow(RF(0),SPM[4][LS][MS]))-1.;
			//B is the cosine of the deflection angle for the VSS model (eqn (11.8)
			A=Math.sqrt(1.-B*B);
			C=2.*PI*RF(0);
			OC=Math.cos(C);
			SC=Math.sin(C);
			D=Math.sqrt(VRC[2]*VRC[2]+VRC[3]*VRC[3]);
			if (D>1.E-6) {
				VRCP[1]=B*VRC[1]+A*SC*D;
				VRCP[2]=B*VRC[2]+A*(VR*VRC[3]*OC-VRC[1]*VRC[2]*SC)/D;
				VRCP[3]=B*VRC[3]-A*(VR*VRC[2]*OC+VRC[1]*VRC[3]*SC)/D;
						}
			else {				
				VRCP[1]=B*VRC[1];
				VRCP[2]=A*OC*VRC[1];
				VRCP[3]=A*SC*VRC[1];
				}
			//the post-collision rel. velocity components are based on eqn (2.22)
			//VRCP(1 to 3) are the components of the post-collision relative vel.
		}
		
		for (int K=1;K<=3;K++)
		{
			PV[K][L]=VCCM[K]+VRCP[K]*RMM;
			PV[K][M]=VCCM[K]-VRCP[K]*RML;
		}
	}
	
	//samples the molecules in the flow
	public void SAMPLE0()
	{
		NSMP=NSMP+1;
		for (NN=1;NN<=MNSG;NN++)
			for (int N=1;N<=MNC;N++)
			{
				L=IC[2][N][NN];
				if (L>0)
				{
					for (int J=1;J<=L;J++)
					{
						int K=IC[1][N][NN]+J;
						M=IR[K];
						int I=IPS[M];
						CS[1][N][I]=CS[1][N][I]+1;
						for (int LL=1;LL<=3;LL++)
						{
							CS[LL+1][N][I]=CS[LL+1][N][I]+PV[LL][M];
							CS[LL+4][N][I]=CS[LL+4][N][I]+PV[LL][M]*PV[LL][M];
						}
					}
				}             
			}
	}

	//output a progressive set of results to file DSMC0.OUT
	public void OUT0()
	{
		double VEL[] = new double[3+1];
		double TCOL[][] = new double[MNSP+1][MNSP+1];
		double SVEL[][] = new double[3+1][MNC+1];

		PrintWriter pw=null;
		try
		{
			pw = new PrintWriter(new FileWriter("DSMC0.OUT"));
		}
		catch (Exception e)
		{
			System.err.println("Error opening file");
		}
		
		pw.printf(" FROM ZERO TIME TO TIME %f\n",TIME);
		pw.println(" COLLISIONS:-");
		for (M=1;M<=MNSP;M++)
		{
			for (L=1;L<=MNSP;L++)
				pw.printf("%9d\t",(int)COL[M][L]);
			pw.println();
		}

		pw.println(" TOTAL NUMBER OF SAMPLES "+NSMP);
		pw.println(NM+ " MOLECULES");
		pw.println(MOVT+" TOTAL MOLECULAR MOVES");
		pw.println((int)(SELT) + " SELECTIONS " + 
			   (int)(NCOL) + " COLLISION EVENTS, RATIO  "+ NCOL/(double)SELT);
		if (NCOL>0) pw.println (" MEAN COLLISION SEPARATION "+ SEPT/(double)NCOL);

		pw.println ("SAMPLES");
		pw.println (" CELL     N SP 1    N SP 2     ETC ");
		for (int N=1;N<=MNC;N++)
		{
			pw.printf("%6d\t",N);
			for (L=1;L<=MNSP;L++)
				pw.printf("%6.0f\t",CS[1][N][L]);
			pw.println();
		}

		pw.println(" FLOWFIELD PROPERTIES");
		pw.println(" CELL   X COORD     DENSITY       U         V        W         TEMP");

		double TOT=0.;
		double SMU[]=new double[3+1];

		//first the mixture properties
		for (int N=1;N<=MNC;N++)
		{
			double A=FNUM/(CG[3][N]*NSMP);
			double SN=0.;
			double SM=0.;
			for (int K=1;K<=3;K++)
				SMU[K]=0.;

			double SMCC=0.;
			for (L=1;L<=MNSP;L++)
			{
				SN=SN+CS[1][N][L];		//SN is the number sum
				SM=SM+SP[5][L]*CS[1][N][L];	//SM is the sum of molecular masses
				for (int K=1;K<=3;K++)
					SMU[K]=SMU[K]+SP[5][L]*CS[K+1][N][L];	//SMU(1 to 3) are the sum of mu, mv, mw
				SMCC=SMCC+(CS[5][N][L]+CS[6][N][L]+CS[7][N][L])*SP[5][L];	//SMCC is the sum of m(u**2+v**2+w**2)
			}

			double DENN=SN*A;	//DENN is the number density, see eqn (1.34)
			double DEN=DENN*SM/SN;	//DEN is the density, see eqn (1.42)
			for (int K=1;K<=3;K++)
			{
				//VEL and SVEL are the stream velocity components, see eqn (1.43)
				VEL[K]=SMU[K]/SM;
				SVEL[K][N]=VEL[K];
			}

			double UU=VEL[1]*VEL[1]+VEL[2]*VEL[2]+VEL[3]*VEL[3];
			double TT=(SMCC-SM*UU)/(3.*BOLTZ*SN);	//TT is the temperature, see eqn (1.51)
			TOT=TOT+TT;
			double XC=0.5*(CG[1][N]+CG[2][N]);	//XC is the x coordinate of the midpoint of the cell

			pw.printf("%6d %9.4f %9.4g %9.4f %9.4f %9.4f %9.4f\n",N,XC,DEN,VEL[1],VEL[2],VEL[3],TT);
		}

		pw.println();

		for (L=1;L<=MNSP;L++)
		{
			//now the properties of the separate species
			pw.println();
			pw.printf(" SPECIES %d\n",L);
			pw.println(" CELL   X COORD      N DENS     DENSITY U DIF VEL V DIF VEL W DIF VEL	TEMP ");
			for (int N=1;N<=MNC;N++)
			{
				double A=FNUM/(CG[3][N]*NSMP);
				double DENN=CS[1][N][L]*A;	//DENN is the partial number density
				double DEN=SP[5][L]*DENN;	//DEN is the partial density, see eqn (1.13)
				for (int K=1;K<=3;K++)
					VEL[K]=CS[K+1][N][L]/CS[1][N][L];	//VEL defines the average velocity of the species L molecules

				double UU=VEL[1]*VEL[1]+VEL[2]*VEL[2]+VEL[3]*VEL[3];
				double TT=(SP[5][L]/(3.*BOLTZ))*((CS[5][N][L]+CS[6][N][L]+CS[7][N][L])/CS[1][N][L]-UU);	//TT is the temperature, see eqn (1.29)

				for (int K=1;K<=3;K++)
					VEL[K]=VEL[K]-SVEL[K][N];	//VEL now defines the diffusion velocity of species L, see eqn (1,45)

				double XC=0.5*(CG[1][N]+CG[2][N]);
				pw.printf("%6d %9.4f %9.4g %9.4g %9.4f %9.4f %9.4f %9.4f\n",N,XC,DENN,DEN,VEL[1],VEL[2],VEL[3],TT);
			}
		}
	
		//compare with the theoretical collision number for actual temperarure
		double AVTMP=TOT/MNC;

		pw.println();
		pw.printf(" AVERAGE TEMPERATURE %g\n",AVTMP);
		pw.println();

		pw.println (" COLLISIONS:-");
		for (L=1;L<=MNSP;L++)
		{
			for (M=1;M<=MNSP;M++)
				pw.printf("%8.0f\t",COL[M][L]);
			pw.println();
		}

		pw.println();
		pw.println(" RATIO OF COLLISION NUMBER TO THEORETICAL VALUE");
		pw.println();

		for (M=1;M<=MNSP;M++)
		{
			double SML=0;
			for (int K=1;K<=MNC;K++)
				SML=SML+CS[1][K][M]/NSMP;

			for (L=1;L<=MNSP;L++)
			{
				double SLL=0;
				for (int K=1;K<=MNC;K++)
					SLL=SLL+CS[1][K][L]/NSMP;

				//TCOL is the equilibrium collision rate, see eqn (4.78)
				TCOL[M][L]=2.*TIME*FNUM*SML*SLL*(1./(CG[2][MNC]-CG[1][1]))
				   *SPM[1][M][L]*Math.pow((AVTMP/SPM[2][M][L]),1.-SPM[3][M][L])
				   *Math.sqrt(2.*BOLTZ*SPM[2][L][M]/(PI*SPM[5][L][M]));
			}
		}

		for (L=1;L<=MNSP;L++)
		{
			for (M=1;M<=MNSP;M++)
				pw.printf("%8.4f\t",COL[M][L]/TCOL[M][L]);
			pw.println();
		}

		pw.close();
	}

	//generates two random velocity components U an V in an equilibrium
	//gas with most probable speed VMP  (based on eqns (C10) and (C12))
	public double RVELC(double V, double VMP)
	{
		double A=Math.sqrt(-Math.log(RF(0)));	/*this was LOG, is that LN or LOG10?*/
		double B=6.283185308*RF(0);

		/*this is what bird does but then uses only U...*/	  
		double U=A*Math.sin(B)*VMP;
		V=A*Math.cos(B)*VMP;
		return U;
	}

	/*  FUNCTION GAM(X)
	//calculates the Gamma function of X
	* tested against built in function in Excel...
	*/
	double GAM(double X)
	{
		double A=1.;
		double Y=X;
		if (Y<1.0)
			A=A/Y;
		else
		{
			do {
				Y=Y-1;
				if (Y>=1.)
					A=A*Y;	
			} while (Y>=1.);
				
		}
		
		double GAM=A*(1.-0.5748646*Y+0.9512363*Y*Y-0.6998588*Y*Y*Y+
				0.4245549*Y*Y*Y*Y-0.1010678*Y*Y*Y*Y*Y);
		return GAM;
	}

	/*random number generator...*/
	/*this seems to work fine since no difference in result when replaced with System.random()..*/
	int MA[] = new int[55+1];
	int INEXT,INEXTP;
	int IFF;
	Random rng=new Random(0);
	
	public double RF0(int IDUM)
	{
		return rng.nextDouble();
	}
	
	public double RF(int IDUM)
	{
		//generates a uniformly distributed random fraction between 0 and 1
		//--IDUM will generally be 0, but negative values may be used to
		//----re-initialize the seed
		int MBIG=1000000000;
		int MSEED=161803398;
		int MZ=0;
		double FAC=1.0E-9;
		double RF;
		int MJ;
		int MK;
		
		if (IDUM<0 || IFF==0) 
		{
			IFF=1;
			MJ=MSEED-Math.abs(IDUM);
			MJ=MJ%MBIG;
			MA[55]=MJ;
			MK=1;
			
			for (int I=1;I<=54;I++)
			{
				int II=(21*I)%55;
				MA[II]=MK;
				MK=MJ-MK;
				if (MK<MZ) MK=MK+MBIG;
				
				MJ=MA[II];
			}
	
			for (int K=1;K<=4;K++)
				for (int I=1;I<=55;I++)
				{
					MA[I]=MA[I]-MA[1+(I+30)%55];
					if (MA[I]<MZ) MA[I]=MA[I]+MBIG;
				}
	
			INEXT=0;
			INEXTP=31;
		} /*if*/
		
		while (true)
		{
			INEXT=INEXT+1;
			if (INEXT==56) INEXT=1;
			INEXTP=INEXTP+1;
			if (INEXTP==56) INEXTP=1;
			MJ = MA[INEXT]-MA[INEXTP];
			if (MJ<MZ) MJ=MJ+MBIG;
			MA[INEXT]=MJ;
			RF=MJ*FAC;
			if (RF>1.E-8 && RF<0.99999999) return RF;
		}
}
  
  
/* 
//defines the data for a particular run of DSMC0.FOR
**/
	public void DATA0()
	{
		ND=1.E20;        //FND  is the number densty
		TMP=300.;           //--FTMP is the temperature
		FSP[1]=.6;         ////FSP(N) is the number fraction of species N
		FSP[2]=.2;
		FSP[3]=.1;
		FSP[4]=.08;
		FSP[5]=.02;

		FNUM=1.0E17;      //FNUM  is the number of real molecules represented by a simulated mol.
		DTM=.25E-4;       //--DTM is the time step
		NSC=8;            //--NSC is the number of sub-cells in each cell (8)
		XF=0.;            ////the simulated region is from x=XF to x=XR
		XR=1.;
		ISPD=0;           ////the cross-collision data is set to the mean values

		//SP(1,N) is the molecular diameter of species N
		//SP(2,N) is the reference temperature
		//SP(3,N) is the viscosity-temperatire index
		//SP(4,N) is the reciprocal of the VSS scattering parameter
		//SP(5,N) is the molecular mass of species N
		//ISP(N) is the group for species N
	  	SP[1][1]=3.5E-10;
		SP[2][1]=273.;
		SP[3][1]=0.75;
		SP[4][1]=1.;
		SP[5][1]=5.E-26;
		ISP[1]=1;
		SP[1][2]=4.E-10;
		SP[2][2]=273.;
		SP[3][2]=0.75;
		SP[4][2]=1.;
		SP[5][2]=4.5E-26;
		ISP[2]=1;
		SP[1][3]=3.E-10;
		SP[2][3]=273.;
		SP[3][3]=0.75;
		SP[4][3]=1.;
		SP[5][3]=2.5E-26;
		ISP[3]=1;
		SP[1][4]=3.E-10;
		SP[2][4]=273.;
		SP[3][4]=0.75;
		SP[4][4]=1.;
		SP[5][4]=2.E-26;
		ISP[4]=1;
		SP[1][5]=4.E-10;
		SP[2][5]=273.;
		SP[3][5]=0.75;
		SP[4][5]=1.;
		SP[5][5]=4.E-26;
		ISP[5]=1;

		NIS=4;	//NIS is the number of time steps between samples
		NSP=40;	//NSP is the number of samples between restart and output file updates
		NPT=500;	//NPT is the number of file updates to STOP		
	}
}