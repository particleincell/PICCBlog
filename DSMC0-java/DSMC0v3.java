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
 
public class DSMC0v3 
{
    public static void main(String[] args) 
    {
        new Simulation().run();
    }
}
class Simulation
{
	final int np_max=1000;		//MNM  is the maximum number of molecules
	final int NumCells=50;		//MNC  is the maximum number of sub
	final int NumSubCells=400;		//MNSC is the maximum number of sub-cells
	final int NumSpecies=5;		//MNSP is the maximum number of molecular species
	final int NumSpGroups=1;		//MNSG is the number of species groups for collision sampling

	final int SampleNumIt=4;	//NIS is the number of time steps between samples 4
	final int RestartNumIt=40;	//NSP is the number of samples between restart and output file updates 40
	final int NumFileUpdates=500;	//NPT is the number of file updates to STOP 500

	double COL[][] = new double[NumSpecies][NumSpecies];	//COL(M,N) is the number of collisions between species N-M molecules
	double TotNumMoves;    //MOVT the total number of molecular moves
	long TotNumCols;    //NCOL is the total number of collisions
	double TotNumSelections;	//SELT the total number of pair selections
	double SumColPairSeparations;	//SEPT the sum of collision pair separations
	SampledInfo sampled_info[][] = new SampledInfo[NumCells][NumSpecies];	//CS(N,M,L) sampled information on species L in cell M
												//---- N=1 number sum
												//----N=2,3,4 sum of u,v,w
												//----N=5,6,7 sum of u*u,v*v,w*w
	class SampledInfo
	{
		double num_sum;
		double v_sum[]=new double[3];
		double v2_sum[] = new double[3];
	}
	/*MOLS block*/
	int np;                          //NM is the number of molecules
	
	class Part
	{
		double x;						//PP(M)  is the x coordinate molecule M					
		double v[] = new double[3];		//PV(1 to 3,M)  u,v,w velocity components of molecule M
	}
	
	Part part[] = new Part[np_max];
	
	int SubCellNumber[] = new int[np_max];     //IPL(M) sub-cell number for molecule M
	int SpeciesCodeNumber[] = new int[np_max];     //IPS(M) species code number
	int CrossRef[] = new int[np_max];      //IR(M)  cross-reference array (molecule numbers in order of sub-cells)

	/*cells block*/
	double CellVolume[] = new double[NumCells];      //CC(M) is the cell volume
	CellGeom cell_geom[] = new CellGeom[NumCells];   //--CG(N,M) is the geometry related information on cell M
										  //----N=1 the minimum x coordinate
										  //----N=2 the maximum x coordinate
										  //----N=3 the cell width
	class CellGeom
	{
		double x_min,x_max;
		double cell_width;
	}
	
	CellInfo cell_info[][] = new CellInfo[NumCells][NumSpGroups]; //--IC(N,M,L) information on the molecules of species group L in cell M
				
	class CellInfo
	{
		int start;								//----N=1 (start address -1) of the molecule numbers in the array IR
		int count;								//----N=2 the number of molecules in the cell
	}
	
	int Sub2Cell[] = new int[NumSubCells];    //ISC(M) the cell in which the sub-cell lies
	SpGroupColData spg_col_data[][][] = new SpGroupColData[NumCells][NumSpGroups][NumSpGroups]; //--CCG(N,M,L,K) is for collisions between species groups L-K in cell M
	class SpGroupColData
	{
		double sigma_g_max;				//----N=1 is the maximum value of (relative speed)*(coll. cross-section)
		double rem;								//----N=2 is the remainder when the selection number is rounded
	}
	
	SpGroupCellInfo spg_cell_info[][] = new SpGroupCellInfo[NumSubCells][NumSpGroups]; //--ISCG(N,M,L) is the information on species group L in sub-cell M
	class SpGroupCellInfo
	{
		int start;			   //----N=1 (start address -1) of the molecule numbers in the array IR
		int count;			   //----N=2 the number of molecules in the sub-cell
	}
	
	GroupInfo group_info[] = new GroupInfo[NumSpGroups];    //--IG(2,M) information on group L molecules
	class GroupInfo
	{	
		int start;			//----N=1 (start address -1) of the molecule numbers in the array IR
		int count;			//----N=2 the number of molecules in the cell
	}
	

	/*gas block*/
	class SpeciesInfo
	{
		double diam;
		double ref_temp;
		double visc_temp_index;
		double vss_scat_inv;
		double mass;
	};
	
	SpeciesInfo species_info[] = new SpeciesInfo[NumSpecies];   //--SP(N,M) information on species M
										   //----N=1 the reference cross-section (diameter in the data)
										   //----N=2 the reference temperature
										   //----N=3 the viscosity-temperature power law
										   //----N=4 the reciprocal of the VSS scattering parameter
										   //----N=5 the molecular mass
	SpeciesInteraction species_interaction[][] = new SpeciesInteraction[NumSpecies][NumSpecies]; ////SPM(N,M,L) information on the interaction between L-M molecules
	class SpeciesInteraction
	{
		double sigma;			//--N=1  the reference cross-section (diameter in the data)
		double ref_temp;		//--N=2  the reference temperature
		double visc_temp_index;				//--N=3  the viscosity-temperature power law
		double vss_scat_inv;				//--N=4  the reciprocal of the VSS scattering parameter
		double reduced_mass;	//--N=5  the reduced mass
		double gamma;			//--N=6  the Gamma function of (5/2 - viscosity-temperature power law)
	}
	
	int Species2SamplingGroup[] = new int[NumSpecies];    //ISP(M) the colision sampling group in which species M lies


	/*SAMP block*/
	double TIME;            //TIME time
	int FileOutputIt;         ////NPR the number of output/restart file update cycles
	int TotNumSamples;        ////NSMP the total number of samples
	double StreamNumDen;         ////FND the stream number density
	double StreamTemp;        ////FTMP the stream temperature
	double SpeciesFraction[] = new double[NumSpecies];    ////FSP(M) the fraction of species M in the stream
	int UnlikeColsFlag;		////ISPD relates to the setting of data for colls. between unlike mols.
					////--set to 0 if data is set automatically to the mean values
					////--set to 1 if the values are set explicitly in the data

	/*comp block*/
	double Spwt;    ////FNUM  is the number of real molecules represented by a simulated mol.
	double DeltaT;     ////DTM is the time step
	int SteadyFlowEstimateCycles;     ////NPS is the estimated number of samples to steady flow
	
	/*geom block*/
	double CellWidth;      ////CW is the cell width
	int SubCellsPerCell;     ////NSC is the number of sub-cells per cell
	double XMIN;      ////XF is the minimum x coordinate
	double XMAX;      ////XR is the maximum x coordinate

	/*const block*/
	double PI;      ////PI is pi and SPI is the square root of pi
	double SPI;
	double BOLTZ;   ////BOLTZ is the Boltzmann constant
	
	/*elast block*/  
	double VRC[]=new double[3];
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
		//System.out.printf("RF %.10g %.10g %.10g %.10g\n",RF(0),RF(0),RF(0),RF(0));
		INIT0();
		SAMPLE_INIT();

		while (FileOutputIt<NumFileUpdates)
		{
			FileOutputIt++;
			for (int restart_it=1;restart_it<=RestartNumIt;restart_it++)
			{
				for (int sample_it=1;sample_it<=SampleNumIt;sample_it++)
				{
					TIME=TIME+DeltaT;
	
					System.out.printf(" DSMC0:- Move %d %d   of  %d %d\t%d Collisions\n",sample_it,restart_it,SampleNumIt,RestartNumIt,TotNumCols);
	
					MOVE0();
					INDEXM();
					COLLM();
					
				}

				SAMPLE0();
			}	/*for*/

			/*save restart...*/
			
			System.out.println(TotNumCols+" Collisions");
	
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
		if (NumSpecies==1) UnlikeColsFlag=0;
		for (int N=0;N<NumSpecies;N++)
		{
			for (M=0;M<NumSpecies;M++)
			{
				if ((UnlikeColsFlag==0)||(N==M)) 
				{
					//--the collision cross section is assumed to be given by eqn (1.35)
					species_interaction[N][M].sigma=0.25*PI*(species_info[N].diam+species_info[M].diam)*(species_info[N].diam+species_info[M].diam);
				
					//--mean values are used for ISPD=0
					species_interaction[N][M].ref_temp=0.5*(species_info[N].ref_temp+species_info[M].ref_temp);
					species_interaction[N][M].visc_temp_index=0.5*(species_info[N].visc_temp_index+species_info[M].visc_temp_index);
					species_interaction[N][M].vss_scat_inv=0.5*(species_info[N].vss_scat_inv+species_info[M].vss_scat_inv);
				}
				else
				{
					//--the cross-collision diameter is converted to the cross-section
					species_interaction[N][M].sigma=PI*species_interaction[N][M].sigma*species_interaction[N][M].sigma;
				}	/*end if*/
				
				//*--the reduced mass is defined in eqn (2.7)
				species_interaction[N][M].reduced_mass=(species_info[N].mass/(species_info[N].mass+species_info[M].mass))*species_info[M].mass;
				species_interaction[N][M].gamma=GAM(2.5-species_interaction[N][M].visc_temp_index);
			}	/*for*/
		}	/*for*/

		//*--initialise variables
		TIME=0.;
		np=0;
	  
		cell_geom[0].x_min=XMIN;
		CellWidth=(XMAX-XMIN)/NumCells;
		for (M=0;M<NumCells;M++)
		{
			if (M>0) cell_geom[M].x_min=cell_geom[M-1].x_max;
			cell_geom[M].x_max=cell_geom[M].x_min+CellWidth;
			cell_geom[M].cell_width=CellWidth;
			CellVolume[M]=CellWidth;
			for (L=0;L<NumSpGroups;L++)
				for (int K=0;K<NumSpGroups;K++)
				{
					//*--the maximum value of the (rel. speed)*(cross-section) is set to a
					//*--reasonable, but low, initial value and will be increased as necessary
					spg_col_data[M][L][K].sigma_g_max=species_interaction[0][0].sigma*300.*Math.sqrt(StreamTemp/300.);
					spg_col_data[M][L][K].rem=RF(0);
					
				}
		}	/*for*/

		//--set sub-cells
		for (int N=0;N<NumCells;N++)
			for (M=0;M<SubCellsPerCell;M++)
			{
				L=N*SubCellsPerCell+M;
				Sub2Cell[L]=N;
			}

		//--generate initial gas in equilibrium at temperature FTMP
		for (L=0;L<NumSpecies;L++)
		{
			double REM=0;
			//--VMP is the most probable speed in species L, see eqns (4.1) and (4.7)
			double VMP=Math.sqrt(2.*BOLTZ*StreamTemp/species_info[L].mass);	
			for (int N=0;N<NumCells;N++)
			{
				//--A is the number of simulated molecules of species L in cell N to
				//--simulate the required concentrations at a total number density of FND
				double A=StreamNumDen*cell_geom[N].cell_width*SpeciesFraction[L]/Spwt+REM;
				if (N<NumCells) 
				{
					MM=(int)A;		//! added (int)
					REM=(A-MM);		//--the remainder REM is carried forward to the next cell
				} else
				{
					MM=(int)(A+0.5);	/*this was NINT, round to nearest whole number*/
				}
				
				for (M=0;M<MM;M++)
				{
					if (np<np_max) /*this was <=*/
					{
						SpeciesCodeNumber[np]=L;
						part[np].x=cell_geom[N].x_min+RF(0)*(cell_geom[N].x_max-cell_geom[N].x_min);
						/*!added (int)*/
						SubCellNumber[np]=(int)((part[np].x-cell_geom[N].x_min)*
								(SubCellsPerCell-0.001)/cell_geom[N].cell_width+SubCellsPerCell*N);
						//--species, position, and sub-cell number have been set
					  
						//--velocity components have been set
						for (int i=0;i<3;i++)
							part[np].v[i] = RVELC(A,VMP);
						
						//--round-off error could have taken NM to MNM+1
						np=np+1;
					}	/*end if*/
				} /*for*/
			}	/*for*/
		}
		
		System.out.printf ("%d MOLECULES\n",np);
	}

/*   SAMPI0.FOR
*
**/
	public void SAMPLE_INIT()
	{
		FileOutputIt=0;
		TotNumCols=0;
		TotNumSamples=0;
		TotNumMoves=0.;
		TotNumSelections=0.;
		SumColPairSeparations=0.;
		
		for (L=0;L<NumSpecies;L++)
			for (int N=0;N<NumCells;N++)
			{
				sampled_info[N][L].num_sum=1.E-6;
				for (int i=0;i<3;i++)
				{
					sampled_info[N][L].v_sum[i]=0.0;
					sampled_info[N][L].v2_sum[i]=0.0;
				}
			}
	
		for (M=0;M<NumSpecies;M++)
			for (int N=0;N<NumSpecies;N++)
				COL[M][N]=0.;
	}
	
	public void MOVE0()
	{
		for (int N=0;N<np;N++)
		{
			TotNumMoves++;
			MSC=SubCellNumber[N];
			int MC=Sub2Cell[MSC];			//MC is the initial cell number
			double XI=part[N].x;
			double DX=part[N].v[0]*DeltaT;
			double X=XI+DX;
			//molecule N at XI is moved by DX to X
			if (X<XMIN)
			{
				//specular reflection from the minimum x boundary at x=XF (eqn (11.7))
				X=2.*XMIN-X;
				part[N].v[0]=-part[N].v[0];
			}
			
			if (X>XMAX)
			{
				//specular reflection from the maximum x boundary at x=XR (eqn (11.7))
				X=2.*XMAX-X;
				part[N].v[0]=-part[N].v[0];
			}
			
			if (X<cell_geom[MC].x_min||X>cell_geom[MC].x_max)
			{
				//the molecule has moved from the initial cell
				MC = (int)((X-XMIN)/CellWidth+(0.99999-1));
				if (MC<=0) MC=0;
				//MC is the new cell number (note avoidance of round-off error)
			}
			
			MSC=(int)(((X-cell_geom[MC].x_min)/cell_geom[MC].cell_width)*(SubCellsPerCell-.001)+SubCellsPerCell*MC);
			//MSC is the new sub-cell number
			SubCellNumber[N]=MSC;
			part[N].x=X;
		}
	}
  
	//the NM molecule numbers are arranged in order of the molecule groups
	//and, within the groups, in order of the cells and, within the cells,
	//in order of the sub-cells
        public void INDEXM()
	{
		for (MM=0;MM<NumSpGroups;MM++)
		{
			group_info[MM].count=0;
			for (NN=0;NN<NumCells;NN++)
				cell_info[NN][MM].count=0;
      
			for (NN=0;NN<NumSubCells;NN++)
				spg_cell_info[NN][MM].count=0;
		}
	  
		for (int N=0;N<np;N++)
		{
			LS=SpeciesCodeNumber[N];
			int MG=Species2SamplingGroup[LS];
			group_info[MG].count++;
			MSC=SubCellNumber[N];
			spg_cell_info[MSC][MG].count++;
			int MC=Sub2Cell[MSC];
			cell_info[MC][MG].count++;
		}
		//--number in molecule groups in the cells and sub-cells have been counte
		
		M=0;
		for (L=0;L<NumSpGroups;L++)
		{
			group_info[L].start=M;	//--the (start address -1) has been set for the groups
			M=M+group_info[L].count;
		}
		
		for (L=0;L<NumSpGroups;L++)
		{
			//--the (start address -1) has been set for the cells
			M=group_info[L].start;
			for (int N=0;N<NumCells;N++)
			{
				cell_info[N][L].start=M;
				M=M+cell_info[N][L].count;
			}
					
			//--the (start address -1) has been set for the sub-cells
			M=group_info[L].start;
			for (int N=0;N<NumSubCells;N++)
			{
				spg_cell_info[N][L].start=M;
				M=M+spg_cell_info[N][L].count;
				spg_cell_info[N][L].count=0;
			}
		}
		
		for (int N=0;N<np;N++)
		{
			//*--the molecule number N has been set in the cross-reference array
			LS=SpeciesCodeNumber[N];
			int MG=Species2SamplingGroup[LS];
			MSC=SubCellNumber[N];
			int K=spg_cell_info[MSC][MG].start+spg_cell_info[MSC][MG].count;
			CrossRef[K]=N;
			spg_cell_info[MSC][MG].count++;
		}
                
    //      for (int i=0;i<CrossRef.length;i++)
    //          System.out.printf("%d ",CrossRef[i]);
    //      System.out.printf("\n");
	}
	
	//calculates collisions appropriate to DTM in a monatomic gas mixture
	//VRC(3) are the pre-collision components of the relative velocity
	public void COLLM()
	{
		for (int N=0;N<NumCells;N++)
		{
			//--consider collisions in cell N
			for (NN=0;NN<NumSpGroups;NN++)
			{
					for (MM=0;MM<NumSpGroups;MM++)
					{
						double SN=0.;
						for (int K=0;K<NumSpecies;K++)
						  if (Species2SamplingGroup[K]==MM) SN=SN+sampled_info[N][K].num_sum;
			
						//*--AVN is the average number of group MM molecules in the cell
						double AVN;
						if (SN>1.0)
						  AVN=SN/(float)(TotNumSamples);
						else
						  AVN=cell_info[N][MM].count;
								
						//*--ASEL is the number of pairs to be selected, see eqn (11.5)
						double ASEL=0.5*cell_info[N][NN].count*AVN*Spwt*
								spg_col_data[N][NN][MM].sigma_g_max*DeltaT/CellVolume[N]+spg_col_data[N][NN][MM].rem;
					
						int NSEL=(int)ASEL;
						spg_col_data[N][NN][MM].rem=ASEL-NSEL;
	//		System.out.printf("%g %d\n",ASEL,NSEL);			
						if (NSEL>0)
						{
						//*--if there are insufficient molecules to calculate collisions,
						//*--the number NSEL is added to the remainer CCG(2,N,NN,MM)
							if (((NN!=MM)&&(cell_info[N][NN].count<1 || cell_info[N][MM].count<1))
								||((NN==MM)&&(cell_info[N][NN].count<2))) 
							{
								spg_col_data[N][NN][MM].rem=spg_col_data[N][NN][MM].rem+NSEL;
							}
							else
							{
								double CVR_MAX=spg_col_data[N][NN][MM].sigma_g_max;
								TotNumSelections=TotNumSelections+NSEL;
								for (int ISEL=0;ISEL<NSEL;ISEL++)
								{
									/*note, this sets L and M as well as LS and MS*/
									SELECT(N);
		
									//*--if necessary, the maximum product in CVM is upgraded
									if (CVR>CVR_MAX) CVR_MAX=CVR;
							  
									//*--the collision is accepted with the probability of eqn (11.6)
									if (RF(0)<CVR/spg_col_data[N][NN][MM].sigma_g_max) 
									{
										TotNumCols++;
										SumColPairSeparations=SumColPairSeparations+Math.abs(part[L].x-part[M].x);
										COL[LS][MS]++;
										COL[MS][LS]++;
			
										ELASTIC();
									}/*if*/
								}	/*for*/
								spg_col_data[N][NN][MM].sigma_g_max=CVR_MAX;
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
		int K=(int)(RF(0)*(cell_info[N][NN].count-0.0001))+cell_info[N][NN].start;
		L=CrossRef[K];
		M = L;
		
		//the first molecule L has been chosen at random from group NN in cell
		do {
			MSC=SubCellNumber[L];
			if ((NN==MM && spg_cell_info[MSC][MM].count==1) || (NN!=MM&&spg_cell_info[MSC][MM].count==0))
			{
				//if MSC has no type MM molecule find the nearest sub-cell with one
				int NST=1;	/*these were 1, keep 1 or 0?*/
				int NSG=1;
				do 
				{
					int INC=NSG*NST;
					NSG=-NSG;
					NST++;
					MSC=MSC+INC;
				}	while (MSC<0 || MSC>=NumSubCells ||
						   Sub2Cell[MSC] != N || spg_cell_info[MSC][MM].count<1);
			}
			
			//the second molecule M is now chosen at random from the group MM
			//molecules that are in the sub-cell MSC
			K=(int)(RF(0)*(spg_cell_info[MSC][MM].count-0.0001))+spg_cell_info[MSC][MM].start;
			M=CrossRef[K];
			//choose a new second molecule if the first is again chosen
		} while (L==M);		
			
		for (int i=0;i<3;i++)
		{
			VRC[i] = part[L].v[i]-part[M].v[i];
			//VRC(1 to 3) are the components of the relative velocity
			
		}

		VRR=VRC[0]*VRC[0]+VRC[1]*VRC[1]+VRC[2]*VRC[2];
		VR=Math.sqrt(VRR);//VR is the relative speed
		LS=SpeciesCodeNumber[L];
		MS=SpeciesCodeNumber[M];
		CVR=VR*species_interaction[LS][MS].sigma*
				Math.pow(2.*BOLTZ*species_interaction[LS][MS].ref_temp/(species_interaction[LS][MS].reduced_mass*VRR),
				species_interaction[LS][MS].visc_temp_index-0.5)/species_interaction[LS][MS].gamma;
		//the collision cross-section is based on eqn (4.63)
	}
  
	//generates the post-collision velocity components
	public void ELASTIC()
	{
		double VRCP[] = new double[3];		//VRCP(3) are the post-collision components of the relative velocity
		double VCCM[] = new double[3];		//VCCM(3) are the components of the centre of mass velocity
		double RML=species_interaction[LS][MS].reduced_mass/species_info[MS].mass;
		double RMM=species_interaction[LS][MS].reduced_mass/species_info[LS].mass;
		double A,B,C,D;
		double OC,SC;
		
		for (int i=0;i<3;i++)
		{
			VCCM[i]=RML*part[L].v[i]+RMM*part[M].v[i];
		}
		
		//VCCM defines the components of the centre of mass velocity (eqn 2.1)
		if (Math.abs(species_interaction[LS][MS].vss_scat_inv-1.)<1.E-3)
		{
			//use the VHS logic
			B=2.*RF(0)-1.0;	//B is the cosine of a random elevation angle
			A=Math.sqrt(1.-B*B);
			VRCP[0]=B*VR;
			C=2.*PI*RF(0); //C is a random azimuth angle
			VRCP[1]=A*Math.cos(C)*VR;
			VRCP[2]=A*Math.sin(C)*VR;		
		}
		else
		{
			//use the VSS logic
			B=2.*(Math.pow(RF(0),species_interaction[LS][MS].vss_scat_inv))-1.;
			//B is the cosine of the deflection angle for the VSS model (eqn (11.8)
			A=Math.sqrt(1.-B*B);
			C=2.*PI*RF(0);
			OC=Math.cos(C);
			SC=Math.sin(C);
			D=Math.sqrt(VRC[1]*VRC[1]+VRC[2]*VRC[2]);
			if (D>1.E-6) {
				VRCP[0]=B*VRC[0]+A*SC*D;
				VRCP[1]=B*VRC[1]+A*(VR*VRC[2]*OC-VRC[0]*VRC[1]*SC)/D;
				VRCP[2]=B*VRC[2]-A*(VR*VRC[1]*OC+VRC[0]*VRC[2]*SC)/D;
						}
			else {				
				VRCP[0]=B*VRC[1];
				VRCP[1]=A*OC*VRC[1];
				VRCP[2]=A*SC*VRC[1];
				}
			//the post-collision rel. velocity components are based on eqn (2.22)
			//VRCP(1 to 3) are the components of the post-collision relative vel.
		}
		
		for (int i=0;i<3;i++)
		{
			part[L].v[i]=VCCM[i]+VRCP[i]*RMM;
			part[M].v[i]=VCCM[i]-VRCP[i]*RML;
		}
	}
	
	//samples the molecules in the flow
	public void SAMPLE0()
	{
		TotNumSamples=TotNumSamples+1;
		for (NN=0;NN<NumSpGroups;NN++)
			for (int N=0;N<NumCells;N++)
			{
				L=cell_info[N][NN].count;
				if (L>0)
				{
					for (int J=0;J<L;J++)
					{
						int K=cell_info[N][NN].start+J;
						M=CrossRef[K];
						int I=SpeciesCodeNumber[M];
						sampled_info[N][I].num_sum++;
						for (int i=0;i<3;i++)
						{
							sampled_info[N][I].v_sum[i]+=part[M].v[i];
							sampled_info[N][I].v2_sum[i]+=part[M].v[i]*part[M].v[i];
						}
					}
				}             
			}
	}

	//output a progressive set of results to file DSMC0.OUT
	public void OUT0()
	{
		double VEL[] = new double[3];
		double TCOL[][] = new double[NumSpecies][NumSpecies];
		double SVEL[][] = new double[3][NumCells];

		PrintWriter pw=null;
		try
		{
			pw = new PrintWriter(new FileWriter("DSMC1.OUT"));
		}
		catch (Exception e)
		{
			System.err.println("Error opening file");
		}
		
		pw.println(" FROM ZERO TIME TO TIME "+TIME);
		pw.println(" COLLISIONS:-");
		for (M=0;M<NumSpecies;M++)
		{
			for (L=0;L<NumSpecies;L++)
				pw.printf("%9d\t",(int)COL[M][L]);
			pw.println();
		}

		pw.println(" TOTAL NUMBER OF SAMPLES "+TotNumSamples);
		pw.println(np+ " MOLECULES");
		pw.println(TotNumMoves+" TOTAL MOLECULAR MOVES");
		pw.println((int)(TotNumSelections) + " SELECTIONS " + 
			   (int)(TotNumCols) + " COLLISION EVENTS, RATIO  "+ TotNumCols/(double)TotNumSelections);
		if (TotNumCols>0) pw.println (" MEAN COLLISION SEPARATION "+ SumColPairSeparations/(double)TotNumCols);

		pw.println ("SAMPLES");
		pw.println (" CELL     N SP 1    N SP 2     ETC ");
		for (int N=0;N<NumCells;N++)
		{
			pw.printf("%6d\t",N);
			for (L=0;L<NumSpecies;L++)
				pw.printf("%6.0f\t",sampled_info[N][L].num_sum);
			pw.println();
		}

		pw.println(" FLOWFIELD PROPERTIES");
		pw.println(" CELL   X COORD     DENSITY       U         V        W         TEMP");

		double TOT=0.;
		double SMU[]=new double[3];

		//first the mixture properties
		for (int N=0;N<NumCells;N++)
		{
			double A=Spwt/(cell_geom[N].cell_width*TotNumSamples);
			double SN=0.;
			double SM=0.;
			for (int K=0;K<3;K++)
				SMU[K]=0.;

			double SMCC=0.;
			for (L=0;L<NumSpecies;L++)
			{
				SN=SN+sampled_info[N][L].num_sum;		//SN is the number sum
				SM=SM+species_info[L].mass*sampled_info[N][L].num_sum;	//SM is the sum of molecular masses
				for (int K=0;K<3;K++)
					SMU[K]=SMU[K]+species_info[L].mass*sampled_info[N][L].v_sum[K];	//SMU(1 to 3) are the sum of mu, mv, mw
				SMCC=SMCC+(sampled_info[N][L].v2_sum[0]+sampled_info[N][L].v2_sum[1]+sampled_info[N][L].v2_sum[2])
						*species_info[L].mass;	//SMCC is the sum of m(u**2+v**2+w**2)
			}

			double DENN=SN*A;	//DENN is the number density, see eqn (1.34)
			double DEN=DENN*SM/SN;	//DEN is the density, see eqn (1.42)
			for (int K=0;K<=2;K++)
			{
				//VEL and SVEL are the stream velocity components, see eqn (1.43)
				VEL[K]=SMU[K]/SM;
				SVEL[K][N]=VEL[K];
			}

			double UU=VEL[0]*VEL[0]+VEL[1]*VEL[1]+VEL[2]*VEL[2];
			double TT=(SMCC-SM*UU)/(3.*BOLTZ*SN);	//TT is the temperature, see eqn (1.51)
			TOT=TOT+TT;
			double XC=0.5*(cell_geom[N].x_min+cell_geom[N].x_max);	//XC is the x coordinate of the midpoint of the cell

			pw.printf("%6d %9.4f %9.4g %9.4f %9.4f %9.4f %9.4f\n",N,XC,DEN,VEL[0],VEL[1],VEL[2],TT);
		}

		pw.println();

		for (L=0;L<NumSpecies;L++)
		{
			//now the properties of the separate species
			pw.println();
			pw.printf(" SPECIES %d\n",L);
			pw.println(" CELL   X COORD      N DENS     DENSITY U DIF VEL V DIF VEL W DIF VEL	TEMP ");
			for (int N=0;N<NumCells;N++)
			{
				double A=Spwt/(cell_geom[N].cell_width*TotNumSamples);
				double DENN=sampled_info[N][L].num_sum*A;	//DENN is the partial number density
				double DEN=species_info[L].mass*DENN;	//DEN is the partial density, see eqn (1.13)
				for (int K=0;K<3;K++)
					VEL[K]=sampled_info[N][L].v_sum[K]/sampled_info[N][L].num_sum;	//VEL defines the average velocity of the species L molecules

				double UU=VEL[0]*VEL[0]+VEL[1]*VEL[1]+VEL[2]*VEL[2];
				double TT=(species_info[L].mass/(3.*BOLTZ))*
						((sampled_info[N][L].v2_sum[0]+sampled_info[N][L].v2_sum[1]+sampled_info[N][L].v2_sum[2])/
						sampled_info[N][L].num_sum-UU);	//TT is the temperature, see eqn (1.29)

				for (int K=0;K<3;K++)
					VEL[K]=VEL[K]-SVEL[K][N];	//VEL now defines the diffusion velocity of species L, see eqn (1,45)

				double XC=0.5*(cell_geom[N].x_min+cell_geom[N].x_max);
				pw.printf("%6d %9.4f %9.4g %9.4g %9.4f %9.4f %9.4f %9.4f\n",N,XC,DENN,DEN,VEL[0],VEL[1],VEL[2],TT);
			}
		}
	
		//compare with the theoretical collision number for actual temperarure
		double AVTMP=TOT/NumCells;

		pw.println();
		pw.printf(" AVERAGE TEMPERATURE %g\n",AVTMP);
		pw.println();

		pw.println (" COLLISIONS:-");
		for (L=0;L<NumSpecies;L++)
		{
			for (M=0;M<NumSpecies;M++)
				pw.printf("%8.0f\t",COL[M][L]);
			pw.println();
		}

		pw.println();
		pw.println(" RATIO OF COLLISION NUMBER TO THEORETICAL VALUE");
		pw.println();

		for (M=0;M<NumSpecies;M++)
		{
			double SML=0;
			for (int K=0;K<NumCells;K++)
				SML=SML+sampled_info[K][M].num_sum/TotNumSamples;

			for (L=0;L<NumSpecies;L++)
			{
				double SLL=0;
				for (int K=0;K<NumCells;K++)
					SLL=SLL+sampled_info[K][L].num_sum/TotNumSamples;

				//TCOL is the equilibrium collision rate, see eqn (4.78)
				TCOL[M][L]=2.*TIME*Spwt*SML*SLL*(1./(cell_geom[NumCells-1].x_max-cell_geom[0].x_min))
				   *species_interaction[M][L].sigma*
						Math.pow((AVTMP/species_interaction[M][L].ref_temp),1.-species_interaction[M][L].visc_temp_index)
				   *Math.sqrt(2.*BOLTZ*species_interaction[L][M].ref_temp/(PI*species_interaction[L][M].reduced_mass));
				
			}
		}

		for (L=0;L<NumSpecies;L++)
		{
			for (M=0;M<NumSpecies;M++)
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
		StreamNumDen=1.E20;        //FND  is the number densty
		StreamTemp=300.;           //--FTMP is the temperature
		SpeciesFraction[0]=.6;         ////FSP(N) is the number fraction of species N
		SpeciesFraction[1]=.2;
		SpeciesFraction[2]=.1;
		SpeciesFraction[3]=.08;
		SpeciesFraction[4]=.02;

		Spwt=1.0E17;      //FNUM  is the number of real molecules represented by a simulated mol.
		DeltaT=.25E-4;       //--DTM is the time step
		SubCellsPerCell=8;            //--NSC is the number of sub-cells in each cell
		XMIN=0.;            ////the simulated region is from x=XF to x=XR
		XMAX=1.;
		UnlikeColsFlag=0;           ////the cross-collision data is set to the mean values

		for (int i=0;i<NumSpecies;i++) species_info[i] = new SpeciesInfo();
		for (int i=0;i<NumCells;i++)
			for (int j=0;j<NumSpecies;j++)	sampled_info[i][j] = new SampledInfo();
		
		for (int i=0;i<NumCells;i++)
			for (int j=0;j<NumSpGroups;j++)	cell_info[i][j] = new CellInfo();
			
		for (int i=0;i<NumCells;i++) cell_geom[i] = new CellGeom();
		for (int i=0;i<np_max;i++) part[i]=new Part();
		
		for (int i=0;i<NumCells;i++)
			for (int j=0;j<NumSpGroups;j++)
				for (int k=0;k<NumSpGroups;k++)
					spg_col_data[i][j][k] = new SpGroupColData();

		for (int i=0;i<NumSubCells;i++)
			for (int j=0;j<NumSpGroups;j++)
				spg_cell_info[i][j] = new SpGroupCellInfo(); //--ISCG(N,M,L) is the information on species group L in sub-cell M
		
		for (int i=0;i<NumSpGroups;i++)
			group_info[i] = new GroupInfo(); 
	
		for(int i=0;i<NumSpecies;i++)
			for(int j=0;j<NumSpecies;j++)
				species_interaction[i][j] = new SpeciesInteraction(); 
	
		//SP(1,N) is the molecular diameter of species N
		//SP(2,N) is the reference temperature
		//SP(3,N) is the viscosity-temperatire index
		//SP(4,N) is the reciprocal of the VSS scattering parameter
		//SP(5,N) is the molecular mass of species N
		//ISP(N) is the group for species N
		species_info[0].diam=3.5E-10;
		species_info[0].ref_temp=273.;
		species_info[0].visc_temp_index=0.75;
		species_info[0].vss_scat_inv=1.;
		species_info[0].mass=5.E-26;
		Species2SamplingGroup[0]=0;
		
		species_info[1].diam=4.E-10;
		species_info[1].ref_temp=273.;
		species_info[1].visc_temp_index=0.75;
		species_info[1].vss_scat_inv=1.;
		species_info[1].mass=4.5E-26;
		Species2SamplingGroup[1]=0;
		
		species_info[2].diam=3.E-10;
		species_info[2].ref_temp=273.;
		species_info[2].visc_temp_index=0.75;
		species_info[2].vss_scat_inv=1.;
		species_info[2].mass=2.5E-26;
		Species2SamplingGroup[2]=0;
		
		species_info[3].diam=3.E-10;
		species_info[3].ref_temp=273.;
		species_info[3].visc_temp_index=0.75;
		species_info[3].vss_scat_inv=1.;
		species_info[3].mass=2.E-26;
		Species2SamplingGroup[3]=0;
		
		species_info[4].diam=4.E-10;
		species_info[4].ref_temp=273.;
		species_info[4].visc_temp_index=0.75;
		species_info[4].vss_scat_inv=1.;
		species_info[4].mass=4.E-26;
		Species2SamplingGroup[4]=0;
	}
}