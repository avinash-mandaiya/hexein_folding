#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include "PHS.h"

int SeqLen, NumSol, LTSize, NDTNum, NDUNum, *Seq;

char *SeqID;
int iter,maxiter;
int (*GridInfo)[9], Vec[6][2], **GridLabel;
int Vec[6][2] = {{-1,1},{-1,0},{0,-1},{1,-1},{1,0},{0,1}};
int *PermArray;
int MotifNum, (*Motifs)[7];
double **CostM, **SeqPos, (*ColPos)[3], **NumVarSP, (*NumVarCP)[3];
int *TruePos;

double (**PS), (**PSA), (**PSR), (**PSB);
double (***CC), (***CCA), (***CCR), (***CCB);
double (*CG)[3], (*CGA)[3], (*CGR)[3], (*CGB)[3];
double (*LM)[7][3], (*LMA)[7][3], (*LMR)[7][3], (*LMB)[7][3];

double etaPS,*etaCC,*etaLM,etaCG;
double tPSerr,tCCerr,tLMerr,tCGerr,*CCerr,*LMerr, toterr, avgerr;
double beta,epsilon;

char errfile[50],statsfile[50],solfile[50],SeqFile[50];


static inline double sq(double diff)
        {
        return diff*diff;
        }


double urand(double a, double b)
        {                               
        double num;
        num = (b-a)*(((double)rand())/RAND_MAX)+a;
        return num;
        } 


void SetGrid()
	{
	int i, Count = 0, j,m,n,p,q,sizeM,pn,qn;

	NDTNum = 3*LTSize*(LTSize-1)+1;
	NDUNum = 3*(LTSize-2)*(LTSize-1)+1;

	GridLabel = malloc((2*LTSize-1)*sizeof(int*));
	for (i=0;i<2*LTSize-1;++i) 
		GridLabel[i] = malloc((2*LTSize-1)*sizeof(int));
	
	for (i=0;i<2*LTSize-1;++i) 
		for (j=0;j<2*LTSize-1;++j) 
			GridLabel[i][j] = -1;


	GridInfo = malloc(NDTNum*sizeof(int*[9]));

	GridInfo[0][0] = 0;
	GridInfo[0][1] = 0;
	GridInfo[0][2] = 0;

	GridLabel[GridInfo[Count][1]+LTSize-1][GridInfo[Count][2]+LTSize-1] = Count;

	Count++;

	for (i=1;i<LTSize;++i)
		{
		GridInfo[Count][0] = i;
		GridInfo[Count][1] = i;
		GridInfo[Count][2] = 0;

		GridLabel[GridInfo[Count][1]+LTSize-1][GridInfo[Count][2]+LTSize-1] = Count;

		Count++;

		for (j=0;j<6;++j)
			{
			sizeM = (j != 5) ? i : i-1;

			for (m=0;m<sizeM;++m)
				{
				GridInfo[Count][0] = i;
				GridInfo[Count][1] = GridInfo[Count-1][1]+Vec[j][0];
				GridInfo[Count][2] = GridInfo[Count-1][2]+Vec[j][1];

				GridLabel[GridInfo[Count][1]+LTSize-1][GridInfo[Count][2]+LTSize-1] = Count;

				Count++;
				}
			}
		}

	// Neighbor Information

	for (i=0; i<NDUNum; ++i)
		{
		p = GridInfo[i][1];
		q = GridInfo[i][2];

		for(j=0;j<6;++j)
			{
			n = (j+4)%6;
			pn = Vec[n][0]+p+LTSize-1; 
			qn = Vec[n][1]+q+LTSize-1; 
			GridInfo[i][3+j] = GridLabel[pn][qn];
			}
		}
	}


void makevars()
        {
        int i,j;

	PermArray = malloc(NDUNum*sizeof(int));

        PS = malloc(NDUNum*sizeof(double*));
        PSA = malloc(NDUNum*sizeof(double*));
        PSR = malloc(NDUNum*sizeof(double*));
        PSB = malloc(NDUNum*sizeof(double*));

	for (i=0;i<NDUNum;++i)
		{
		PS[i] = malloc(NDUNum*sizeof(double));
		PSA[i] = malloc(NDUNum*sizeof(double));
		PSR[i] = malloc(NDUNum*sizeof(double));
		PSB[i] = malloc(NDUNum*sizeof(double));
		}

        CG = malloc(NDUNum*sizeof(double*[3]));
        CGA = malloc(NDUNum*sizeof(double*[3]));
        CGR = malloc(NDUNum*sizeof(double*[3]));
        CGB = malloc(NDUNum*sizeof(double*[3]));

        CC = malloc(NDUNum*sizeof(double**));
        CCA = malloc(NDUNum*sizeof(double**));
        CCR = malloc(NDUNum*sizeof(double**));
        CCB = malloc(NDUNum*sizeof(double**));

	for (i=0;i<NDUNum;++i)
		{
		CC[i] = malloc(7*sizeof(double*));
		CCA[i] = malloc(7*sizeof(double*));
		CCR[i] = malloc(7*sizeof(double*));
		CCB[i] = malloc(7*sizeof(double*));

		for (j=0;j<7;++j)
			{
			CC[i][j] = malloc(NDUNum*sizeof(double));
			CCA[i][j] = malloc(NDUNum*sizeof(double));
			CCR[i][j] = malloc(NDUNum*sizeof(double));
			CCB[i][j] = malloc(NDUNum*sizeof(double));
			}
		}

        LM = malloc(NDUNum*sizeof(double*[7][3]));
        LMA = malloc(NDUNum*sizeof(double*[7][3]));
        LMR = malloc(NDUNum*sizeof(double*[7][3]));
        LMB = malloc(NDUNum*sizeof(double*[7][3]));


        CostM = malloc(NDUNum*sizeof(double*));

	for (i=0;i<NDUNum;++i)
		CostM[i] = malloc(NDUNum*sizeof(double));

       	SeqPos = malloc(NDUNum*sizeof(double*));
       	NumVarSP = malloc(NDUNum*sizeof(double*));

	for (i=0;i<NDUNum;++i)
		{
		SeqPos[i] = malloc(NDTNum*sizeof(double));
		NumVarSP[i] = malloc(NDTNum*sizeof(double));
		}

        ColPos = malloc(NDTNum*sizeof(double*[3]));
        NumVarCP = malloc(NDTNum*sizeof(double*[3]));

	etaCC = malloc(NDUNum*sizeof(double));
	etaLM = malloc(NDUNum*sizeof(double));

	CCerr = malloc(NDUNum*sizeof(double));
	LMerr = malloc(NDUNum*sizeof(double));

	TruePos = malloc(SeqLen*sizeof(int));
	}


void changeVar(double (**PSo), double (***CCo), double (*CGo)[3], double (*LMo)[7][3])
	{
	int i,l,m,n,p;
	for (l=0;l<NDUNum;++l)
		for (i=0;i<NDUNum;++i)
			PSo[i][l] = SeqPos[l][i];

	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			{
			p = (m == 0) ? i : GridInfo[i][3+m-1];

			for (l=0;l<NDUNum;++l)
				CCo[i][m][l] = SeqPos[l][p];
			}

	for (i=0;i<NDUNum;++i)
		for (m=0;m<3;++m)
			CGo[i][m] = ColPos[i][m];

	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			{
			p = (m == 0) ? i : GridInfo[i][3+m-1];

			for (n=0;n<3;++n)
				LMo[i][m][n] = ColPos[p][n];
			}
	}

			
void projA(double (**PSo), double (***CCo), double (*CGo)[3], double (*LMo)[7][3])
	{
	int i,j,m,n,l,p;
	double delta1,tempweight1;

	// Unique Matching Constraint

	for (i=0;i<NDUNum;++i)
		{
		delta1 = 0.;

		for (m=0;m<3;++m)
			CGA[i][m] = 0.;

		for (j=0;j<NDUNum;++j)
			{
			PSA[i][j] = 0.;
			delta1 += sq(PSo[i][j]);
			}

		for (j=0;j<NDUNum;++j)
			CostM[i][j] = sq(etaPS)*(delta1 - sq(PSo[i][j]) + sq(1.-PSo[i][j]));
		
		for (j=0;j<NDUNum;++j)
			{
			if (j < SeqLen)
				CostM[i][j] += sq(etaCG)*(sq(1.-CGo[i][Seq[j]]) + sq(CGo[i][1-Seq[j]]) + sq(CGo[i][2]));
			else
				CostM[i][j] += sq(etaCG)*(sq(1.-CGo[i][2]) + sq(CGo[i][1]) + sq(CGo[i][0]));
			}
		}

	perm_proj(CostM,PermArray);

	for (i=0;i<NDUNum;++i)
		{
		j = PermArray[i];
		PSA[i][j] = 1.;
		if (j < SeqLen)
			CGA[i][Seq[j]] = 1.;
		else
			CGA[i][2] = 1.;
		}

	// Sequentiality Constraint

	int minhot[7], truemin[7], minsq[7], mm, mp;
	double minhotD[7], minsqcost, delta, curmincost, MINCOST; 

	for (i=0;i<NDUNum;++i)
		{
		// Projecting all the neighbors to their nearest 1-hot vector

		for (m=1;m<7;++m)
			{
			p = GridInfo[i][2+m];
			if (GridInfo[p][0] == LTSize-1)
				{
				minhot[m] = -1;
				minhotD[m] = -FLT_MAX;
				continue;
				}

			minhotD[m] = FLT_MAX;
			for (j=0;j<NDUNum;++j)
				{
				delta = sq(1.-CCo[i][m][j]) - sq(CCo[i][m][j]);

				if (delta < minhotD[m])
					{
					minhotD[m] = delta;
					minhot[m] = j;
					}
				}
			}

		// Projecting the central site to a solvent

		curmincost = FLT_MAX;
		for (l=SeqLen;l<NDUNum;++l)
			{
			delta = sq(1.-CCo[i][0][l]) - sq(CCo[i][0][l]);

			if (delta < curmincost)
				{
				curmincost = delta;
				truemin[0] = l;
				}
			}

		MINCOST = curmincost;

		for (m=1;m<7;++m)
			truemin[m] = minhot[m];	

		// Looping over all the possible projection to residues for the central site in addition to sequentiality

		for (l=0;l<SeqLen;++l)
			{
			minsq[0] = l;  

			for (m=1;m<7;++m)
				minsq[m] = minhot[m];

			if (l == 0)
				{
				minsqcost = FLT_MAX;
				for (m=1;m<7;++m)
					{
					delta = sq(1.-CCo[i][m][l+1]) - sq(CCo[i][m][l+1]) - minhotD[m];

					if (delta < minsqcost)
						{
						minsqcost = delta;
						mp = m;
						}
					}

				minsqcost +=  sq(1.-CCo[i][0][l]) - sq(CCo[i][0][l]);

				if (minsqcost < MINCOST)
					{
					MINCOST = minsqcost;

					truemin[0] = l;

					for (m=1;m<7;++m)
						truemin[m] = minhot[m];

					truemin[mp] = l+1;
					}
				}

			else if (l == SeqLen - 1)
				{
				minsqcost = FLT_MAX;
				for (m=1;m<7;++m)
					{
					delta = sq(1.-CCo[i][m][l-1]) - sq(CCo[i][m][l-1]) - minhotD[m];

					if (delta < minsqcost)
						{
						minsqcost = delta;
						mm = m;
						}
					}

				minsqcost +=  sq(1.-CCo[i][0][l]) - sq(CCo[i][0][l]);

				if (minsqcost < MINCOST)
					{
					MINCOST = minsqcost;

					truemin[0] = l;

					for (m=1;m<7;++m)
						truemin[m] = minhot[m];	

					truemin[mm] = l-1;
					}
				}

			else
				{
				minsqcost = FLT_MAX;

				for (int m1=1;m1<7;++m1)
					for (int m2=1;m2<7;++m2)
						{
						if (m1 == m2)
							continue;

						delta = sq(1.-CCo[i][m1][l-1]) - sq(CCo[i][m1][l-1]) - minhotD[m1];
						delta += sq(1.-CCo[i][m2][l+1]) - sq(CCo[i][m2][l+1]) - minhotD[m2];

						if (delta < minsqcost)
							{
							minsqcost = delta;
							mm = m1;
							mp = m2;
							}
						}
						
				minsqcost +=  sq(1.-CCo[i][0][l]) - sq(CCo[i][0][l]);

				if (minsqcost < MINCOST)
					{
					MINCOST = minsqcost;

					truemin[0] = l;

					for (m=1;m<7;++m)
						truemin[m] = minhot[m];	

					truemin[mm] = l-1;
					truemin[mp] = l+1;
					}

				}
			}

		for (m=0;m<7;++m)
			for (j=0;j<NDUNum;++j)
				CCA[i][m][j] = 0.;

		for (m=0;m<7;++m)
			if (truemin[m] != -1)
				CCA[i][m][truemin[m]] = 1.;
		}

	// Motif Constraint

	int tempproj[7];

	for (i=0;i<NDUNum;++i)
		{
		tempweight1 = FLT_MAX;
		for (j=0;j<MotifNum;++j)
			{
			delta1 = 0;
			for (m=0;m<7;++m)
				for (n=0;n<3;++n)
					{
					if (n == Motifs[j][m])
						delta1 += sq(1-LMo[i][m][n]);
					else
						delta1 += sq(LMo[i][m][n]);
					}

			if (delta1 < tempweight1)
				{
				tempweight1 = delta1;

				for (m=0;m<7;++m)
					tempproj[m] = Motifs[j][m]; 
				}
			}

		for (m=0;m<7;++m)
			for (n=0;n<3;++n)
				{
				if (n == tempproj[m])
					LMA[i][m][n] = 1.;
				else
					LMA[i][m][n] = 0.;
				}
		}

	}


void reflect (double (**PSo), double (***CCo), double (*CGo)[3], double (*LMo)[7][3])
	{
	int i,j,m,n;

	for (i=0;i<NDUNum;++i)
		for (j=0;j<NDUNum;++j)
			PSR[i][j] = 2.*PSo[i][j] - PS[i][j];
	
	
	for (i=0;i<NDUNum;++i)
		for (j=0;j<3;++j)
			CGR[i][j] = 2.*CGo[i][j] - CG[i][j];

	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			for (j=0;j<NDUNum;++j)
				CCR[i][m][j] = 2.*CCo[i][m][j] - CC[i][m][j];

	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			for (n=0;n<3;++n)
				LMR[i][m][n] = 2.*LMo[i][m][n] - LM[i][m][n];
	}


void projB(double (**PSo), double (***CCo), double (*CGo)[3], double (*LMo)[7][3])
	{
	int i,l,m,n,p;

	for (l=0;l<NDUNum;++l)
		for (i=0;i<NDTNum;++i)
			{
			SeqPos[l][i] = 0.;
			NumVarSP[l][i] = 0.;
			}

	for (i=0;i<NDTNum;++i)
		for (m=0;m<3;++m)
			{
			ColPos[i][m] = 0.;
			NumVarCP[i][m] = 0.;
			}

	for (l=0;l<NDUNum;++l)
		for (i=0;i<NDUNum;++i)
			{
			SeqPos[l][i] += sq(etaPS)*PSo[i][l];
			NumVarSP[l][i] += sq(etaPS);
			}

	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			{
			p = (m == 0) ? i : GridInfo[i][3+m-1];

			for (l=0;l<NDUNum;++l)
				{
				SeqPos[l][p] += sq(etaCC[i])*CCo[i][m][l];
				NumVarSP[l][p] += sq(etaCC[i]);
				}
			}

	for (i=0;i<NDUNum;++i)
		for (m=0;m<3;++m)
			{
			ColPos[i][m] += sq(etaCG)*CGo[i][m];
			NumVarCP[i][m] += sq(etaCG);
			}

	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			{
			p = (m == 0) ? i : GridInfo[i][3+m-1];

			for (n=0;n<3;++n)
				{
				ColPos[p][n] += sq(etaLM[i])*LMo[i][m][n];
				NumVarCP[p][n] += sq(etaLM[i]);
				}
			}

	for (l=0;l<NDUNum;++l)
		for (i=0;i<NDTNum;++i)
			SeqPos[l][i] /= NumVarSP[l][i];

	for (i=0;i<NDUNum;++i)
		for (m=0;m<3;++m)
			ColPos[i][m] /= NumVarCP[i][m];

	for (i=NDUNum;i<NDTNum;++i)
		{
		ColPos[i][0] = 0.;
		ColPos[i][1] = 0.;
		ColPos[i][2] = 1.;

		for (l=0;l<SeqLen;++l)
			SeqPos[l][i] = 0.;
		}
/*
	int aa0 = 24;
	int aa1 = 23;

	for (i=0;i<NDTNum;++i)
		{
		SeqPos[aa0][i] = 0.;
		SeqPos[aa1][i] = 0.;
		}
	
	SeqPos[aa0][0] = 1.;
	SeqPos[aa1][1] = 1.;

	ColPos[0][0] = 0.;
	ColPos[0][1] = 1.;
	ColPos[0][2] = 0.;

	ColPos[1][0] = 0.;
	ColPos[1][1] = 1.;
	ColPos[1][2] = 0.;
*/
	changeVar(PSB,CCB,CGB,LMB);

	}

void RRR()
        {
	int i,j,l,m,n,totVar,totCons, NumVarCG,NumVarPS,NumVarLM;
        double diff,avgerr;

	projA(PS,CC,CG,LM);
	reflect(PSA,CCA,CGA,LMA);
	projB(PSR,CCR,CGR,LMR);

        tPSerr = 0.;
        tCCerr = 0.;
        tCGerr = 0.;
        tLMerr = 0.;
        toterr = 0.;
        avgerr = 0.;

	totVar = 0;
	totCons = 0;

	NumVarCG = 0;
	NumVarPS = 0;

	for (i=0;i<NDUNum;++i)
		{
		for (j=0;j<NDUNum;++j)
			{
			diff = PSB[i][j] - PSA[i][j];
			PS[i][j] += beta*diff;
			tPSerr += sq(diff);
			NumVarPS++;
			}

		for (j=0;j<3;++j)
			{
			diff = CGB[i][j] - CGA[i][j];
			CG[i][j] += beta*diff;
			tCGerr += sq(diff);
			NumVarCG++;
			}
		}

	totVar += NumVarPS;
	toterr += tPSerr;

	totVar += NumVarCG;
	toterr += tCGerr;

	tPSerr /= NumVarPS;
	tCGerr /= NumVarCG;

	avgerr += tPSerr;
	totCons++;

	avgerr += tCGerr;
	totCons++;

	// Sequentiality Constraint

	for (i=0;i<NDUNum;++i)
		{
		CCerr[i] = 0.;

		for (m=0;m<7;++m)
			for (l=0;l<NDUNum;++l)
				{
				diff = CCB[i][m][l] - CCA[i][m][l];
				CC[i][m][l] += beta*diff;
				CCerr[i] += sq(diff);
				}
		
		tCCerr += CCerr[i];
		CCerr[i] /= (NDUNum*7.);
		avgerr += CCerr[i];

		}

	totVar += NDUNum*NDUNum*7;
	toterr += tCCerr;

	totCons += NDUNum;
	tCCerr /= NDUNum*NDUNum*7;



	NumVarLM = 0;

	for (i=0;i<NDUNum;++i)
		{
		LMerr[i] = 0.;
		for (m=0;m<7;++m)
			for (n=0;n<3;++n)
				{
				diff = LMB[i][m][n] - LMA[i][m][n];
				LM[i][m][n] += beta*diff;
				LMerr[i] += sq(diff);
				NumVarLM++;
				}

		tLMerr += LMerr[i];
		LMerr[i] /= 21.;
		avgerr += LMerr[i];
		totCons++;
		}

	totVar += NumVarLM;
	toterr += tLMerr;

	tLMerr /= NumVarLM;


        toterr /= totVar;
        avgerr /= totCons;


	// weight tuning 

	etaPS += epsilon*((tPSerr/avgerr) - etaPS);

	etaCG += epsilon*((tCGerr/avgerr) - etaCG);

	for (i=0;i<NDUNum;++i)
		etaCC[i] += epsilon*((CCerr[i]/avgerr) - etaCC[i]);

	for (i=0;i<NDUNum;++i)
		etaLM[i] += epsilon*((LMerr[i]/avgerr) - etaLM[i]);

	
	tPSerr = sqrt(tPSerr);
	tCGerr = sqrt(tCGerr);
	tCCerr = sqrt(tCCerr);
	tLMerr = sqrt(tLMerr);

	toterr = sqrt(toterr);
	}


void printsol(char *solfile)
        {   
        FILE *fp;
        int i,j,k,m,p,tempproj,l;
	int ColSol[NDTNum], SeqSol[SeqLen];
	int ColState[NDTNum], SeqState[NDTNum];
	double StateErr[NDTNum];
	double StateVars[NDTNum];
	double Colerr[NDTNum];
	double WCol[NDTNum];
	double tempweight, projcost, n;

	for (i=0;i<NDTNum;++i)
		ColState[i] = 2;

        for(i=0;i<NDUNum;++i)
		{
		tempweight = FLT_MAX;
		for(j=0;j<NDUNum;++j)
			{
			projcost = sq(1 - SeqPos[j][i]) - sq(SeqPos[j][i]);

			if (projcost < tempweight)
				{
				tempweight = projcost;
				tempproj = j;
				}
			}

		SeqState[i] = tempproj;
		ColState[i] = (tempproj < SeqLen) ? Seq[tempproj] : 2;
		}

	for (i=0;i<NDTNum;++i)
		{
		StateErr[i] = 0.;
		StateVars[i] = 0.;
		}
		

	for (i=0;i<NDUNum;++i)
		{
		for (j=0;j<NDUNum;++j)
			{
			n = (SeqState[i] == j) ? 1. : 0.;
			StateErr[i] += sq(n - PSA[i][j]);
			StateVars[i] += 1.;
			}

		for (j=0;j<3;++j)
			{
			n = (ColState[i] == j) ? 1. : 0.;
			StateErr[i] += sq(n - CGA[i][j]);
			StateVars[i] += 1.;
			}
		}

	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			{
			p = (m == 0) ? i : GridInfo[i][3+m-1];
			for (l=0;l<NDUNum;++l)
				{
				if (GridInfo[p][0] == LTSize-1)
					n = 0.;
				else
					n = (SeqState[p] == l) ? 1. : 0.;

				StateErr[p] += sq(n - CCA[i][m][l]);
				StateVars[p] += 1.;
				}
			}

	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			{
			p = (m == 0) ? i : GridInfo[i][3+m-1];
			for (l=0;l<3;++l)
				{
				n = (ColState[p] == l) ? 1. : 0.;
				StateErr[p] += sq(n - LMA[i][m][l]);
				StateVars[p] += 1.;
				}
			}

	for (i=0;i<NDTNum;++i)
		{
		StateErr[i] /= StateVars[i];
		StateErr[i] = sqrt(StateErr[i]);
		}

/*
	// double Seqerr[SeqLen], WSeq[SeqLen];

        for(i=0;i<NDTNum;++i)
		{   
		tempweight = FLT_MAX;
		for(k=0;k<3;++k)
			{
			projcost = sq(1 - ColPos[i][k]) - sq(ColPos[i][k]);

			if (projcost < tempweight)
				{
				tempweight = projcost;
				tempproj = k;
				}
			}

		ColSol[i] = tempproj;
		}   

        for(i=0;i<NDUNum;++i)
		{
		tempweight = FLT_MAX;
		for(j=0;j<NDTNum;++j)
			{
			projcost = sq(1 - SeqPos[i][j]) - sq(SeqPos[i][j]);

			if (projcost < tempweight)
				{
				tempweight = projcost;
				tempproj = j;
				}
			}

		SeqSol[i] = tempproj;
		}
*/


/*
        for(i=0;i<SeqLen;++i)
		{
		Seqerr[i] = 0.;
		WSeq[i] = 0.;
		}

	for (i=0;i<NDUNum;++i)
		{
		Colerr[i] = 0.;
		WCol[i] = 0.;
		}


	for (l=0;l<SeqLen;++l)
		{
		for (i=0;i<NDUNum;++i)
			{
			if (SeqSol[l] == i)
				Seqerr[l] += sq(etaPS)*sq(1-PSB[i][l]);
			else
				Seqerr[l] += sq(etaPS)*sq(PSB[i][l]);
			WSeq[l] += sq(etaPS);
			}
		}

	for (l=1;l<SeqLen-1;++l)
		for (i=0;i<NDUNum;++i)
			{
			SeqPos[l][i] += sq(etaCC[l])*CCo[l][i][0][0];
			NumVarSP[l][i] += sq(etaCC[l]);

			for (m=0;m<6;++m)
				{
				p = GridInfo[i][3+m]; 
				SeqPos[l-1][p] += sq(etaCC[l])*CCo[l][i][1][m+1];
				NumVarSP[l-1][p] += sq(etaCC[l]);
				SeqPos[l+1][p] += sq(etaCC[l])*CCo[l][i][2][m+1];
				NumVarSP[l+1][p] += sq(etaCC[l]);
				}
			}

        for(i=0;i<NDTNum;++i)
		{
		Colerr[i] = 0.;
		WCol[i] = 0.;
		}

	for (i=0;i<NDUNum;++i)
		for (m=0;m<3;++m)
			{
			Colerr[i] += sq(etaPS)*sq(CGB[i][m]-CGA[i][m]);
			WCol[i] += sq(etaPS);
			}

	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			{
			p = (m == 0) ? i : GridInfo[i][3+m-1];

			for (n=0;n<3;++n)
				{
				Colerr[p] += sq(etaLM[i])*sq(LMB[i][m][n]-LMA[i][m][n]);
				WCol[p] += sq(etaLM[i]);
				}
			}

        for(i=0;i<NDTNum;++i)
		{
		Colerr[i] /= WCol[i];
		Colerr[i] = sqrt(Colerr[i]);
		}
*/

        fp=fopen(solfile,"a");
        for(i=0;i<NDTNum;++i)
                fprintf(fp,"%d ",ColState[i]);
        fprintf(fp,"\n");

        for(i=0;i<NDUNum;++i)
                fprintf(fp,"%d ",SeqState[i]+1);
        fprintf(fp,"\n");

        for(i=0;i<NDTNum;++i)
                fprintf(fp,"%lf ",StateErr[i]);
        fprintf(fp,"\n\n");


        fclose(fp);

        }


int solve(int maxiter,int iterstride,double stoperr)
        {
        FILE *fperr;
	FILE *fpsol;

        fperr=fopen(errfile,"w");

        for(iter=0;iter<=maxiter;++iter)
                {
                RRR();

                if(iter%iterstride==0 || toterr<stoperr)
                        {
                        fprintf(fperr,"%d \t %.6e\t%.6e %.6e %.6e %.6e\n",iter,toterr,tPSerr,tCGerr,tCCerr,tLMerr);
                        printsol(solfile);
                        }   

                if(toterr<stoperr)
                        {   
                        fclose(fperr);
                        return iter;
                        }
                }

        fclose(fperr);

        return 0;
        }


int getseq(char *seq_file)
        {
        int i;
        FILE *fp;

        fp=fopen(seq_file,"r");
        if(!fp)
                {
                printf("sequence file not found\n");
                return 0;
                }

        fscanf(fp,"%d%*[\n]",&SeqLen);

        Seq=malloc(SeqLen*sizeof(int));

	for(i=0;i<SeqLen;++i)
		fscanf(fp,"%d%*[ ]%*[\n]",&Seq[i]);

        fclose(fp);

	return 1;
	}


int getMotifs(char *motif_file)
        {
        int i,j;
        FILE *fp;

        fp=fopen(motif_file,"r");
        if(!fp)
                {
                printf("motif_file not found\n");
                return 0;
                }

        fscanf(fp,"%d%*[\n]",&MotifNum);
	
        Motifs=malloc(MotifNum*sizeof(int*[7]));

	for(i=0;i<MotifNum;++i)
		for(j=0;j<7;++j)
			fscanf(fp,"%d%*[ ]%*[\n]",&Motifs[i][j]);

        fclose(fp);

	return 1;
	}


void init()
	{
        int i,m,l;

	for (l=0;l<NDUNum;++l)
		for (i=0;i<NDTNum;++i)
			SeqPos[l][i] = urand(0,1);

	for (i=0;i<NDUNum;++i)
		for (m=0;m<3;++m)
			ColPos[i][m] = urand(0,1);

	for (i=NDUNum;i<NDTNum;++i)
		{
		ColPos[i][0] = 0.;
		ColPos[i][1] = 0.;
		ColPos[i][2] = 1.;
		}

	changeVar(PS,CC,CG,LM);

	etaPS = 1.;

	etaCG = 1.;

	for (i=0;i<NDUNum;++i)
		etaCC[i] = 1.;

	for (i=0;i<NDUNum;++i)
		etaLM[i] = 1.;
	}


void initSol(double lambda)
	{
        int i,m,l;
        FILE *fp;
        char tsolfile[50];

	snprintf(tsolfile, 50, "%s_tsol.txt", SeqID);

        fp = fopen(tsolfile,"r");
        if(!fp)
                printf("true solution file not found\n");

	for (l=0;l<NDUNum;++l)
		fscanf(fp,"%d%*[ ]%*[\n]",&TruePos[l]);

        fclose(fp);

	for (l=0;l<NDUNum;++l)
		for (i=0;i<NDTNum;++i)
			SeqPos[l][i] = 0.;

	for (i=0;i<NDTNum;++i)
		{
		ColPos[i][0] = 0.;
		ColPos[i][1] = 0.;
		ColPos[i][2] = 1.;
		}

	for (l=0;l<SeqLen;++l)
		{
		SeqPos[l][TruePos[l]] = 1.;	
		ColPos[TruePos[l]][Seq[l]] = 1.;
		ColPos[TruePos[l]][2] = 0.;
		}

	for (l=SeqLen;l<NDUNum;++l)
		SeqPos[l][TruePos[l]] = 1.;	

	for (l=0;l<NDUNum;++l)
		for (i=0;i<NDTNum;++i)
			SeqPos[l][i] += lambda*urand(-1,1);

	for (i=0;i<NDTNum;++i)
		for (m=0;m<3;++m)	
			ColPos[i][m] += lambda*urand(-1,1);


	changeVar(PS,CC,CG,LM);

	etaPS = 1.;

	etaCG = 1.;

	for (i=0;i<NDUNum;++i)
		etaCC[i] = 1.;

	for (i=0;i<NDUNum;++i)
		etaLM[i] = 1.;
	}


int main(int argc,char* argv[])
        {
        char *id,MotifFile[20];
        int c,iterstride,i,j,m,seed,numrun;
        double deltat, stoperr, iterpersec;

        FILE *fp;
        clock_t start;

        if(argc==11)
                {
                SeqID = argv[1];
                LTSize = atoi(argv[2]);
                id = argv[3];
                beta = atof(argv[4]);
                maxiter = atoi(argv[5]);
                iterstride = atoi(argv[6]);
                stoperr = atof(argv[7]);
                epsilon = atof(argv[8]);
                seed = atoi(argv[9]);
		numrun = atoi(argv[10]);
                }
        else
                {
                printf("expected more/less arguments \n");
                return 1;
                }


        snprintf(errfile,50,"%s.err",id);
        snprintf(statsfile,50,"%s.stats",id);
        snprintf(solfile,50,"%s.sol",id);

        snprintf(SeqFile,50,"%s.txt",SeqID);

        if(!getseq(SeqFile))
                return 1;

	snprintf(MotifFile,20,"PHS_motifs.txt");
        if(!getMotifs(MotifFile))
                return 1;

        fp=fopen(solfile,"w");
        fclose(fp);

        //fp=fopen(CAfile,"w");
        //fclose(fp);

        fp=fopen(statsfile,"w");
        for(c=0;c<argc;++c)
                fprintf(fp,"%s ",argv[c]);
        fprintf(fp,"\n\n");

        fprintf(fp,"Sequence Length: %d\n\n",SeqLen);
        fprintf(fp,"Sequence :\n");

        for(c=0;c<SeqLen;++c)
                fprintf(fp,"%d ",Seq[c]);
        fprintf(fp,"\n\n");

        fclose(fp);

	SetGrid();

        makevars();

        //initSol(0.0);

	//changeVar(PS,CC,CG,LM);
	//projA(PS,CC,CG,LM);
	//reflect(PSA,CCA,CGA,LMA);
	//projB(PSR,CCR,CGR,LMR);
	//RRR();
	//iter=solve(maxiter,iterstride,stoperr);

/*
	for(i=0;i<numrun;++i)
		{
		srand(i);
		init();

		start=clock();
		iter=solve(maxiter,iterstride,stoperr);

		deltat=((double)(clock()-start))/CLOCKS_PER_SEC;

		fp=fopen(statsfile,"a");
		if(iter)
			{
			iterpersec=iter/deltat;
			fprintf(fp,"%d \t %lf \t %lf\n",iter,deltat,iterpersec);
			}
		else
			{
			iterpersec=maxiter/deltat;
			fprintf(fp,"-1 \t %lf \t %lf\n",deltat,iterpersec);
			}

		fclose(fp);
		}
*/

	srand(seed);
	init();

	start=clock();
	iter=solve(maxiter,iterstride,stoperr);

	deltat=((double)(clock()-start))/CLOCKS_PER_SEC;

	fp=fopen(statsfile,"a");
	if(iter)
		{
		iterpersec=iter/deltat;
		fprintf(fp,"%d \t %lf \t %lf\n",iter,deltat,iterpersec);
		}
	else
		{
		iterpersec=maxiter/deltat;
		fprintf(fp,"-1 \t %lf \t %lf\n",deltat,iterpersec);
		}


        fprintf(fp,"\n\nColor:\n\n");

        for(i=0;i<NDUNum;++i)
                fprintf(fp,"%d \t %.2lf \t %.2lf \t %.2lf\n",i,CGB[i][0],CGB[i][1],CGB[i][2]);

        fprintf(fp,"\n\nSequence:\n\n");

        for(i=0;i<NDUNum;++i)
                fprintf(fp,"%d \t %d \t %.2lf\n",i,PermArray[i],PSB[i][PermArray[i]]);

        fprintf(fp,"\n\nConnectivity Constraint:\n\n");


	for (i=0;i<NDUNum;++i)
		for (m=0;m<7;++m)
			for (j=0;j<NDUNum;++j)
				if (CCA[i][m][j] > 0.5)
					fprintf(fp,"%d \t %d \t %d \t %.2lf \t %.2lf\n",i,m,j,CCA[i][m][j],CCB[i][m][j]);

	fprintf(fp,"\n\nColor Positions:\n\n");

        for(i=0;i<NDTNum;++i)
                fprintf(fp,"%d \t %.2lf \t %.2lf \t %.2lf\n",i,ColPos[i][0],ColPos[i][1],ColPos[i][2]);

        fprintf(fp,"\n\nSequence Positions:\n\n");

        for(i=0;i<SeqLen;++i)
		for(j=0;j<NDTNum;++j)
			if (SeqPos[i][j] > 0.5)
				fprintf(fp,"%d \t %d \t %.2lf\n",i,j,SeqPos[i][j]);

        return 0;
	}
