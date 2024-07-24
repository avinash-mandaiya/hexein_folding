#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#define BMAX 1e10
#include "mincost.h"
#include "PHS.h"

int *row,*col,*col2,*i_chain,*j_chain,**mark;
double b_min;

void setup() ;
int cov_col() ;
int cov_row() ;
void add_sub() ;
void chain() ;

void makevarspp2()
	{
	int chainlen,i,j;
	chainlen = 2*NDUNum;
	row=malloc(NDUNum*sizeof(int));	
	col=malloc(NDUNum*sizeof(int));	
	col2=malloc(NDUNum*sizeof(int));	
	i_chain=malloc(chainlen*sizeof(int));	
	j_chain=malloc(chainlen*sizeof(int));	
	mark=malloc(NDUNum * sizeof(int*));	
	for(i = 0; i<NDUNum; ++i){
		mark[i]=malloc(NDUNum * sizeof(int));
		for (j = 0 ; j < NDUNum ; ++j){
			mark[i][j] = 0;
			}
		}
	}

void freemem()
	{
	free(row);
	free(col);
	free(col2);
	free(i_chain);
	free(j_chain);
	for (int i=0;i<NDUNum;++i)
		free(mark[i]);
	free(mark);	
	}

void perm_proj( double **CostM, int *PermArray)
	{
	makevarspp2();
	int i,j,k;
	k = 0 ;

	setup() ;
	
	while ( !cov_col() )
		{
		while ( !cov_row() )
			add_sub() ;
			
		chain() ;
		}
		
	for (i = 0 ; i < NDUNum ; ++i)
		for (j = 0 ; j < NDUNum ; ++j)
			if (mark[i][j] == 2)
				{
				PermArray[i] = j ;
				break ;
				}
	freemem();
	}
	
	
void setup()
//setup finds the minimum weight in all the rows and subtract each row element
//with that weight and saves the index of those minimum CostM in each rows in 
//the matrix "mark" by switching 0 to 1.
	{
	int i,j,k;
	
	for (i = 0 ; i < NDUNum ; ++i)
		{
		b_min = CostM[i][0] ;
		k = 0 ;
		for (j = 1 ; j < NDUNum ; ++j)
			if ( CostM[i][j] < b_min )
				{
				b_min = CostM[i][j] ;
				k = j ;
				}
		
		for (j = 0 ; j < NDUNum ; ++j)
			CostM[i][j] -= b_min ;
			
		mark[i][k] = 1 ;
		}
	
	for (j = 1 ; j < NDUNum ; ++j)
		col2[j] = 0 ;
		
	for (i = 0 ; i < NDUNum ; ++i)
		{
		for (j = 0 ; j < NDUNum ; ++j)
			{
			if (col2[j])
				continue ;
				
			if (mark[i][j] == 1)
				{
				mark[i][j] = 2 ;
				col2[j] = 1 ;
				goto nexti ;
				}
			}
		nexti : ;
		}
	}
	
	
int cov_col()
	{
	int i,j,k;
	
	for (i = 0 ; i < NDUNum ; ++i)
		row[i] = 0 ;
		
	for (j = 0 ; j < NDUNum ; ++j)
		col[j] = 0 ;
		
	k = 0 ;
	for (i = 0 ; i < NDUNum ; ++i)
		for (j = 0 ; j < NDUNum ; ++j)
			{
			if (col[j])
				continue ;
				
			if (mark[i][j] == 2)
				{
				col[j] = 1 ;
				++k ;
				}
			}
			
	if (k == NDUNum)
		return 1 ;
	else
		return 0 ;
	}
	
	
int cov_row()
	{
	int i,j,k;
	
	int j2, j3, i1, j1 ;
	
	start:
	
	for (i = 0 ; i < NDUNum ; ++i)
		{
		if (row[i])
			continue ;
			
		j2 = -1 ;
		j3 = -1 ;
		for (j = 0 ; j < NDUNum ; ++j)
			{
			if (j2 == -1 && mark[i][j] == 2)
				{
				j2 = j ;
				continue ;
				}
				
			if (j3 == -1 && !col[j] && mark[i][j] > 0)
				{
				mark[i][j] = 3 ;
				j3 = j ;
				}
			}
			
		if (j3 == -1)
			continue ;
		else if (j2 == -1)
			{
			i_chain[0] = i ;
			j_chain[0] = j3 ;
			return 1 ;
			}
		else
			{
			row[i] = 1 ;
			col[j2] = 0 ;
			goto start ;
			}
		}
		
	b_min = BMAX ;
	for (i = 0 ; i < NDUNum ; ++i)
		{
		if (row[i])
			continue ;
			
		for (j = 0 ; j < NDUNum ; ++j)
			if (!col[j] && CostM[i][j] < b_min)
				{
				b_min = CostM[i][j] ;
				i1 = i ;
				j1 = j ;
				}
		}
		
	mark[i1][j1] = 1 ;
	
	return 0 ;
	}
		

void add_sub()
	{
	int i,j;
	
	for (i = 0 ; i < NDUNum ; ++i)
		if (row[i])
			{
			for (j = 0 ; j < NDUNum ; ++j)
				if (col[j])
					{
					CostM[i][j] += b_min ;
					if (mark[i][j] > 0)
						mark[i][j] = 0 ;
					}
			}
		else
			{
			for (j = 0 ; j < NDUNum ; ++j)
				if (!col[j])
					CostM[i][j] -= b_min ;
			}
				
	return ;
	}
	
	
void chain()
	{
	int i,j,k;
	int kc = 0 ;
	
	start: 
	
	for (i = 0 ; i < NDUNum ; ++i)
		if (mark[i][j_chain[kc]] == 2)
			{
			++kc ;
			j_chain[kc] = j_chain[kc-1] ;
			i_chain[kc] = i ;
			
			for (j = 0 ; j < NDUNum ; ++j)
				if (mark[i][j] == 3)
					{
					++kc ;
					i_chain[kc] = i ;
					j_chain[kc] = j ;
					
					goto start ;
					}
			}
			
	mark[i_chain[0]][j_chain[0]] = 2 ;
	for (k = 1 ; k <= kc ; )
		{
		mark[i_chain[k]][j_chain[k]] = 1 ;
		++k ;
		mark[i_chain[k]][j_chain[k]] = 2 ;
		++k ;
		}
		
	for (i = 0 ; i < NDUNum ; ++i)
	for (j = 0 ; j < NDUNum ; ++j)
		if (mark[i][j] == 3)
			mark[i][j] = 1 ;
			
	for (i = 0 ; i < NDUNum ; ++i)
		row[i] = 0 ;
		
	for (j = 0 ; j < NDUNum ; ++j)
		col[j] = 0 ;
	
	return ;
	}
	
/*	
int main(int argc,char* argv[])
	{
	int i,j,*PermArray;
	
	srand(time(NULL));

	if(argc==2)
                {
                NDUNum=atoi(argv[1]);
                }
        else
                {
                printf("expected 1 argument: NDUNum");
                return 1;
                }
	CostM = malloc(NDUNum*sizeof(double));
	PermArray = malloc(NDUNum*sizeof(int));
	for(i = 0; i<NDUNum; i++)
		CostM[i] = malloc(NDUNum*sizeof(double));
	
	for (i = 0 ; i < NDUNum ; ++i)
		{
		for (j = 0 ; j < NDUNum ; ++j){
			CostM[i][j]=rand()%100;
			printf("%lf\t",CostM[i][j]);}
		PermArray[i] = 0;		
		printf("\n");
		}
	printf("\n");
	
	perm_proj(CostM, PermArray) ;
	
	for (i = 0 ; i < NDUNum ; ++i)
		printf("%d ", PermArray[i]) ;
	printf("\n") ;
	
	return 0 ;
	}
			
*/			
