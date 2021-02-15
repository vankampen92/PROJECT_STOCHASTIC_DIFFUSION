#include <MODEL.h>


double In_Mu(Parameter_Table * Table, int m, int Sp, int J, const double * y)
{
  /* Input: 
     ------
     . m: Equation index
     . Sp: Species 
     . J: Index of the population receiveing immigrants. 
         (K: Index of the population from which immigrants come)
     . y: State vector 
     
     Output:
     -------
     . Mu: Total number of immigrations arriving to local popualtion J from all other 
     neighboring populations per unit time
  */
  /* Warning: This assumes the Patch System has been previously set up */
  Community ** P = Table->Patch_System;
  
  int K;
  double Mu;
  int i,j, k, n;

  /*
     J: Index for the J-th population 
     K: Index for the K-th population 
  */
  
  k = m%Table->LOCAL_STATE_VARIABLES;
  assert(k == Sp); 

  Mu = 0.0;

  if(Table->TYPE_of_NETWORK == 1) 
    assert( Table->No_of_NEIGHBORS == P[J]->No_NEI ); 
              
  for( j=0; j<P[J]->No_NEI; j++) {
    
    K = P[J]->Patch_Connections[j];
    n   = k + Table->LOCAL_STATE_VARIABLES*K;

    if (Table->TYPE_of_MODEL == 0 || Table->TYPE_of_MODEL == 1 )
      assert(Table->Mu == P[j]->In_Migration_Vector[Sp][j]);

    if (Sp == 0)
      assert(Table->Mu   == P[j]->In_Migration_Vector[Sp][j]);
    
    if (Sp == 1 && Table->TYPE_of_MODEL == 2)
      assert(Table->Mu_C == P[j]->In_Migration_Vector[Sp][j]);

    Mu += P[J]->In_Migration_Vector[Sp][j] * y[n]; 

    /* Mu += Table->Metapop_Connectivity_Matrix[Sp][J][K] * y[n]; */
      
  }
  
  return (Mu); 
}

double Out_Mu_Per_Capita(Parameter_Table * Table, int Sp, int J)
{
  
  /* Input: 
     ------
     . Sp: Species or Type
     . J: Index of the population exporting individuals over all its local neighboring 
     populations. 
     
     Output:
     -------
     . Mu: Per Capita emigration rate from the J-th population over all other local population 
  */
  double Mu;
  int j, K; 
  Community ** P = Table->Patch_System;

  if(Table->TYPE_of_NETWORK == 1) 
    assert( Table->No_of_NEIGHBORS == P[J]->No_NEI ); 
  
  Mu = 0.0;
  for( j=0; j < P[J]->No_NEI; j++) {

    K = P[J]->Patch_Connections[j];
    
    Mu += P[J]->Out_Migration_Vector[Sp][j]; 
    
    /* if(Table->TYPE_of_NETWORK == 0)                       */
    /*   Mu += Table->Metapop_Connectivity_Matrix[Sp][K][J]; */
    /* else                                                  */
    /*   Mu += Table->Metapop_Connectivity_Matrix[Sp][J][K]; */
    
  }

  if (Table->TYPE_of_MODEL == 0 || Table->TYPE_of_MODEL == 1) 
    assert( (P[J]->No_NEI * Table->Mu) == Mu );

  if (Table->TYPE_of_MODEL == 2) { 
    if (Sp == 0)      assert( (P[J]->No_NEI * Table->Mu)   == Mu );
    else if (Sp == 1) assert( (P[J]->No_NEI * Table->Mu_C) == Mu );
    else if (Sp == 2) assert( (P[J]->No_NEI * 0.0 * Table->Mu_C) == Mu );
    else              assert(                          0.0 == Mu ); 
  }
  
  return (Mu);
}
