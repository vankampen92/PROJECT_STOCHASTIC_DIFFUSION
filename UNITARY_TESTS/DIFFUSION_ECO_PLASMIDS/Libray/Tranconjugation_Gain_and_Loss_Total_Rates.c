#include <MODEL.h>

void Transconjugation_Gain_and_Loss_Total_Rates( const double * Y, double * Gain, double * Loss, 
                                                 Parameter_Table * Table )
{
  /* Input: 
        . Table
        . Y, state vector 
     Output
        . Gain[k]: production rate of transconjugants of Strain ID 'k'
        . Loss[k]: loss rate of recipients as they turn into transconjugants. 
  */
  int k, l, n, N;
  double G;
  int Strain_ID_R, Strain_ID_D;     
  double K_R = (double)Table->K_R; 
  
  double * L = (double *)calloc(Table->No_of_RESOURCES, sizeof(double));   
  
  for(k = 0; k < Table->No_of_RESOURCES; k++) { 

    G = 0.0; 
    for(l = 0; l < Table->Putative_Recipient_List_Indeces[k][Table->No_of_PROFILES]; l++) {
        Strain_ID_R = Table->Putative_Recipient_List_Indeces[k][l];
        
        N = Table->DoRe[k][l][0]->N; 
        for(n = 0; n < N; n++) {
           
          Strain_ID_D     = Table->DoRe[k][l][n]->ne;    
          
          G              += Table->Eta_RP[Strain_ID_D] * Table->DoRe[k][l][n]->ra * Y[Strain_ID_R]/K_R * Y[Strain_ID_D];
          L[Strain_ID_R] += Table->Eta_RP[Strain_ID_D] * Table->DoRe[k][l][n]->ra * Y[Strain_ID_R]/K_R * Y[Strain_ID_D];
        }   
    }
    Gain[k] = G;
  }

  for(k = 0; k < Table->No_of_RESOURCES; k++) 
    Loss[k] = L[k];   
  
  free(L); 
}
