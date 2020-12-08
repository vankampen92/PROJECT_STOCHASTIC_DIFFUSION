/* Important: Under Construction: Check updating, please!!!!!!!!!!!!!!!!!!!  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                             David Alonso, 2010 (c)                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <MODEL.h>

void Temporal_Dynamics_Update( Community ** My_Community,
			       Parameter_Table * Table, Stochastic_Rate * Rate,
			       int Type_of_Event, int Patch)
{
  /* This function calculates the stochastic rates after the execution of a stochastic event
     in terms of the old ones, with no recalculation. This is a way to optimize the algorithm.
     Optimization is worth it to try in dynamically, very sparsely coupled systems. 
  */
  int i,n;        /* n is the old value of population at "Patch" before the change occurred */
  Community * pVil;
  double y, w, ss, ii, xx, ww, z0, z1, dz; 
  double a, b, c; 
  int No_of_CELLS;
  double Extinction_Rate;
  double Colonization_Rate;
  Parameter_Model * P = Table->P; 

  /// M_O_D_E_L___V_A_R_I_A_B_L_E_S___C_O_D_E ( Table );
  int S, I, C;
  int X, W, K;
  S = Table->S; I = Table->I;
  X = Table->X; W = Table->W; 
  C = Table->C; K = Table->K;
  
  No_of_CELLS = P->No_of_CELLS;
  pVil          = My_Community[Patch];
  
  ss = (double)pVil->n[S];
  ii = (double)pVil->n[I];
  xx = (double)pVil->n[X];
  ww = (double)pVil->n[W];

  a  = P->M_a;  b = P->M_b;  c = P->M_c;  
  
  switch(Type_of_Event)
    {
    case  0: /* 0: Death (S --> S-1) */
      ss = ss + 1.0;
      pVil->rToI[0] -= P->H_Delta;  pVil->ratePatch -= P->H_Delta;  Rate->Total_Rate -= P->H_Delta;
      pVil->rToI[1] -= P->H_Birth;  pVil->ratePatch -= P->H_Birth;  Rate->Total_Rate -= P->H_Birth;
      pVil->rToI[3] -= P->Imm;      pVil->ratePatch -= P->Imm;      Rate->Total_Rate -= P->Imm;

      w  = ww / (ss + ii);  y = ii / (ss + ii); 
      dz = a * b * w * ii / (ii + ss - 1.0);
      pVil->rToI[3] -= dz;          pVil->ratePatch -= dz;          Rate->Total_Rate -= dz;

      dz = a * c * xx * y / (ii + ss - 1.0);
      pVil->rToI[5] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;

      break;
    case  1: /* 1: Birth (S --> S+1) */ 
      ss = ss - 1.0;
      pVil->rToI[0] += P->H_Delta;  pVil->ratePatch += P->H_Delta;  Rate->Total_Rate += P->H_Delta;
      pVil->rToI[1] += P->H_Birth;  pVil->ratePatch += P->H_Birth;  Rate->Total_Rate += P->H_Birth;
      pVil->rToI[3] += P->Imm;      pVil->ratePatch += P->Imm;      Rate->Total_Rate += P->Imm;

      w  = ww / (ss + ii);  y = ii / (ss + ii);
      
      dz = a * b * w * ii / (ii + ss + 1.0);
      pVil->rToI[3] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;

      dz = a * c * xx * y / (ii + ss + 1.0);
      pVil->rToI[5] -= dz;          pVil->ratePatch -= dz;           Rate->Total_Rate -= dz;

      break;
    case  2: /* 2: Death (I --> I-1) */ 
      ii = ii + 1.0; 
      pVil->rToI[2] -= P->H_Delta;  pVil->ratePatch -= P->H_Delta;  Rate->Total_Rate -= P->H_Delta;
      pVil->rToI[1] -= P->H_Birth;  pVil->ratePatch -= P->H_Birth;  Rate->Total_Rate -= P->H_Birth;
      
      
      P->H_Recov = Queu_Function_Recovery(Table, ww, P->H_Recov_0, P->W_Recov);
      pVil->rToI[4] -= P->H_Recov;  pVil->ratePatch -= P->H_Recov;  Rate->Total_Rate -= P->H_Recov;
      
      w  = ww / (ss + ii);  y = ii / (ss + ii);
      
      dz = a * b * w * ss / (ii + ss - 1.0);
      pVil->rToI[3] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;
      
      dz = a * c * xx * ss / (ss + ii) / (ii + ss - 1.0);
      pVil->rToI[5] -= dz;          pVil->ratePatch -= dz;          Rate->Total_Rate -= dz;
      
      break;
    case  3: /* 3: Infection (S --> S-1 and I --> I+1) */ 
      ss = ss + 1.0; ii = ii - 1.0; 
      pVil->rToI[0] -= P->H_Delta;  pVil->ratePatch -= P->H_Delta;  Rate->Total_Rate -= P->H_Delta;
      pVil->rToI[2] += P->H_Delta;  pVil->ratePatch += P->H_Delta;  Rate->Total_Rate += P->H_Delta;

      P->H_Recov = Queu_Function_Recovery(Table, ww, P->H_Recov_0, P->W_Recov);
      pVil->rToI[4] += P->H_Recov;  pVil->ratePatch += P->H_Recov;  Rate->Total_Rate += P->H_Recov;
      
      w  = ww / (ss + ii);  y = ii / (ss + ii);
      dz = a * b * w + P->Imm;
      pVil->rToI[3] -= dz;          pVil->ratePatch -= dz;          Rate->Total_Rate -= dz;
      
      dz = a * c * xx / (ii + ss);  
      pVil->rToI[5] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;
      
      break;
    case  4: /* 4: Recovery (S --> S+1  and I --> I-1) */ 
      ss = ss - 1.0; ii = ii + 1.0;
      pVil->rToI[2] -= P->H_Delta;  pVil->ratePatch -= P->H_Delta;  Rate->Total_Rate -= P->H_Delta;
      pVil->rToI[0] += P->H_Delta;  pVil->ratePatch += P->H_Delta;  Rate->Total_Rate += P->H_Delta;

      P->H_Recov = Queu_Function_Recovery(Table, ww, P->H_Recov_0, P->W_Recov);
      pVil->rToI[4] -= P->H_Recov;  pVil->ratePatch -= P->H_Recov;  Rate->Total_Rate -= P->H_Recov;

      w  = ww / (ss + ii);  y = ii / (ss + ii);
      dz = a * b * w + P->Imm;
      pVil->rToI[3] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;

      dz = a * c * xx / (ii + ss);  
      pVil->rToI[5] -= dz;          pVil->ratePatch -= dz;          Rate->Total_Rate -= dz;
      
      break;
    case  5: /* 5: Mosquito Infection (X --> X-1 and W ---> W+1) */ 
      xx = xx + 1.0;   ww = ww - 1.0;
      pVil->rToI[7] -= P->M_Delta;  pVil->ratePatch -= P->M_Delta;  Rate->Total_Rate -= P->M_Delta;
      pVil->rToI[8] += P->M_Delta;  pVil->ratePatch += P->M_Delta;  Rate->Total_Rate += P->M_Delta;
      
      z1 = Queu_Function_Recovery(Table, ww+1.0, P->H_Recov_0, P->W_Recov);
      z0 = Queu_Function_Recovery(Table, ww, P->H_Recov_0, P->W_Recov);
      dz = (z1 - z0)*ii;  
      pVil->rToI[4] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;
      
      y = ii / (ss + ii);
      dz = a * b * ss / (ss + ii); 
      pVil->rToI[3] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;
      
      dz = a * c * y;
      pVil->rToI[5] -= dz;          pVil->ratePatch -= dz;          Rate->Total_Rate -= dz;

      break; 
    case  6: /* 6: Mosquito Birth (X --> X+1) */  /* M_Delta = M_Fecundity */ 
      xx = xx - 1.0;
      pVil->rToI[7] += P->M_Delta;  pVil->ratePatch += P->M_Delta;  Rate->Total_Rate += P->M_Delta;
      pVil->rToI[6] += P->M_Delta;  pVil->ratePatch += P->M_Delta;  Rate->Total_Rate += P->M_Delta;

      y = ii / (ss + ii);
      dz = a * c * y;
      pVil->rToI[5] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;
      
      break;
    case  7: /* 7: Healthy Mosquito Death (X --> X-1) */
      xx = xx + 1.0;
      pVil->rToI[7] -= P->M_Delta;  pVil->ratePatch -= P->M_Delta;  Rate->Total_Rate -= P->M_Delta;
      pVil->rToI[6] -= P->M_Delta;  pVil->ratePatch -= P->M_Delta;  Rate->Total_Rate -= P->M_Delta;

      y = ii / (ss + ii);
      dz = a * c * y;
      pVil->rToI[5] -= dz;          pVil->ratePatch -= dz;          Rate->Total_Rate -= dz;

      break;
    case  8: /* 8: Infectious Mosquito Death (W --> W-1) */  
      ww = ww + 1.0;
      pVil->rToI[6] -= P->M_Delta;  pVil->ratePatch -= P->M_Delta;  Rate->Total_Rate -= P->M_Delta;
      pVil->rToI[8] -= P->M_Delta;  pVil->ratePatch -= P->M_Delta;  Rate->Total_Rate -= P->M_Delta;

      z1 = Queu_Function_Recovery(Table, ww-1.0, P->H_Recov_0, P->W_Recov);
      z0 = Queu_Function_Recovery(Table, ww, P->H_Recov_0, P->W_Recov);
      dz = (z1 - z0)*ii;  
      pVil->rToI[4] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;
      
      dz = a * b * ss / (ss + ii); 
      pVil->rToI[3] += dz;          pVil->ratePatch += dz;          Rate->Total_Rate += dz;
      
      break;
  
    default:
    /* Something is very very wrong!!! */
      printf("The number of the event to occur should be between 0 and 8\n");
      printf("Event to Occur = %d\n", Type_of_Event);
      Press_Key();
    }

  if(Rate->Total_Rate <= 0.0){
      printf("\n");
      printf(" R is the total temporal rate of system configuration change\n");
      printf(" R = %g\n", Rate->Total_Rate );
      printf(" As R is zero or negative, no change is possible\n");
      printf(" R shouldn't be negative. If it is, there are definitely some errors in the code\n");
      printf("\n");
      if( Rate->Total_Rate < 0.0) exit(0);
  }
}
