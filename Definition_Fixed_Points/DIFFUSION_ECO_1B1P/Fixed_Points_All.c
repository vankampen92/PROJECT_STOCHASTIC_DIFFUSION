#include <MODEL.h>

void Plasmid_Free_Solution(Parameter_Table * Table, 
                            double * Vector_Stationarity_Lower,
                            double * Vector_Stationarity_Inter,
                            double * Vector_Stationarity_Upper);
                            
void Fixed_Points_All( Parameter_Table * Table,
		                   double * Vector_Stationarity_Lower,
		                   double * Vector_Stationarity_Inter,
		                   double * Vector_Stationarity_Upper,
		                   double Epsilon)
{

  if (Table->No_of_CELLS == 1)
    Fixed_Points_OneCell(Table, Vector_Stationarity_Inter, 
                                Vector_Stationarity_Lower, 
                                Vector_Stationarity_Upper);
  else{
    printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
    printf("Thefore, the dynamics occurs in one single patch or local population\n");
    printf("However, No of CELLS is larger than 1 (N = %d)\n", Table->No_of_CELLS);
    printf("The program will safely exit\n");
    Print_Press_Key(1,0,".");
    exit(0);
  } 
}

void Fixed_Points_OneCell(Parameter_Table * Table, 
                          double * Vector_Stationarity_Lower,
		                      double * Vector_Stationarity_Inter,
		                      double * Vector_Stationarity_Upper )
  
{ 
  int i;
  double x, x_0, x_0_MIN, x_0_MAX, x_1, x_1_MIN, x_1_MAX, x_p, x_n, x_1_p, x_1_n; 
  double K_R = Table->K_R;

  if(Table->No_of_CELLS > 1) {
      printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
      printf("Thefore, the dynamics should occurs in one single patch or local population\n");
      printf("However, No of CELLS has been set to be larger than 1 (M = %d)\n", 
              Table->No_of_CELLS);
      printf("Press Key...\n");
      Print_Press_Key(1,1,"The program will safely exit\n");
      exit(0);
  }

  double D_0, D_1, beta_0, beta_1, epsilon, chi, p_C, Dis;
  double a, b, c;

  /* Relative rates to the encounter rate, Gamma (coded as Table->Lambda_R_1) */
  D_0 = Table->Delta_AP[0] / Table->Lambda_R_1; /* D_0 = (Delta_R_0 + Delta_R_1)/Gamma */
  D_1 = Table->Delta_AP[1] / Table->Lambda_R_1; /* D_1 = (Delta_R_0 + Delta_R_1 * (1.0 - Table->Nu_C_0))/Gamma */

  beta_0 = Table->Beta_AP[0] / Table->Lambda_R_1;
  beta_1 = Table->Beta_AP[1] / Table->Lambda_R_1;

  epsilon = Table->p_1;
  chi = Table->Chi_C_0;
  p_C = Table->p_2;

  // Calculate 'a'
  a = (1 - D_1 / (beta_1 * (1 - epsilon))) * 
    (epsilon / (1 - epsilon) * D_1 + chi + p_C) - 
    beta_0 * D_1 / (beta_1 * (1 - epsilon));

  // Calculate 'b'
  b = -(1 + (p_C - chi) / (beta_1 * (1 - epsilon))) * 
    (2 * epsilon / (1 - epsilon) * D_1 + p_C + chi) - 
    D_0 - beta_0 * (p_C - chi) / (beta_1 * (1 - epsilon));

  // Calculate 'c'
  c = (epsilon / (1 - epsilon) * D_1) * 
    (1 - D_1 / (beta_1 * (1 - epsilon)));

  Dis = b * b - 4 * a * c;
 
  if(Dis < 0) {
    printf("The discriminant is negative\n");
    printf("Mixed solution is not possbile\n");
    
    if (beta_0 < D_0){
      /* The system goes to full collapse: extinction */
      for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
        Vector_Stationarity_Inter[i] =                  0.0;
        Vector_Stationarity_Lower[i] =                  0.0;
        Vector_Stationarity_Upper[i] =                  0.0;
        Table->Vector_Model_Variables_Stationarity[i] = 0.0;
      }

      // Type_of_Coexistence = 0; /* Extinction */
    }
    else {
      /* One population: only Plasmid Free */
      Plasmid_Free_Solution(Table, Vector_Stationarity_Lower,
                                   Vector_Stationarity_Inter,
                                   Vector_Stationarity_Upper);

      // Type_of_Coexistence = 1; /* One solution: plasmid-free equilibrium */                             
    }
  }
  else {
    /* The discriminant is positive */
    /* Possibility of coexistence of two solutions: 
       mixed solution and plasmid-free equilibrium 
    */
    
    /* Two populations: Plasmid Free and Plasmid Carrying:
       mixed equilibrium 
    */
    x_p = (-b + sqrt(Dis)) / (2.0 * a);
    x_n = (-b - sqrt(Dis)) / (2.0 * a);
    
    if (x_p > 0.0 || x_n > 0.0) { 
      if (x_p *  x_n > 0.0) {
        /* Possibility of Coexistence of Two Stationary Solutions */
        x_0 = x_n;
        // Translate the given equation into C code
        x_1_n = (beta_1 * (1 - epsilon) - D_1) / (beta_1 * (1 - epsilon)) - 
              (beta_1 * (1 - epsilon) + p_C - chi) / (beta_1 * (1 - epsilon)) * x_0;

        x_0 = x_p;
        // Translate the given equation into C code
        x_1_p = (beta_1 * (1 - epsilon) - D_1) / (beta_1 * (1 - epsilon)) - 
                    (beta_1 * (1 - epsilon) + p_C - chi) / (beta_1 * (1 - epsilon)) * x_0;         
                    
        x_1_MIN = MIN(x_1_n, x_1_p);
        x_1_MAX = MAX(x_1_n, x_1_p);

        x_0_MIN = MIN(x_n, x_p);
        x_0_MAX = MAX(x_n, x_p);

        if (x_1_n * x_1_p > 0.0) {
          /* Two coexisting solutions: mixed populations */      
          Vector_Stationarity_Inter[0] =                 K_R * x_0_MIN;
          Vector_Stationarity_Lower[0] =                 K_R * x_0_MIN;
          Vector_Stationarity_Upper[0] =                 K_R * x_0_MAX;
          Table->Vector_Model_Variables_Stationarity[0] = K_R * x_0_MIN;

          Vector_Stationarity_Inter[1] =                  K_R * x_1_MIN;
          Vector_Stationarity_Lower[1] =                  K_R * x_1_MIN;
          Vector_Stationarity_Upper[1] =                  K_R * x_1_MAX;
          Table->Vector_Model_Variables_Stationarity[1] = K_R * x_1_MIN;

          // Type_of_Coexistence = 3; /* Two coexistence solutions */
        }
        else {
          /* One single solution: mixed populations */
          if( x_1_n > 0.0) {
            x_0 = x_n;
            x_1 = x_1_n;
          }
          else {
            x_0 = x_p;
            x_1 = x_1_p;
          }
          Vector_Stationarity_Inter[0] =                  K_R * x_0;
          Vector_Stationarity_Lower[0] =                  K_R * x_0;
          Vector_Stationarity_Upper[0] =                  K_R * x_0;
          Table->Vector_Model_Variables_Stationarity[0] = K_R * x_0;

          Vector_Stationarity_Inter[1] =                  K_R * x_1;
          Vector_Stationarity_Lower[1] =                  K_R * x_1;
          Vector_Stationarity_Upper[1] =                  K_R * x_1;
          Table->Vector_Model_Variables_Stationarity[1] = K_R * x_1;
          
          // fType_of_Coexistence = 2; /* One solution: mixed populations */
        }   
      }
      else {
        assert( x_p * x_n < 0.0); 
        x_0_MIN = MIN(x_n, x_p);
        x_0_MAX = MAX(x_n, x_p);

        x_0 = x_0_MAX;
        x_1 = (beta_1 * (1 - epsilon) - D_1) / (beta_1 * (1 - epsilon)) - 
                (beta_1 * (1 - epsilon) + p_C - chi) / (beta_1 * (1 - epsilon)) * x_0;

        if (x_1 > 0) {
            /* One single solution: mixed populations */
          Vector_Stationarity_Inter[0] =                  K_R * x_0;
          Vector_Stationarity_Lower[0] =                  K_R * x_0;
          Vector_Stationarity_Upper[0] =                  K_R * x_0;
          Table->Vector_Model_Variables_Stationarity[0] = K_R * x_0;

          Vector_Stationarity_Inter[1] =                  K_R * x_1;
          Vector_Stationarity_Lower[1] =                  K_R * x_1;
          Vector_Stationarity_Upper[1] =                  K_R * x_1;
          Table->Vector_Model_Variables_Stationarity[1] = K_R * x_1;

          // Type_of_Coexistence = 2; /* One solution: mixed populations */
        }
        else {
          /* One population: only Plasmid Free */
          Plasmid_Free_Solution(Table, Vector_Stationarity_Lower,
                                       Vector_Stationarity_Inter,
                                       Vector_Stationarity_Upper);

          // Type_of_Coexistence = 1; /* One single population: plasmid free */                             
        }      
      }
    }
    else { /* x_n and x_p are both negative */
      /* One population: only Plasmid Free */
      Plasmid_Free_Solution(Table, Vector_Stationarity_Lower,
                                   Vector_Stationarity_Inter,
                                   Vector_Stationarity_Upper);
                                  
      // Type_of_Coexistence = 1; /* One single population: plasmid free */                             
      
    }
  }
}

void Plasmid_Free_Solution(Parameter_Table * Table, 
                           double * Vector_Stationarity_Lower,
                           double * Vector_Stationarity_Inter,
                           double * Vector_Stationarity_Upper)
{
  double x;
  
  x = 1.0 - Table->Delta_AP[0] / Table->Beta_AP[0];
  double K_R = (double)Table->K_R;

  assert(x > 0.0); 

  Vector_Stationarity_Inter[0] =                  K_R * x;
  Vector_Stationarity_Lower[0] =                  K_R * x;
  Vector_Stationarity_Upper[0] =                  K_R * x;
  Table->Vector_Model_Variables_Stationarity[0] = K_R * x;

  Vector_Stationarity_Inter[1] =                  0.0;
  Vector_Stationarity_Lower[1] =                  0.0;
  Vector_Stationarity_Upper[1] =                  0.0;
  Table->Vector_Model_Variables_Stationarity[1] = 0.0;
  
}  

double Calculate_Type_of_Coexistence(Parameter_Table * Table )
{ 
  int i;
  double x, x_0, x_0_MIN, x_0_MAX, x_1, x_1_MIN, x_1_MAX, x_p, x_n, x_1_n, x_1_p; 
  double Type_of_Coexistence_Double;
  int Type_of_Coexistence;

  if(Table->No_of_CELLS > 1) {
  printf("Here, Type of Coexistence can only be assessed analytically if No_of_CELLS is 1\n");
  printf("Thefore, the dynamics should occur in one single patch or local population\n");
  printf("However, No of CELLS has been set to be larger than 1 (M = %d)\n", 
       Table->No_of_CELLS);
  printf("Press Key...\n");
  Print_Press_Key(1,1,"The program will safely exit\n");
  exit(0);
  }

  double D_0, D_1, beta_0, beta_1, epsilon, chi, p_C, Dis;
  double a, b, c;

  /* Relative rates to the encounter rate, Gamma (coded as Table->Lambda_R_1) */
  D_0 = Table->Delta_AP[0] / Table->Lambda_R_1; /* D_0 = (Delta_R_0 + Delta_R_1)/Gamma */
  D_1 = Table->Delta_AP[1] / Table->Lambda_R_1; /* D_1 = (Delta_R_0 + Delta_R_1 * (1.0 - Table->Nu_C_0))/Gamma */

  beta_0 = Table->Beta_AP[0] / Table->Lambda_R_1;
  beta_1 = Table->Beta_AP[1] / Table->Lambda_R_1;

  epsilon = Table->p_1;
  chi = Table->Chi_C_0;
  p_C = Table->p_2;

  // Calculate 'a'
  a = (1 - D_1 / (beta_1 * (1 - epsilon))) * 
    (epsilon / (1 - epsilon) * D_1 + chi + p_C) - 
    beta_0 * D_1 / (beta_1 * (1 - epsilon));

  // Calculate 'b'
  b = -(1 + (p_C - chi) / (beta_1 * (1 - epsilon))) * 
    (2 * epsilon / (1 - epsilon) * D_1 + p_C + chi) - 
    D_0 - beta_0 * (p_C - chi) / (beta_1 * (1 - epsilon));

  // Calculate 'c'
  c = (epsilon / (1 - epsilon) * D_1) * 
    (1 - D_1 / (beta_1 * (1 - epsilon)));

  Dis = b * b - 4 * a * c;

  if(Dis < 0) {
    printf("The discriminant is negative\n");
    printf("Mixed solution is not possbile\n");

    if (beta_0 < D_0){
      /* The system goes to full collapse: extinction */

      Type_of_Coexistence = 0;
    }
    else {
      /* One population: only Plasmid Free */
      
      Type_of_Coexistence = 1;                             
    }
  }
  else {
    /* The discriminant is positive */
    /* Possibility of coexistence of two solutions: 
      mixed solution and plasmid-free equilibrium 
    */
  
    /* Two populations: Plasmid Free and Plasmid Carrying:
      mixed equilibrium 
    */
    x_p = (-b + sqrt(Dis)) / (2.0 * a);
    x_n = (-b - sqrt(Dis)) / (2.0 * a);
  
    if (x_p > 0.0 || x_n > 0.0) { 
      if (x_p *  x_n > 0.0) {
      /* Possibility of Coexistence of Two Stationary Solutions */
      x_0 = x_n;
      // Translate the given equation into C code
      x_1_n = (beta_1 * (1 - epsilon) - D_1) / (beta_1 * (1 - epsilon)) - 
        (beta_1 * (1 - epsilon) + p_C - chi) / (beta_1 * (1 - epsilon)) * x_0;

      x_0 = x_p;
      // Translate the given equation into C code
      x_1_p = (beta_1 * (1 - epsilon) - D_1) / (beta_1 * (1 - epsilon)) - 
        (beta_1 * (1 - epsilon) + p_C - chi) / (beta_1 * (1 - epsilon)) * x_0;         
          
      x_1_MIN = MIN(x_1_n, x_1_p);
      x_1_MAX = MAX(x_1_n, x_1_p);

      x_0_MIN = MIN(x_n, x_p);
      x_0_MAX = MAX(x_n, x_p);

        if (x_1_n * x_1_p > 0.0) {
        /* Two coexisting solutions: mixed populations */      
       
          Type_of_Coexistence = 3; /* Two coexistence solutions */
        }
        else {
          /* One single solution: mixed populations */
          if( x_1_n > 0.0) {
            x_0 = x_n;
            x_1 = x_1_n;
          }
          else {
            x_0 = x_p;
            x_1 = x_1_p;
          }
        
          Type_of_Coexistence = 2; /* One solution: mixed populations */
        }   
      }
      else {
        assert( x_p * x_n < 0.0); 
        x_0_MIN = MIN(x_n, x_p);
        x_0_MAX = MAX(x_n, x_p);

        x_0 = x_0_MAX;
        x_1 = (beta_1 * (1 - epsilon) - D_1) / (beta_1 * (1 - epsilon)) - 
        (beta_1 * (1 - epsilon) + p_C - chi) / (beta_1 * (1 - epsilon)) * x_0;

        if (x_1 > 0) {
          /* One single solution: mixed populations */

          Type_of_Coexistence = 2; /* One solution: mixed populations */
        }  
        else {
          /* One population: only Plasmid Free */
        
          Type_of_Coexistence = 1; /* One single population: plasmid free */                             
        }       
      }
    }
    else { /* x_n and x_p are both negative */
    /* One population: only Plasmid Free */
                  
      Type_of_Coexistence = 1; /* One single population: plasmid free */                             
    }
  }

  Type_of_Coexistence_Double = (double)Type_of_Coexistence;
  return(Type_of_Coexistence_Double);
}

double Calculate_Stability_Stationary_Point( Parameter_Table * Table)
{
  double Type_of_Stability;

  Type_of_Stability = Calculate_Type_of_Stability_2D( Table );

  return(Type_of_Stability);
}