#include <MODEL.h>

void Calculate_Bistable_Branches(double , double , double , double , double , 
                                double * , double * , double * , Parameter_Table * );

void Fixed_Points_All( Parameter_Table * Table,
		                   double * Vector_Stationarity_Lower,
		                   double * Vector_Stationarity_Inter,
		                   double * Vector_Stationarity_Upper,
		                   double Epsilon)
{
  int k; 
  
  if (Table->No_of_CELLS == 1) {    
  
    
    Fixed_Points_OneCell_Bis(Table, 
                                   Vector_Stationarity_Inter, 
                                   Vector_Stationarity_Lower, 
                                   Vector_Stationarity_Upper);

    printf("Lower Fixed Point (model variables):\t");
      for (k=0; k < Table->MODEL_STATE_VARIABLES; k++)
        printf("y_LOWER[%d] = %g\t", k, Vector_Stationarity_Lower[k]);
    printf("\n");
    printf("Inter Fixed Point (model variables):\t");
      for (k=0; k < Table->MODEL_STATE_VARIABLES; k++)
        printf("y_INTER[%d] = %g\t", k, Vector_Stationarity_Inter[k]);
    printf("\n");
    printf("Upper Fixed Point (model variables):\t");
      for (k=0; k < Table->MODEL_STATE_VARIABLES; k++)
        printf("y_UPPER[%d] = %g\t", k, Vector_Stationarity_Upper[k]);
    printf("\n");
      
    Stationary_Solution_dvdt_Check (Table);

    getchar();
  }
  else {
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
  double Type_of_Coexistence; 

  if(Table->No_of_CELLS > 1) {
    printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
    printf("Thefore, the dynamics should occurs in one single patch or local population\n");
    printf("However, No of CELLS has been set to be larger than 1 (M = %d)\n", 
            Table->No_of_CELLS);
    printf("Press Key...\n");
    Print_Press_Key(1,1,"The program will safely exit\n");
    exit(0);
  }

  Type_of_Coexistence = Calculate_Type_of_Coexistence(Table);  
}

double Calculate_Type_of_Coexistence( Parameter_Table * Table)
{ 
  int i;
  double x, x_0, x_p, x_n, x_1_p, x_1_n, x_0_p, x_0_n; 
  double K_R = Table->K_R;
  double Type_of_Coexistence;

  double * Vector_Stationarity_Lower =  Table->Vector_Model_Variables_MultiStability[0]; 
	double * Vector_Stationarity_Inter =  Table->Vector_Model_Variables_MultiStability[1];
	double * Vector_Stationarity_Upper =  Table->Vector_Model_Variables_MultiStability[2];

  if(Table->No_of_CELLS > 1) {
    printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
    printf("Thefore, the dynamics should occurs in one single patch or local population\n");
    printf("However, No of CELLS has been set to be larger than 1 (M = %d)\n", 
           Table->No_of_CELLS);
    printf("Press Key...\n");
    Print_Press_Key(1,1,"The program will safely exit\n");
    exit(0);
  }

  double D_0, D_1, C_0, C_1, B_0, B_1, beta_0, beta_1, epsilon, chi, p_C, Dis;
  double a, b, c;

  /* Relative rates to the encounter rate, Gamma (coded as Table->Lambda_R_1) */
  D_0 = (Table->Delta_R_0 + Table->Delta_R_1) / Table->Lambda_R_1; /* D_0 = (Delta_R_0 + Delta_R_1)/Gamma */
  D_1 = (Table->Delta_R_0 + Table->Delta_R_1*(1.0 - Table->Nu_C_0)) / Table->Lambda_R_1; /* D_1 = (Delta_R_0 + Delta_R_1 * (1.0 - Table->Nu_C_0))/Gamma, where Nu_C_0 is the resistence */

  beta_0 = Table->Beta_R / Table->Lambda_R_1;
  beta_1 = Table->Beta_R * (1.0 - Table->Alpha_C_0) / Table->Lambda_R_1;

  epsilon = Table->p_1;
  chi = Table->Chi_C_0;
  p_C = Table->p_2;

  B_0 = beta_0 - D_0;
  B_1 = beta_1 * (1.0 - epsilon) - D_1;
  C_0 = beta_0 + beta_1*epsilon + p_C + chi;
  C_1 = beta_1 * ( 1.0 - epsilon) + p_C - chi;

  // Calculate 'a'
  a = -beta_0 + (C_0 - epsilon / (1.0 - epsilon)) * (1.0 + (p_C - chi) / (beta_1 * (1.0 - epsilon)));

  // Calculate 'b'
  b = B_0 - (epsilon / (1.0 - epsilon)) * C_1 * (1.0 + 2.0 * (1.0 - D_1 / (beta_1*(1.0 - epsilon))) ) - C_0 * (1.0 - D_1 / (beta_1*(1.0 - epsilon)));
  
  // Calculate 'c'
  c = (D_1 * epsilon / (1.0 - epsilon)) * (1.0 - D_1 / (beta_1 * (1.0 - epsilon)));

  Dis = b * b - 4.0 * a * c;

  x = 1.0 - D_0 / beta_0;
  
  if(Dis < 0) {
    printf("The discriminant is negative\n");
    printf("Mixed solution is not possible\n");
    
    Type_of_Coexistence = Calculate_Plasmid_Free_Solution(beta_0, D_0, 
                                                          Vector_Stationarity_Lower,
                                                          Vector_Stationarity_Inter,
                                                          Vector_Stationarity_Upper, Table); 
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

    x_0_n = x_n;
    x_1_n =  (beta_1 * (1 - epsilon) - D_1) / (beta_1 * (1 - epsilon)) - 
           (beta_1 * (1 - epsilon) + p_C - chi) / (beta_1 * (1 - epsilon)) * x_0_n; 
   
    x_0_p = x_p;
    x_1_p = (beta_1 * (1 - epsilon) - D_1) / (beta_1 * (1 - epsilon)) - 
           (beta_1 * (1 - epsilon) + p_C - chi) / (beta_1 * (1 - epsilon)) * x_0_p;

    if ( (x_0_n >= 0.0 && x_0_n <= 1.0) && 
         (x_0_p >= 0.0 && x_0_p <= 1.0) && 
         (x_1_n >= 0.0 && x_1_n <= 1.0) && 
         (x_1_p >= 0.0 && x_1_p <= 1.0) ) { 
      /* Bounded solutions between 0 and 1.0 */
      Calculate_Bistable_Branches(x, x_0_n, x_0_p, x_1_n, x_1_p, 
                                  Vector_Stationarity_Lower,
                                  Vector_Stationarity_Inter,
                                  Vector_Stationarity_Upper, Table);

      Table->Vector_Model_Variables_Stationarity[0] = Vector_Stationarity_Upper[0];
      Table->Vector_Model_Variables_Stationarity[1] = Vector_Stationarity_Upper[1];

      Type_of_Coexistence = 3; /* Two coexistence solutions */
      printf("Bistability (Type_of_Coexistence = %g): Two coexisting solutions\n", 
      Type_of_Coexistence);   
    }
    else if ( (x_0_n >= 0.0 && x_0_n <= 1.0) && (x_1_n >= 0.0 && x_1_n <= 1.0) ){
      /* One single solution: mixed populations */ 
      Vector_Stationarity_Inter[0] =                  K_R * x_0_n;
      Vector_Stationarity_Lower[0] =                  K_R * x_0_n;
      Vector_Stationarity_Upper[0] =                  K_R * x_0_n;
      Table->Vector_Model_Variables_Stationarity[0] = K_R * x_0_n;

      Vector_Stationarity_Inter[1] =                  K_R * x_1_n;
      Vector_Stationarity_Lower[1] =                  K_R * x_1_n;
      Vector_Stationarity_Upper[1] =                  K_R * x_1_n;
      Table->Vector_Model_Variables_Stationarity[1] = K_R * x_1_n;
      
      Type_of_Coexistence = 2; /* One solution: mixed populations */    
      printf("One solution (Type_of_Coexistence = %g): coexistence of plamid free and plasmid carrying populations\n", 
                 Type_of_Coexistence);  
    }
    else if ( (x_0_p >= 0.0 && x_0_p <= 1.0) && (x_1_p >= 0.0 && x_1_p <= 1.0) ){
       /* One single solution: mixed populations */ 
       Vector_Stationarity_Inter[0] =                  K_R * x_0_p;
       Vector_Stationarity_Lower[0] =                  K_R * x_0_p;
       Vector_Stationarity_Upper[0] =                  K_R * x_0_p;
       Table->Vector_Model_Variables_Stationarity[0] = K_R * x_0_p;

       Vector_Stationarity_Inter[1] =                  K_R * x_1_p;
       Vector_Stationarity_Lower[1] =                  K_R * x_1_p;
       Vector_Stationarity_Upper[1] =                  K_R * x_1_p;
       Table->Vector_Model_Variables_Stationarity[1] = K_R * x_1_p;
       
       Type_of_Coexistence = 2; /* One solution: mixed populations */
       printf("One solution (Type_of_Coexistence = %g): coexistence of plamid free and plasmid carrying populations\n", 
              Type_of_Coexistence);  
    }
    else {
        Type_of_Coexistence = Calculate_Plasmid_Free_Solution(beta_0, D_0, 
                                                              Vector_Stationarity_Lower,
                                                              Vector_Stationarity_Inter,
                                                              Vector_Stationarity_Upper, 
                                                              Table);
    }       
  }

  return(Type_of_Coexistence);
}

void Fixed_Points_OneCell_Bis (Parameter_Table * Table, 
                              double * Vector_Stationarity_Lower,
                              double * Vector_Stationarity_Inter,
                              double * Vector_Stationarity_Upper )

{   
  double Type_of_Coexistence; 

  if(Table->No_of_CELLS > 1) {
    printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
    printf("Thefore, the dynamics should occurs in one single patch or local population\n");
    printf("However, No of CELLS has been set to be larger than 1 (M = %d)\n", 
            Table->No_of_CELLS);
    printf("Press Key...\n");
    Print_Press_Key(1,1,"The program will safely exit\n");
  exit(0);
  }

  Type_of_Coexistence = Calculate_Type_of_Coexistence_Bis (Table);  
}


double Calculate_Type_of_Coexistence_Bis (Parameter_Table * Table) 
                               
  
{ 
  int i;
  double x, x_0, x_p, x_n, x_1_p, x_1_n, x_0_n, x_0_p; 
  double K_R = Table->K_R;
  double Type_of_Coexistence;

  double * Vector_Stationarity_Lower =  Table->Vector_Model_Variables_MultiStability[0]; 
  double * Vector_Stationarity_Inter =  Table->Vector_Model_Variables_MultiStability[1];
  double * Vector_Stationarity_Upper =  Table->Vector_Model_Variables_MultiStability[2];

  if(Table->No_of_CELLS > 1) {
    printf("Here, fixed points are only analytically possible if No_of_CELLS is 1\n");
    printf("Thefore, the dynamics should occurs in one single patch or local population\n");
    printf("However, No of CELLS has been set to be larger than 1 (M = %d)\n", 
              Table->No_of_CELLS);
    printf("Press Key...\n");
    Print_Press_Key(1,1,"The program will safely exit\n");
    exit(0);
  }

  double D_0, D_1, beta_0, beta_1, epsilon, cost, chi, p_C, Dis;
  double A_0, A_1, B_0, B_1;
  double a, b, c;

  /* Relative rates to the encounter rate, Gamma (coded as Table->Lambda_R_1) */
  /* D_0 = (Delta_R_0 + Delta_R_1)/Gamma */
  D_0 = (Table->Delta_R_0 + Table->Delta_R_1) / Table->Lambda_R_1; 
  /* D_1 = (Delta_R_0 + Delta_R_1 * (1.0 - Table->Nu_C_0))/Gamma, where Nu_C_0 is the resistence */
  D_1 = (Table->Delta_R_0 + Table->Delta_R_1*(1.0 - Table->Nu_C_0)) / Table->Lambda_R_1; 

  beta_0 = Table->Beta_R / Table->Lambda_R_1;
  beta_1 = Table->Beta_R * (1.0 - Table->Alpha_C_0) / Table->Lambda_R_1;

  epsilon = Table->p_1;         /* Segregation Error */ 
  chi     = Table->Chi_C_0;     /* Probability of Death */
  p_C     = Table->p_2;         /* Probability of Competition Induced Death */
  cost    = Table->Alpha_C_0;   /* Cost of Resistance */ 
  
  A_0 = 1.0 - D_0/beta_0;
  A_1 = 1.0 - D_1/ (beta_1 * (1.0 - epsilon));

  B_0 = 1.0 + (p_C + chi) / beta_0;
  B_1 = 1.0 + (p_C - chi) / (beta_1 * (1.0 - epsilon));

  // Calculate 'a'
  a =   ( 1.0/B_1 - (B_0+epsilon*(1.0-cost)) + B_1*epsilon*(1.0-cost) );
  
  // Calculate 'b'
  b =   ( A_0/A_1 - 2.0/B_1 + B_0+epsilon*(1.0-cost) - B_1/A_1*epsilon*(1.0-cost) ); 

  // Calculate 'c'
  c = - ( A_0/A_1 - 1.0/B_1 );

  Dis = b * b - 4.0 * a * c;
  
  x = 1.0 - D_0 / beta_0;
 
  if(Dis < 0) {
    printf("The discriminant is negative\n");
    printf("Mixed solution is not possible\n");
    
    Type_of_Coexistence = Calculate_Plasmid_Free_Solution(beta_0, D_0, 
                                                          Vector_Stationarity_Lower,
                                                          Vector_Stationarity_Inter,
                                                          Vector_Stationarity_Upper, Table);     
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

    x_1_n = A_1 * x_n;
    x_0_n = A_1/B_1 * (1.0 - x_n); 
    
    x_1_p = A_1 * x_p;
    x_0_p = A_1/B_1 * (1.0 - x_p); 

   if ( (x_0_n >= 0.0 && x_0_n <= 1.0) && 
        (x_0_p >= 0.0 && x_0_p <= 1.0) && 
        (x_1_n >= 0.0 && x_1_n <= 1.0) && 
        (x_1_p >= 0.0 && x_1_p <= 1.0) ) { 
     /* Bounded solutions between 0 and 1.0 */
     /* Two coexisting solutions: mixed populations */      
      Calculate_Bistable_Branches(x, x_0_n, x_0_p, x_1_n, x_1_p, 
                                  Vector_Stationarity_Lower,
                                  Vector_Stationarity_Inter,
                                  Vector_Stationarity_Upper, Table);

      Table->Vector_Model_Variables_Stationarity[0] = Vector_Stationarity_Upper[0];
      Table->Vector_Model_Variables_Stationarity[1] = Vector_Stationarity_Upper[1];

      Type_of_Coexistence = 3; /* Two coexistence solutions */
      printf("Bistability (Type_of_Coexistence = %g): Two coexisting solutions\n", 
             Type_of_Coexistence); 
   }
   else if ( (x_0_n >= 0.0 && x_0_n <= 1.0) && (x_1_n >= 0.0 && x_1_n <= 1.0) ){
     /* One single solution: mixed populations */ 
     Vector_Stationarity_Inter[0] =                  K_R * x_0_n;
     Vector_Stationarity_Lower[0] =                  K_R * x_0_n;
     Vector_Stationarity_Upper[0] =                  K_R * x_0_n;
     Table->Vector_Model_Variables_Stationarity[0] = K_R * x_0_n;

     Vector_Stationarity_Inter[1] =                  K_R * x_1_n;
     Vector_Stationarity_Lower[1] =                  K_R * x_1_n;
     Vector_Stationarity_Upper[1] =                  K_R * x_1_n;
     Table->Vector_Model_Variables_Stationarity[1] = K_R * x_1_n;
     
     Type_of_Coexistence = 2; /* One solution: mixed populations */    
     printf("One solution (Type_of_Coexistence = %g): coexistence of plamid free and plasmid carrying populations\n", 
                Type_of_Coexistence);  
   }
   else if ( (x_0_p >= 0.0 && x_0_p <= 1.0) && (x_1_p >= 0.0 && x_1_p <= 1.0) ){
      /* One single solution: mixed populations */ 
      Vector_Stationarity_Inter[0] =                  K_R * x_0_p;
      Vector_Stationarity_Lower[0] =                  K_R * x_0_p;
      Vector_Stationarity_Upper[0] =                  K_R * x_0_p;
      Table->Vector_Model_Variables_Stationarity[0] = K_R * x_0_p;

      Vector_Stationarity_Inter[1] =                  K_R * x_1_p;
      Vector_Stationarity_Lower[1] =                  K_R * x_1_p;
      Vector_Stationarity_Upper[1] =                  K_R * x_1_p;
      Table->Vector_Model_Variables_Stationarity[1] = K_R * x_1_p;
      
      Type_of_Coexistence = 2; /* One solution: mixed populations */
      printf("One solution (Type_of_Coexistence = %g): coexistence of plamid free and plasmid carrying populations\n", 
             Type_of_Coexistence);  
   }
   else {
       Type_of_Coexistence = Calculate_Plasmid_Free_Solution(beta_0, D_0, 
                                                             Vector_Stationarity_Lower,
                                                             Vector_Stationarity_Inter,
                                                             Vector_Stationarity_Upper, 
                                                             Table);
   }
  }

  return(Type_of_Coexistence);
}

double Calculate_Plasmid_Free_Solution(double beta_0, double D_0,
                                       double * Vector_Stationarity_Lower,
                                       double * Vector_Stationarity_Inter,
                                       double * Vector_Stationarity_Upper, 
                                       Parameter_Table * Table) 
{
  int i;
  double Type_of_Coexistence;

  if (beta_0 <= D_0){
  /* The system goes to full collapse: extinction */
    for(i=0; i<Table->MODEL_STATE_VARIABLES; i++) {
      Vector_Stationarity_Inter[i] =                  0.0;
      Vector_Stationarity_Lower[i] =                  0.0;
      Vector_Stationarity_Upper[i] =                  0.0;
      Table->Vector_Model_Variables_Stationarity[i] = 0.0;
    }
  Type_of_Coexistence = 0.0;
  printf("Extinction (Type_of_Coexistence = %g): Collapse of both populations\n", 
          Type_of_Coexistence);  
  }
  else {
  /* One population: only Plasmid Free */
    Plasmid_Free_Solution(Table, Vector_Stationarity_Lower,
                                 Vector_Stationarity_Inter,
                                 Vector_Stationarity_Upper);

    Type_of_Coexistence = 1.0;
    printf("One solution (Type_of_Coexistence = %g): Plasmid-free equilibrium\n", 
            Type_of_Coexistence);                           
  }

  return(Type_of_Coexistence);
}


void Plasmid_Free_Solution(Parameter_Table * Table, 
                           double * Vector_Stationarity_Lower,
                           double * Vector_Stationarity_Inter,
                           double * Vector_Stationarity_Upper)
{
  double x;
  
  x = 1.0 - (Table->Delta_R_0 + Table->Delta_R_1) / Table->Beta_R;
  
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

double Calculate_Stability_Stationary_Point( Parameter_Table * Table)
{
  double Type_of_Stability;

  Type_of_Stability = Calculate_Type_of_Stability_2D( Table );

  return(Type_of_Stability);
}

void Calculate_Bistable_Branches(double x, double x_0_n, double x_0_p, double x_1_n, double x_1_p, 
                                double * Vector_Stationarity_Lower,
                                double * Vector_Stationarity_Inter,
                                double * Vector_Stationarity_Upper, 
                                Parameter_Table * Table)
{
  /* Since bistability involves the coexistence of two stable solutions, this scenario can 
     be characterized by three branches: two stable branches separated by an instable one: 
     the upper branch (stable), the lower branch (stable), and the intermediate branch 
     (unstable).
  */
  double K_R = (double)Table->K_R;

  /* Plasmid-free strain */
  Vector_Stationarity_Inter[0] =                  K_R * MAX(x_0_p, x_0_n);  /* Unstable */
  Vector_Stationarity_Lower[0] =                  K_R * MIN(x_0_p, x_0_n);  /* Stable */
  Vector_Stationarity_Upper[0] =                  K_R * MAX(x, 0.0);        /* Stable */

  /* Plasmid-carrying strain: (upper and lower branches refer only to the plasmid-free solution) */
  Vector_Stationarity_Inter[1] =                  K_R * MIN(x_1_n, x_1_p);  /* Unstable */  
  Vector_Stationarity_Upper[1] =                  0.0;                      /* Stable */
  Vector_Stationarity_Lower[1] =                  K_R * MAX(x_1_n, x_1_p);  /* Stable */   

  /* Checking Stationarity and Type of Stability of the three branches */
} 

