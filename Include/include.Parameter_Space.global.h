  int No_of_PARAMETERS; /* No of parameters defining the actual dimension of
			                     the subparameter space */

  int I0, I1, I2, I3, I4, I5, I6, I7, I8, I9;
  int d0, d1, d2, d3, d4, d5, d6, d7, d8, d9;
  double MAX_P0, MAX_P1, MAX_P2, MAX_P3, MAX_P4, MAX_P5, MAX_P6, MAX_P7, MAX_P8, MAX_P9;
  double min_P0, min_P1, min_P2, min_P3, min_P4, min_P5, min_P6, min_P7, min_P8, min_P9;
  double Acc_P0, Acc_P1, Acc_P2, Acc_P3, Acc_P4, Acc_P5, Acc_P6, Acc_P7, Acc_P8, Acc_P9;

  int I10, I11, I12, I13, I14, I15, I16, I17, I18, I19;
  int d10, d11, d12, d13, d14, d15, d16, d17, d18, d19;
  double MAX_P10, MAX_P11, MAX_P12, MAX_P13, MAX_P14, MAX_P15, MAX_P16, MAX_P17, MAX_P18, MAX_P19;
  double min_P10, min_P11, min_P12, min_P13, min_P14, min_P15, min_P16, min_P17, min_P18, min_P19;
  double Acc_P10, Acc_P11, Acc_P12, Acc_P13, Acc_P14, Acc_P15, Acc_P16, Acc_P17, Acc_P18, Acc_P19;

  int I20, I21, I22, I23, I24, I25, I26, I27, I28, I29;
  int d20, d21, d22, d23, d24, d25, d26, d27, d28, d29;
  double MAX_P20, MAX_P21, MAX_P22, MAX_P23, MAX_P24, MAX_P25, MAX_P26, MAX_P27, MAX_P28, MAX_P29;
  double min_P20, min_P21, min_P22, min_P23, min_P24, min_P25, min_P26, min_P27, min_P28, min_P29;
  double Acc_P20, Acc_P21, Acc_P22, Acc_P23, Acc_P24, Acc_P25, Acc_P26, Acc_P27, Acc_P28, Acc_P29;

  int I30, I31, I32, I33, I34, I35, I36, I37, I38, I39;
  int d30, d31, d32, d33, d34, d35, d36, d37, d38, d39;
  double MAX_P30, MAX_P31, MAX_P32, MAX_P33, MAX_P34, MAX_P35, MAX_P36, MAX_P37, MAX_P38, MAX_P39;
  double min_P30, min_P31, min_P32, min_P33, min_P34, min_P35, min_P36, min_P37, min_P38, min_P39;
  double Acc_P30, Acc_P31, Acc_P32, Acc_P33, Acc_P34, Acc_P35, Acc_P36, Acc_P37, Acc_P38, Acc_P39;
  
  double TOLERANCE;
  int MAX_No_of_ITERATIONS;

  double ** Ranges;
  double * Acc;     //Accuracy for each parameter dimension
  int    * d;       //Discretization for each parameter dimension
  int    * Index;
