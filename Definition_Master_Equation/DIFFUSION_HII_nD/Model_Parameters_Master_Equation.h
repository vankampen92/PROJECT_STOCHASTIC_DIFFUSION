void assert_HLL_nD_Equal_Nus( Parameter_Table * Table );
unsigned long long int Func_Hyper(unsigned int D, unsigned int k);
void Groups(unsigned int D, int D_Max, int * V, int * s);
unsigned long long int General_Configuration_Table(unsigned int D, unsigned int n_0, 
                                                   int ** Configuration_Table);
int Configuration_to_i_Map(int * Configuration, int D, int n_0);
void i_to_Configuration_Map(int ** Configuration_Table, int i, int * Configuration, int D);
void Writing_General_Configuration_Table(int D,
                                         unsigned long long int S, 
                                         int ** Configuration_Table);
void Comparing_Two_Configuration_Tables(int D, unsigned long long int S, 
                                        int ** Configuration_Table_D5, 
                                        int ** Configuration_Table);
void Pickup_a_Configuration_from_Table(int D, int Index, int ** Configuration_Table);
int Vector_of_Adjacent_Configurations(int * Configuration, int D, int n_0, int d,  
                                      int ** Nei, int * Index_Nei);
void Generating_Network_of_Configurations(configuration ** Co, int ** Configuration_Table, 
                                          int D, int S, int n_0);
