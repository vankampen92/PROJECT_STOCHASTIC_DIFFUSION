#include <MODEL.h>

void Event_Adjacence_List_Initialization(Community ** PATCH,
					 Parameter_Model * P)
{
  int i, j, k, Sp, n, m, no;
  int TOTAL_No_of_EVENTS;
  
  int ** AD; 

  no  = P->No_of_CELLS;
  m   = P->TOTAL_No_of_EVENTS; 
  Sp  = P->No_of_RESOURCES; 

  assert(no == 1);

  for(i=0; i<no; i++){
    TOTAL_No_of_EVENTS = 0;

    AD = PATCH[i]->Event_Adjacence_List;
     
    for(j=0; j<Sp; j++) {
      
      /* Attack   */
      n = 0;
      for(k=0; k<Sp; k++) { 
        AD[j*P->No_of_EVENTS][n++]    = k*P->No_of_EVENTS;
      }
      AD[j*P->No_of_EVENTS][n++]      = j*P->No_of_EVENTS + 1;
      AD[j*P->No_of_EVENTS][n++]      = Sp*P->No_of_EVENTS;
      AD[j*P->No_of_EVENTS][m] = n;
      TOTAL_No_of_EVENTS++;

      /* Handling */
      n=0;
      for(k=0; k<Sp; k++) { 
        AD[j*P->No_of_EVENTS+1][n++]    = k*P->No_of_EVENTS;
      }
      AD[j*P->No_of_EVENTS+1][n++]      = j*P->No_of_EVENTS + 1;
      AD[j*P->No_of_EVENTS+1][n++]      = Sp*P->No_of_EVENTS;  
      AD[j*P->No_of_EVENTS+1][m] = n;
      TOTAL_No_of_EVENTS++;  
    }   

    assert(TOTAL_No_of_EVENTS == (m-2));

    /* Out-Migration */
    n = 0;
    for(k=0; k<Sp; k++) { 
        AD[Sp*P->No_of_EVENTS][n++]  = k*P->No_of_EVENTS;
    }
    AD[Sp*P->No_of_EVENTS][n++]  = Sp*P->No_of_EVENTS;
    AD[Sp*P->No_of_EVENTS][m] = n;
    TOTAL_No_of_EVENTS++;

    /* In-Migration */
    n = 0;
    for(k=0; k<Sp; k++) { 
        AD[Sp*P->No_of_EVENTS][n++]  = k*P->No_of_EVENTS;
    }
    AD[Sp*P->No_of_EVENTS][n++]  = Sp*P->No_of_EVENTS;
    AD[Sp*P->No_of_EVENTS][m] = n;
    TOTAL_No_of_EVENTS++;

    assert(TOTAL_No_of_EVENTS == m);
  }
}



