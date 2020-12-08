/* c: Children (0-5) */
dydt[n0S] +=  In_Mu(Table, n0S, 0, j, y) - Out_Mu_Per_Capita(Table, 0, j) * y[n0S];              /* n0S */

dydt[n0E] +=  In_Mu(Table, n0E, 0, j, y) - Out_Mu_Per_Capita(Table, 0, j) * y[n0E];               /* n0E */

dydt[n0I1] += In_Mu(Table, n0I1, 0, j, y) - Out_Mu_Per_Capita(Table, 0, j) * y[n0I1];             /* n0I1 */

dydt[n0I2] += (1.0-Table->Eps_I)*(In_Mu(Table, n0I2, 0, j, y) - Out_Mu_Per_Capita(Table, 0, j) * y[n0I2]);             /* n0I2 */

dydt[n0A] +=  In_Mu(Table, n0A, 0, j, y) - Out_Mu_Per_Capita(Table, 0, j) * y[n0A];               /* n0A */

dydt[n0Ad] += (1.0-Table->Eps_I)*(In_Mu(Table, n0Ad, 0, j, y) - Out_Mu_Per_Capita(Table, 0, j) * y[n0Ad]);             /* n0Ad */

  //dydt[n0Y] += ;                                                                     /* n0Y */
  
dydt[n0R] +=  In_Mu(Table, n0R, 0, j, y) - Out_Mu_Per_Capita(Table, 0, j) * y[n0R];                /* n0R */

							       
/* s: Students (6-25) */
dydt[n1S] +=  In_Mu(Table, n1S, 1, j, y) - Out_Mu_Per_Capita(Table, 1, j) * y[n1S];              /* n1S */

dydt[n1E] +=  In_Mu(Table, n1E, 1, j, y) - Out_Mu_Per_Capita(Table, 1, j) * y[n1E];               /* n1E */

dydt[n1I1] += In_Mu(Table, n1I1, 1, j, y) - Out_Mu_Per_Capita(Table, 1, j) * y[n1I1];             /* n1I1 */

dydt[n1I2] += (1.0-Table->Eps_I)*(In_Mu(Table, n1I2, 1, j, y) - Out_Mu_Per_Capita(Table, 1, j) * y[n1I2]);             /* n1I2 */

dydt[n1A] +=  In_Mu(Table, n1A, 1, j, y) - Out_Mu_Per_Capita(Table, 1, j) * y[n1A];               /* n1A */

dydt[n1Ad] += (1.0-Table->Eps_I)*(In_Mu(Table, n1Ad, 1, j, y) - Out_Mu_Per_Capita(Table, 1, j) * y[n1Ad]);             /* n1Ad */

  //dydt[n1Y] += ;                                                                     /* n1Y */
  
dydt[n1R] +=  In_Mu(Table, n1R, 1, j, y) - Out_Mu_Per_Capita(Table, 1, j) * y[n1R];                /* n1R */

							       
/* w: Adults (26-65) */
dydt[n2S] +=  In_Mu(Table, n2S, 2, j, y) - Out_Mu_Per_Capita(Table, 2, j) * y[n2S];              /* n2S */

dydt[n2E] +=  In_Mu(Table, n2E, 2, j, y) - Out_Mu_Per_Capita(Table, 2, j) * y[n2E];               /* n2E */

dydt[n2I1] += In_Mu(Table, n2I1, 2, j, y) - Out_Mu_Per_Capita(Table, 2, j) * y[n2I1];             /* n2I2 */

dydt[n2I2] += (1.0-Table->Eps_I)*(In_Mu(Table, n2I2, 2, j, y) - Out_Mu_Per_Capita(Table, 2, j) * y[n2I2]);             /* n2I2 */

dydt[n2A] +=  In_Mu(Table, n2A, 2, j, y) - Out_Mu_Per_Capita(Table, 2, j) * y[n2A];               /* n2A */

dydt[n2Ad] += (1.0-Table->Eps_I)*(In_Mu(Table, n2Ad, 2, j, y) - Out_Mu_Per_Capita(Table, 2, j) * y[n2Ad]);             /* n2Ad */

  //dydt[n2Y] += ;                                                                     /* n2Y */
  
dydt[n2R] +=  In_Mu(Table, n2R, 2, j, y) - Out_Mu_Per_Capita(Table, 2, j) * y[n2R];                /* n2R */
							       

/* a: Seniors (65-100) */
dydt[n3S] +=  In_Mu(Table, n3S, 3, j, y) - Out_Mu_Per_Capita(Table, 3, j) * y[n3S];              /* n3S */

dydt[n3E] +=  In_Mu(Table, n3E, 3, j, y) - Out_Mu_Per_Capita(Table, 3, j) * y[n3E];               /* n3E */

dydt[n3I1] += In_Mu(Table, n3I1, 3, j, y) - Out_Mu_Per_Capita(Table, 3, j) * y[n3I1];             /* n3I3 */

dydt[n3I2] += (1.0-Table->Eps_I)*(In_Mu(Table, n3I2, 3, j, y) - Out_Mu_Per_Capita(Table, 3, j) * y[n3I2]);             /* n3I3 */

dydt[n3A] +=  In_Mu(Table, n3A, 3, j, y) - Out_Mu_Per_Capita(Table, 3, j) * y[n3A];               /* n3A */

dydt[n3Ad] += (1.0-Table->Eps_I)*(In_Mu(Table, n3Ad, 3, j, y) - Out_Mu_Per_Capita(Table, 3, j) * y[n3Ad]);             /* n3Ad */

  //dydt[n3Y] += ;                                                                     /* n3Y */
  
dydt[n3R] +=  In_Mu(Table, n3R, 3, j, y) - Out_Mu_Per_Capita(Table, 3, j) * y[n3R];                /* n3R */



