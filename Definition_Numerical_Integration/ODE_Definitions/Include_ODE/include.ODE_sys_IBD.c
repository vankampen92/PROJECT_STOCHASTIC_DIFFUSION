/* c: Children (0-5) */
dydt[n0S] = - Force_of_Infection_0  * y[n0S];                                                                                                                                /* n0S */

dydt[n0E] =   Force_of_Infection_0  * y[n0S] - Table->Sigma * y[n0E];                                                                                                         /* n0E */

dydt[n0I1] = Table->Sigma * y[n0E] - Table->Gamma_1 * y[n0I1];                                                                                                               /* n0I1 */

dydt[n0I2] = Table->p_0 * Table->Gamma_1 * y[n0I1] - Table->Gamma_2 * y[n0I2]  - Table->Alpha_0 * y[n0I2];                                                                   /* n0I2 */

dydt[n0A] = (1.0 - Table->p_0) * Table->Gamma_1 * y[n0I1] - Table->Gamma_2 * y[n0A] - Table->Kappa * y[n0A];                                                                 /* n0A */

dydt[n0Ad] = Table->Kappa * y[n0A] - Table->Gamma_2 * y[n0Ad];                                                                                                               /* n0Ad */

dydt[n0Y] = Table->Alpha_0 * y[n0I2] - Table->Gamma_2 * y[n0Y] - Table->Delta_0 * y[n0Y];                                                                                    /* n0Y */

dydt[n0R] = Table->Gamma_2 * (y[n0A] + y[n0I2] + y[n0Y] + y[n0Ad]);                                                                                                          /* n0R */

dydt[a0I] = Table->p_0 * Table->Gamma_1 * y[n0I1] + Table->Kappa * y[n0A];                                                                                                   /* a0I */

dydt[a0R] = Table->Gamma_2 * (y[n0I2] + y[n0Y] + y[n0Ad]);                                                                                                                   /* a0R */

dydt[a0D] = Table->Delta_0 * y[n0Y];                                                                                                                                         /* a0D */

/* s: Students (6-25) */
dydt[n1S] = - Force_of_Infection_1 * y[n1S];                                                                                                                                 /* n1S */

dydt[n1E] =   Force_of_Infection_1 * y[n1S] - Table->Sigma * y[n1E];                                                                                                         /* n1E */

dydt[n1I1] = Table->Sigma * y[n1E] - Table->Gamma_1 * y[n1I1];                                                                                                               /* n1I1 */

dydt[n1I2] = Table->p_1 * Table->Gamma_1 * y[n1I1] - Table->Gamma_2 * y[n1I2]  - Table->Alpha_1 * y[n1I2];                                                                   /* n1I2 */

dydt[n1A] = (1.0 - Table->p_1) * Table->Gamma_1 * y[n1I1] - Table->Gamma_2 * y[n1A] - Table->Kappa * y[n1A];                                                                 /* n1A */

dydt[n1Ad] = Table->Kappa * y[n1A] - Table->Gamma_2 * y[n1Ad];                                                                                                               /* n1Ad */
 
dydt[n1Y] = Table->Alpha_1 * y[n1I2] - Table->Gamma_2 * y[n1Y] - Table->Delta_1 * y[n1Y];                                                                                    /* n1Y */

dydt[n1R] = Table->Gamma_2 * (y[n1A] + y[n1I2] + y[n1Y] + y[n1Ad]) ;                                                                                                         /* n1R */

dydt[a1I] = Table->p_1 * Table->Gamma_1 * y[n1I1] + Table->Kappa * y[n1A];                                                                                                   /* a1I */

dydt[a1R] = Table->Gamma_2 * (y[n1I2] + y[n1Y] + y[n1Ad]);                                                                                                                   /* a1R */

dydt[a1D] = Table->Delta_1 * y[n1Y];                                                                                                                                         /* a1D */

/* w: Adults (26-65) */
dydt[n2S] = - Force_of_Infection_2 * y[n2S];                                                                                                                                 /* n2S */

dydt[n2E] =   Force_of_Infection_2 * y[n2S] - Table->Sigma * y[n2E];                                                                                                         /* n2E */

dydt[n2I1] = Table->Sigma * y[n2E] - Table->Gamma_1 * y[n2I1];                                                                                                               /* n2I1 */

dydt[n2I2] = Table->p_2 * Table->Gamma_1 * y[n2I1] - Table->Gamma_2 * y[n2I2]  - Table->Alpha_2 * y[n2I2];                                                                   /* n2I2 */

dydt[n2A] = (1.0 - Table->p_2) * Table->Gamma_1 * y[n2I1] - Table->Gamma_2 * y[n2A] - Table->Kappa * y[n2A];                                                                  /* n2A */

dydt[n2Ad] = Table->Kappa * y[n2A] - Table->Gamma_2 * y[n2Ad];                                                                                                                /* n2Ad */
 
dydt[n2Y] = Table->Alpha_2 * y[n2I2] - Table->Gamma_2 * y[n2Y] - Table->Delta_2 * y[n2Y];                                                                                     /* n2Y */

dydt[n2R] = Table->Gamma_2 * (y[n2A] + y[n2I2] + y[n2Y] + y[n2Ad]);                                                                                                           /* n2R */

dydt[a2I] = Table->p_2 * Table->Gamma_1 * y[n2I1] + Table->Kappa * y[n2A];                                                                                                    /* a2I */

dydt[a2R] = Table->Gamma_2 * (y[n2I2] + y[n2Y] + y[n2Ad]);                                                                                                                    /* a2R */

dydt[a2D] = Table->Delta_2 * y[n2Y];                                                                                                                                          /* a2D */

/* a: Seniors (65-100) */
dydt[n3S] = - Force_of_Infection_3 * y[n3S];                                                                                                                                 /* n3S */

dydt[n3E] =   Force_of_Infection_3 * y[n3S] - Table->Sigma * y[n3E];                                                                                                          /* n3E */

dydt[n3I1] = Table->Sigma * y[n3E] - Table->Gamma_1 * y[n3I1];                                                                                                               /* n3I1 */

dydt[n3I2] = Table->p_3 * Table->Gamma_1 * y[n3I1] - Table->Gamma_2 * y[n3I2]  - Table->Alpha_3 * y[n3I2];                                                                   /* n3I2 */

dydt[n3A] = (1.0 - Table->p_3) * Table->Gamma_1 * y[n3I1] - Table->Gamma_2 * y[n3A] - Table->Kappa * y[n3A];                                                                 /* n3A */

dydt[n3Ad] = Table->Kappa * y[n3A] - Table->Gamma_2 * y[n3Ad];                                                                                                               /* n3Ad */
 
dydt[n3Y] = Table->Alpha_3 * y[n3I2] - Table->Gamma_2 * y[n3Y] - Table->Delta_3 * y[n3Y];                                                                                    /* n3Y */

dydt[n3R] = Table->Gamma_2 * (y[n3A] + y[n3I2] + y[n3Y] + y[n3Ad]) ;                                                                                                         /* n3R */

dydt[a3I] = Table->p_3 * Table->Gamma_1 * y[n3I1] + Table->Kappa * y[n3A];                                                                                                   /* a3I */

dydt[a3R] = Table->Gamma_2 * (y[n3I2] + y[n3Y] + y[n3Ad]) ;                                                                                                                  /* a3R */

dydt[a3D] = Table->Delta_3 * y[n3Y];                                                                                                                                         /* a3D */
