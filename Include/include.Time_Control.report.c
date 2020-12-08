fprintf(fp, "\n");
fprintf(fp, "Time Control Input Arguments\n");
fprintf(fp, "-tn  %d [Number of Time Points       = %d]\n", I_Time, I_Time );
fprintf(fp, "-t0  %g [Initial Time                = %g]\n", Time_0, Time_0 ); 
fprintf(fp, "-t1  %g [Final Time                  = %g]\n", Time_1, Time_1 ); 
fprintf(fp, "-tU  %g [Time Unit (Time_Scale_Unit) = %g]\n", Time_Scale_Unit, Time_Scale_Unit); 
fprintf(fp, "-tR  %d [Realizations (Stochastic)   = %d]\n", Realizations, Realizations); 
fprintf(fp, "-t3  %g [Time Interval (Delta_T)     = %g]\n", Delta_T, Delta_T); 
fprintf(fp, "-t4  %d [TYPE of TIME DEPENDENCE     = %d]\n",
	TYPE_of_TIME_DEPENDENCE, TYPE_of_TIME_DEPENDENCE); 
fprintf(fp, "-tD  %d [DISCARTING EXTINCTIONSE     = %d]\n",
	DISCARTING_EXTINCTIONS, DISCARTING_EXTINCTIONS); 




