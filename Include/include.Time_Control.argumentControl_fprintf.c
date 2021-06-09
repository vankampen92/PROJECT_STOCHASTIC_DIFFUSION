/* Time Control Definition */ 
	    
	    case 't': 
	      if(argv[argcount][2]=='n')
		fprintf(fp,"-tn  %d  ",Table->T->I_Time);
              else if(argv[argcount][2]=='0')
		fprintf(fp,"-t0  %g  ",Table->T->Time_0);
              else if(argv[argcount][2]=='1')
		fprintf(fp,"-t1  %g  ",Table->T->Time_1);
              else if(argv[argcount][2]=='U')
		fprintf(fp,"-tu  %g  ",Table->T->Time_Scale_Unit);
              else if(argv[argcount][2]=='R')
		fprintf(fp,"-tR  %d  ",Table->T->Realizations);
              else if(argv[argcount][2]=='3')
		fprintf(fp,"-t3  %g  ",Table->T->Delta_T);
              else if(argv[argcount][2]=='4')
		fprintf(fp,"-t4  %d  ",Table->T->TYPE_of_TIME_DEPENDENCE);
              else if(argv[argcount][2]=='E')
		fprintf(fp,"-tE  %g  ",Table->T->EPSILON);
              else if(argv[argcount][2]=='D')
		fprintf(fp,"-tD  %d  ", Table->T->DISCARTING_EXTINCTIONS);
              else exit(0);
	      skip++;
	      break;
	    
	    
	    
