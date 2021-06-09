/* Human Parameters */ 
	    
	    case 'T': /* Maximum and Minimum Transmission Rate */
	      if(argv[argcount][2]=='Y')
		fprintf(fp,"-TY  %d  ",Table->Tr->TYPE_of_TREND);
	      else if(argv[argcount][2]=='0')
		fprintf(fp,"-T0  %g  ",Table->Tr->Tr_Time_0);
              else if(argv[argcount][2]=='1')
		fprintf(fp,"-T1  %g  ",Table->Tr->Tr_Time_1);
              else if(argv[argcount][2]=='0')
		fprintf(fp,"-T2  %g  ",Table->Tr->Tr_Time_2);
              else if(argv[argcount][2]=='1')
		fprintf(fp,"-T3  %g  ",Table->Tr->Tr_Time_3);
              else if(argv[argcount][2]=='0')
		fprintf(fp,"-T4  %g  ",Table->Tr->Tr_Time_4);
              else if(argv[argcount][2]=='1')
		fprintf(fp,"-T5  %g  ",Table->Tr->Tr_Time_5);
              else if(argv[argcount][2]=='0')
		fprintf(fp,"-T6  %g  ",Table->Tr->Tr_Time_6);
              else if(argv[argcount][2]=='1')
		fprintf(fp,"-T7  %g  ",Table->Tr->Tr_Time_7);
              else if(argv[argcount][2]=='0')
		fprintf(fp,"-T8  %g  ",Table->Tr->Tr_Time_8);
              else if(argv[argcount][2]=='1')
		fprintf(fp,"-T9  %g  ",Table->Tr->Tr_Time_9);
	      else if(argv[argcount][2]=='i')
		fprintf(fp,"-Ti  %g  ",Table->Tr->Tr_Time_i);
	      else if(argv[argcount][2]=='M')
		fprintf(fp,"-TM  %g  ",Table->Tr->Tr_value_i);
              else if(argv[argcount][2]=='J')
		fprintf(fp,"-TJ  %d  ", Table->Tr->Tr_No_of_Jumps);
	      else if(argv[argcount][2]=='P')
		fprintf(fp,"-TP  %d  ", Table->Tr->Tr_Input_Parameter);
              else if(argv[argcount][2]=='v') {
		if(argv[argcount][3]=='0')
		  fprintf(fp,"-Tv0  %g ", Table->Tr->Tr_value_0);
		else if(argv[argcount][3]=='1')
		  fprintf(fp,"-Tv1  %g ", Table->Tr->Tr_value_1);
		else if(argv[argcount][3]=='2')
		  fprintf(fp,"-Tv2  %g ", Table->Tr->Tr_value_2);
		else if(argv[argcount][3]=='3')
		  fprintf(fp,"-Tv3  %g ", Table->Tr->Tr_value_3);
		else if(argv[argcount][3]=='4')
		  fprintf(fp,"-Tv4  %g ", Table->Tr->Tr_value_4);
		else if(argv[argcount][3]=='5')
		  fprintf(fp,"-Tv5  %g ", Table->Tr->Tr_value_5);
		else if(argv[argcount][3]=='6')
		  fprintf(fp,"-Tv6  %g ", Table->Tr->Tr_value_6);
		else if(argv[argcount][3]=='7')
		  fprintf(fp,"-Tv7  %g ", Table->Tr->Tr_value_7);
		else if(argv[argcount][3]=='8')
		  fprintf(fp,"-Tv8  %g ", Table->Tr->Tr_value_8);
		else if(argv[argcount][3]=='9')
		  fprintf(fp,"-Tv9  %g ", Table->Tr->Tr_value_9);
		else {
		  printf(" Error in Trend Control Input arguments\n");
		  printf(" The program will exit\n"); 
		  exit(0); 
		}
	      }
	      else {
		printf(" Error in Trend Control Input arguments\n");
		printf(" The program will exit\n"); 
		exit(0);
	      }
	      skip++;
	      break;
