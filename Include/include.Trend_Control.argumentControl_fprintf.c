/* Human Parameters */ 
	    
	    case 'T': /* Maximum and Minimum Transmission Rate */
	      if(argv[argcount][2]=='Y')
		fprintf(fp,"-TY  %d  ",Table->Tr->TYPE_of_TREND);
	      else if(argv[argcount][2]=='0')
		fprintf(fp,"-T0  %g  ",Table->Tr->Tr_Time_0);
              else if(argv[argcount][2]=='1')
		fprintf(fp,"-T1  %g  ",Table->Tr->Tr_Time_1);
	      else if(argv[argcount][2]=='i')
		fprintf(fp,"-Ti  %g  ",Table->Tr->Tr_Time_i);
              else if(argv[argcount][2]=='I')
		fprintf(fp,"-TI  %g  ",Table->Tr->Tr_value_0);
              else if(argv[argcount][2]=='F')
		fprintf(fp,"-TF  %g  ",Table->Tr->Tr_value_1);
	      else if(argv[argcount][2]=='M')
		fprintf(fp,"-TM  %g  ",Table->Tr->Tr_value_i);
              else if(argv[argcount][2]=='J')
		fprintf(fp,"-TJ  %d  ", Table->Tr->Tr_No_of_Jumps);
	      else if(argv[argcount][2]=='P')
		fprintf(fp,"-TP  %d  ", Table->Tr->Tr_Input_Parameter);
	      else exit(0);
	      skip++;
	      break;
