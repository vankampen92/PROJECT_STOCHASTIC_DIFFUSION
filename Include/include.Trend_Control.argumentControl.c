 
            case 'T': /* Maximum and Minimum Transmission Rate */
	      if(argv[argcount][2]=='Y')
		sscanf(argv[argcount+1],"%d",&TYPE_of_TREND);
	      else if(argv[argcount][2]=='0')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_0);
              else if(argv[argcount][2]=='1')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_1);
              else if(argv[argcount][2]=='2')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_2);
              else if(argv[argcount][2]=='3')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_3);
              else if(argv[argcount][2]=='4')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_4);
              else if(argv[argcount][2]=='5')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_5);
              else if(argv[argcount][2]=='6')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_6);
              else if(argv[argcount][2]=='7')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_7);
              else if(argv[argcount][2]=='8')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_8);
              else if(argv[argcount][2]=='9')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_9);
	      else if(argv[argcount][2]=='i')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_i);
	      else if(argv[argcount][2]=='M')
		sscanf(argv[argcount+1],"%lf",&Tr_value_i);
              else if(argv[argcount][2]=='J')
		sscanf(argv[argcount+1],"%d", &Tr_No_of_Jumps);
	      else if(argv[argcount][2]=='P')
		sscanf(argv[argcount+1],"%d", &Tr_Input_Parameter);
	      else if(argv[argcount][2]=='v') {
		if(argv[argcount][3]=='0')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_0);
		else if(argv[argcount][3]=='1')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_1);
		else if(argv[argcount][3]=='2')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_2);
		else if(argv[argcount][3]=='3')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_3);
		else if(argv[argcount][3]=='4')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_4);
		else if(argv[argcount][3]=='5')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_5);
		else if(argv[argcount][3]=='6')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_6);
		else if(argv[argcount][3]=='7')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_7);
		else if(argv[argcount][3]=='8')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_8);
		else if(argv[argcount][3]=='9')
		  sscanf(argv[argcount+1],"%lf",&Tr_value_9);
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

	    
