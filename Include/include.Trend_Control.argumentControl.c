 
            case 'T': /* Maximum and Minimum Transmission Rate */
	      if(argv[argcount][2]=='Y')
		sscanf(argv[argcount+1],"%d",&TYPE_of_TREND);
	      else if(argv[argcount][2]=='0')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_0);
              else if(argv[argcount][2]=='1')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_1);
	      else if(argv[argcount][2]=='i')
		sscanf(argv[argcount+1],"%lf",&Tr_Time_i);
              else if(argv[argcount][2]=='I')
		sscanf(argv[argcount+1],"%lf",&Tr_value_0);
              else if(argv[argcount][2]=='F')
		sscanf(argv[argcount+1],"%lf",&Tr_value_1);
	      else if(argv[argcount][2]=='M')
		sscanf(argv[argcount+1],"%lf",&Tr_value_i);
              else if(argv[argcount][2]=='J')
		sscanf(argv[argcount+1],"%d", &Tr_No_of_Jumps);
	      else if(argv[argcount][2]=='P')
		sscanf(argv[argcount+1],"%d", &Tr_Input_Parameter);
	      else exit(0);
	      skip++;
	      break;

	    
