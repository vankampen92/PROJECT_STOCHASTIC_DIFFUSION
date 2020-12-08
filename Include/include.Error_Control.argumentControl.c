/* Scan Parameters */

/* Error Model Parameters */
	case 'e':
              if(argv[argcount][2]=='0') 	 sscanf(argv[argcount+1], "%lf",&Err_0);

	      else if(argv[argcount][2]=='1') {
		if(argv[argcount][3]=='0') sscanf(argv[argcount+1], "%lf", &Err_10);
		else if(argv[argcount][3]=='1')  sscanf(argv[argcount+1], "%lf", &Err_11);
                else if(argv[argcount][3]=='2')  sscanf(argv[argcount+1], "%lf", &Err_12);
		else if(argv[argcount][3]=='3')  sscanf(argv[argcount+1], "%lf", &Err_13);
		else if(argv[argcount][3]=='\0') sscanf(argv[argcount+1], "%lf", &Err_1);
	      }

	      else if(argv[argcount][2]=='2')    sscanf(argv[argcount+1], "%lf",&Err_2);
	      else if(argv[argcount][2]=='3')    sscanf(argv[argcount+1], "%lf",&Err_3);
              else if(argv[argcount][2]=='4')    sscanf(argv[argcount+1], "%lf",&Err_4);
              else if(argv[argcount][2]=='5')    sscanf(argv[argcount+1], "%lf",&Err_5);
	      else if(argv[argcount][2]=='6')    sscanf(argv[argcount+1], "%lf",&Err_6);
	      else if(argv[argcount][2]=='7')    sscanf(argv[argcount+1], "%lf",&Err_7);
              else if(argv[argcount][2]=='8')    sscanf(argv[argcount+1], "%lf",&Err_8);
              else if(argv[argcount][2]=='9')    sscanf(argv[argcount+1], "%lf",&Err_9);

	      else if(argv[argcount][2]=='y')
		sscanf(argv[argcount+1], "%d",&TYPE_of_ERROR_FUNCTION);

              else if(argv[argcount][2]=='n')
		sscanf(argv[argcount+1], "%d",&No_of_ERROR_PARAMETERS);

              else if(argv[argcount][2]=='E')
		sscanf(argv[argcount+1], "%lf",&Error_Parameter_1);

              else if(argv[argcount][2]=='V')
		sscanf(argv[argcount+1], "%lf",&Error_Parameter_0);

              else if(argv[argcount][2]=='P') {
		if(argv[argcount][3]=='0')      sscanf(argv[argcount+1], "%d", &E_Var_0);
		else if(argv[argcount][3]=='1') {
			if(argv[argcount][4]=='0')    sscanf(argv[argcount+1], "%d", &E_Var_10);
			else if(argv[argcount][4]=='1') sscanf(argv[argcount+1], "%d", &E_Var_11);
			else if(argv[argcount][4]=='2') sscanf(argv[argcount+1], "%d", &E_Var_12);
			else if(argv[argcount][4]=='3') sscanf(argv[argcount+1], "%d", &E_Var_13);
			else if(argv[argcount][4]=='4') sscanf(argv[argcount+1], "%d", &E_Var_14);
			else if(argv[argcount][4]=='5') sscanf(argv[argcount+1], "%d", &E_Var_15);
			else if(argv[argcount][4]=='6') sscanf(argv[argcount+1], "%d", &E_Var_16);
			else if(argv[argcount][4]=='\0') sscanf(argv[argcount+1], "%d",&E_Var_1);
		}
    else if(argv[argcount][3]=='2') sscanf(argv[argcount+1], "%d", &E_Var_2);
		else if(argv[argcount][3]=='3') sscanf(argv[argcount+1], "%d", &E_Var_3);
		else if(argv[argcount][3]=='4') sscanf(argv[argcount+1], "%d", &E_Var_4);
		else if(argv[argcount][3]=='5') sscanf(argv[argcount+1], "%d", &E_Var_5);
		else if(argv[argcount][3]=='6') sscanf(argv[argcount+1], "%d", &E_Var_6);
		else if(argv[argcount][3]=='7') sscanf(argv[argcount+1], "%d", &E_Var_7);
		else if(argv[argcount][3]=='8') sscanf(argv[argcount+1], "%d", &E_Var_8);
		else if(argv[argcount][3]=='9') sscanf(argv[argcount+1], "%d", &E_Var_9);
	      }

              else if(argv[argcount][2]=='m') {
		if(argv[argcount][3]=='0')      sscanf(argv[argcount+1], "%lf", &E_min_0);
		else if(argv[argcount][3]=='1') {
			if(argv[argcount][4]=='0')    sscanf(argv[argcount+1], "%lf", &E_min_10);
			else if(argv[argcount][4]=='1') sscanf(argv[argcount+1], "%lf", &E_min_11);
			else if(argv[argcount][4]=='2') sscanf(argv[argcount+1], "%lf", &E_min_12);
			else if(argv[argcount][4]=='3') sscanf(argv[argcount+1], "%lf", &E_min_13);
			else if(argv[argcount][4]=='4') sscanf(argv[argcount+1], "%lf", &E_min_14);
			else if(argv[argcount][4]=='5') sscanf(argv[argcount+1], "%lf", &E_min_15);
			else if(argv[argcount][4]=='6') sscanf(argv[argcount+1], "%lf", &E_min_16);
			else if(argv[argcount][4]=='\0') sscanf(argv[argcount+1], "%lf",&E_min_1);
		}
    else if(argv[argcount][3]=='2') sscanf(argv[argcount+1], "%lf", &E_min_2);
		else if(argv[argcount][3]=='3') sscanf(argv[argcount+1], "%lf", &E_min_3);
		else if(argv[argcount][3]=='4') sscanf(argv[argcount+1], "%lf", &E_min_4);
		else if(argv[argcount][3]=='5') sscanf(argv[argcount+1], "%lf", &E_min_5);
		else if(argv[argcount][3]=='6') sscanf(argv[argcount+1], "%lf", &E_min_6);
		else if(argv[argcount][3]=='7') sscanf(argv[argcount+1], "%lf", &E_min_7);
		else if(argv[argcount][3]=='8') sscanf(argv[argcount+1], "%lf", &E_min_8);
		else if(argv[argcount][3]=='9') sscanf(argv[argcount+1], "%lf", &E_min_9);
	      }

              else if(argv[argcount][2]=='M') {
		if(argv[argcount][3]=='0')      sscanf(argv[argcount+1], "%lf", &E_MAX_0);
		else if(argv[argcount][3]=='1') {
			if(argv[argcount][4]=='0')    sscanf(argv[argcount+1], "%lf", &E_MAX_10);
			else if(argv[argcount][4]=='1') sscanf(argv[argcount+1], "%lf", &E_MAX_11);
			else if(argv[argcount][4]=='2') sscanf(argv[argcount+1], "%lf", &E_MAX_12);
			else if(argv[argcount][4]=='3') sscanf(argv[argcount+1], "%lf", &E_MAX_13);
			else if(argv[argcount][4]=='4') sscanf(argv[argcount+1], "%lf", &E_MAX_14);
			else if(argv[argcount][4]=='5') sscanf(argv[argcount+1], "%lf", &E_MAX_15);
			else if(argv[argcount][4]=='6') sscanf(argv[argcount+1], "%lf", &E_MAX_16);
			else if(argv[argcount][4]=='\0') sscanf(argv[argcount+1], "%lf",&E_MAX_1);
		}
		else if(argv[argcount][3]=='2') sscanf(argv[argcount+1], "%lf", &E_MAX_2);
		else if(argv[argcount][3]=='3') sscanf(argv[argcount+1], "%lf", &E_MAX_3);
		else if(argv[argcount][3]=='4') sscanf(argv[argcount+1], "%lf", &E_MAX_4);
		else if(argv[argcount][3]=='5') sscanf(argv[argcount+1], "%lf", &E_MAX_5);
		else if(argv[argcount][3]=='6') sscanf(argv[argcount+1], "%lf", &E_MAX_6);
		else if(argv[argcount][3]=='7') sscanf(argv[argcount+1], "%lf", &E_MAX_7);
		else if(argv[argcount][3]=='8') sscanf(argv[argcount+1], "%lf", &E_MAX_8);
		else if(argv[argcount][3]=='9') sscanf(argv[argcount+1], "%lf", &E_MAX_9);
	      }

	      else {
		  printf(" Error at reading input arguments:\n");
		  printf("-ey -eP0 ... -eP16 -en -eE -eV -em0 ... -em16 -eM0 -eM16 are the only codes allowed. You have typed %s", argv[argcount]);
		  exit(0);
	      }

        skip++;
        break;
