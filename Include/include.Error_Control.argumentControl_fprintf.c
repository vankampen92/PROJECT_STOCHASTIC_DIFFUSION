/* Scan Parameters */

/* Error Model Parameters */
	case 'e':
              if(argv[argcount][2]=='0') 	 fprintf(fp, "-e0  %g  ",Table->Err_0);

	      else if(argv[argcount][2]=='1') {
		if(argv[argcount][3]=='0') fprintf(fp, "-e10  %g  ", Table->Err_10);
		else if(argv[argcount][3]=='1')  fprintf(fp, "-e11  %g  ", Table->Err_11);
                else if(argv[argcount][3]=='2')  fprintf(fp, "-e12  %g  ", Table->Err_12);
		else if(argv[argcount][3]=='3')  fprintf(fp, "-e13  %g  ", Table->Err_13);
		else if(argv[argcount][3]=='\0') fprintf(fp, "-e1  %g ", Table->Err_1);
	      }

	      else if(argv[argcount][2]=='2')    fprintf(fp, "-e2  %g  ",Table->Err_2);
	      else if(argv[argcount][2]=='3')    fprintf(fp, "-e3  %g  ",Table->Err_3);
              else if(argv[argcount][2]=='4')    fprintf(fp, "-e4  %g  ",Table->Err_4);
              else if(argv[argcount][2]=='5')    fprintf(fp, "-e5  %g  ",Table->Err_5);
	      else if(argv[argcount][2]=='6')    fprintf(fp, "-e6  %g  ",Table->Err_6);
	      else if(argv[argcount][2]=='7')    fprintf(fp, "-e7  %g  ",Table->Err_7);
              else if(argv[argcount][2]=='8')    fprintf(fp, "-e8  %g  ",Table->Err_8);
              else if(argv[argcount][2]=='9')    fprintf(fp, "-e9  %g  ",Table->Err_9);

	      else if(argv[argcount][2]=='y')
		fprintf(fp, "-ey  %d  ",Table->TYPE_of_ERROR_FUNCTION);

              else if(argv[argcount][2]=='n')
		fprintf(fp, "-en  %d  ",Table->No_of_ERROR_PARAMETERS);

              else if(argv[argcount][2]=='E')
		fprintf(fp, "-eE  %g  ",Table->Error_Parameter_1);

              else if(argv[argcount][2]=='V')
		fprintf(fp, "-eV  %g  ",Table->Error_Parameter_0);

              else if(argv[argcount][2]=='P') {
		if(argv[argcount][3]=='0')      fprintf(fp, "-eP0  %d  ", Table->E_Var_0);
		else if(argv[argcount][3]=='1') {
		  if(argv[argcount][4]=='0')      fprintf(fp, "-eP10  %d  ", Table->E_Var_10);
		  else if(argv[argcount][4]=='1') fprintf(fp, "-eP11  %d  ", Table->E_Var_11);
		  else if(argv[argcount][4]=='2') fprintf(fp, "-eP12  %d  ", Table->E_Var_12);
		  else if(argv[argcount][4]=='3') fprintf(fp, "-eP13  %d  ", Table->E_Var_13);
		  else if(argv[argcount][4]=='4') fprintf(fp, "-eP14  %d  ", Table->E_Var_14);
		  else if(argv[argcount][4]=='5') fprintf(fp, "-eP15  %d  ", Table->E_Var_15);
		  else if(argv[argcount][4]=='6') fprintf(fp, "-eP16  %d  ", Table->E_Var_16);
		  else if(argv[argcount][4]=='\0') fprintf(fp, "-eP1  %d  ", Table->E_Var_1);
		}
                else if(argv[argcount][3]=='2') fprintf(fp, "-eP2  %d  ", Table->E_Var_2);
		else if(argv[argcount][3]=='3') fprintf(fp, "-eP3  %d  ", Table->E_Var_3);
		else if(argv[argcount][3]=='4') fprintf(fp, "-eP4  %d  ", Table->E_Var_4);
		else if(argv[argcount][3]=='5') fprintf(fp, "-eP5  %d  ", Table->E_Var_5);
		else if(argv[argcount][3]=='6') fprintf(fp, "-eP6  %d  ", Table->E_Var_6);
		else if(argv[argcount][3]=='7') fprintf(fp, "-eP7  %d  ", Table->E_Var_7);
		else if(argv[argcount][3]=='8') fprintf(fp, "-eP8  %d  ", Table->E_Var_8);
		else if(argv[argcount][3]=='9') fprintf(fp, "-eP9  %d  ", Table->E_Var_9);
	      }

	      else if(argv[argcount][2]=='m') {
		if(argv[argcount][3]=='0')      fprintf(fp, "-em0  %g  ", Table->E_min_0);
		else if(argv[argcount][3]=='1') {
		  if(argv[argcount][4]=='0')      fprintf(fp, "-em10  %g  ", Table->E_min_10);
		  else if(argv[argcount][4]=='1') fprintf(fp, "-em11  %g  ", Table->E_min_11);
		  else if(argv[argcount][4]=='2') fprintf(fp, "-em12  %g  ", Table->E_min_12);
		  else if(argv[argcount][4]=='3') fprintf(fp, "-em13  %g  ", Table->E_min_13);
		  else if(argv[argcount][4]=='4') fprintf(fp, "-em14  %g  ", Table->E_min_14);
		  else if(argv[argcount][4]=='5') fprintf(fp, "-em15  %g  ", Table->E_min_15);
		  else if(argv[argcount][4]=='6') fprintf(fp, "-em16  %g  ", Table->E_min_16);
		  else if(argv[argcount][4]=='\0') fprintf(fp, "-em1  %g  ", Table->E_min_1);
		}
		else if(argv[argcount][3]=='2') fprintf(fp, "-em2  %g  ", Table->E_min_2);
		else if(argv[argcount][3]=='3') fprintf(fp, "-em3  %g  ", Table->E_min_3);
		else if(argv[argcount][3]=='4') fprintf(fp, "-em4  %g  ", Table->E_min_4);
		else if(argv[argcount][3]=='5') fprintf(fp, "-em5  %g  ", Table->E_min_5);
		else if(argv[argcount][3]=='6') fprintf(fp, "-em6  %g  ", Table->E_min_6);
		else if(argv[argcount][3]=='7') fprintf(fp, "-em7  %g  ", Table->E_min_7);
		else if(argv[argcount][3]=='8') fprintf(fp, "-em8  %g  ", Table->E_min_8);
		else if(argv[argcount][3]=='9') fprintf(fp, "-em9  %g  ", Table->E_min_9);
	      }

	      else if(argv[argcount][2]=='M') {
		if(argv[argcount][3]=='0')      fprintf(fp, "-eM0  %g  ", Table->E_MAX_0);
		else if(argv[argcount][3]=='1') {
		  if(argv[argcount][4]=='0')      fprintf(fp, "-eM10  %g  ", Table->E_MAX_10);
		  else if(argv[argcount][4]=='1') fprintf(fp, "-eM11  %g  ", Table->E_MAX_11);
		  else if(argv[argcount][4]=='2') fprintf(fp, "-eM12  %g  ", Table->E_MAX_12);
		  else if(argv[argcount][4]=='3') fprintf(fp, "-eM13  %g  ", Table->E_MAX_13);
		  else if(argv[argcount][4]=='4') fprintf(fp, "-eM14  %g  ", Table->E_MAX_14);
		  else if(argv[argcount][4]=='5') fprintf(fp, "-eM15  %g  ", Table->E_MAX_15);
		  else if(argv[argcount][4]=='6') fprintf(fp, "-eM16  %g  ", Table->E_MAX_16);
		  else if(argv[argcount][4]=='\0') fprintf(fp, "-eM1  %g  ", Table->E_MAX_1);
		}
		else if(argv[argcount][3]=='2') fprintf(fp, "-eM2  %g  ", Table->E_MAX_2);
		else if(argv[argcount][3]=='3') fprintf(fp, "-eM3  %g  ", Table->E_MAX_3);
		else if(argv[argcount][3]=='4') fprintf(fp, "-eM4  %g  ", Table->E_MAX_4);
		else if(argv[argcount][3]=='5') fprintf(fp, "-eM5  %g  ", Table->E_MAX_5);
		else if(argv[argcount][3]=='6') fprintf(fp, "-eM6  %g  ", Table->E_MAX_6);
		else if(argv[argcount][3]=='7') fprintf(fp, "-eM7  %g  ", Table->E_MAX_7);
		else if(argv[argcount][3]=='8') fprintf(fp, "-eM8  %g  ", Table->E_MAX_8);
		else if(argv[argcount][3]=='9') fprintf(fp, "-eM9  %g  ", Table->E_MAX_9);
	      }
	      else {
		printf(" Error at reading input arguments:\n");
								printf("-ey -eP0 ... -eP13 -en -eE -eV -em0 ... -em9 -eM0 -eM9 are the only codes allowed. You have typed -  %s  ", argv[argcount]);
								exit(0);
	      }

        skip++;
        break;
