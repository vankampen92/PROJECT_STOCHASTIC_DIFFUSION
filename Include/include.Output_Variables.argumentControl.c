/* Code indetifying which variable is plotted */

              case 'v': /* * * * * * * * */
	      if(argv[argcount][2]=='0')
		sscanf(argv[argcount+1],"%d",&variable_0);
	      else if(argv[argcount][2]=='1'){
		if(argv[argcount][3]=='0'){
		  sscanf(argv[argcount+1],"%d",&variable_10);
		}
		else if(argv[argcount][3]=='1'){
		  sscanf(argv[argcount+1],"%d",&variable_11);
		}
		else if(argv[argcount][3]=='2'){
		  sscanf(argv[argcount+1],"%d",&variable_12);
		}
		else if(argv[argcount][3]=='3'){
		  sscanf(argv[argcount+1],"%d",&variable_13);
		}
		else if(argv[argcount][3]=='4'){
		  sscanf(argv[argcount+1],"%d",&variable_14);
		}
		else if(argv[argcount][3]=='5'){
		  sscanf(argv[argcount+1],"%d",&variable_15);
		}
                else if(argv[argcount][3]=='6'){
                  sscanf(argv[argcount+1],"%d",&variable_16);
                }
		else if(argv[argcount][3]=='7'){
                  sscanf(argv[argcount+1],"%d",&variable_17);
                }
		else if(argv[argcount][3]=='8'){
                  sscanf(argv[argcount+1],"%d",&variable_18);
                }
		else if(argv[argcount][3]=='9'){
                  sscanf(argv[argcount+1],"%d",&variable_19);
                }
		else{
		  sscanf(argv[argcount+1],"%d", &variable_1);
		}
	      }
	      else if(argv[argcount][2]=='2'){
	        if(argv[argcount][3]=='0'){
                  sscanf(argv[argcount+1],"%d",&variable_20);
                }
                else if(argv[argcount][3]=='1'){
                  sscanf(argv[argcount+1],"%d",&variable_21);
                }
		else if(argv[argcount][3]=='2'){
                  sscanf(argv[argcount+1],"%d",&variable_22);
                }
		else if(argv[argcount][3]=='3'){
                  sscanf(argv[argcount+1],"%d",&variable_23);
                }
                else if(argv[argcount][3]=='4'){
                  sscanf(argv[argcount+1],"%d",&variable_24);
                }
                else if(argv[argcount][3]=='5'){
                  sscanf(argv[argcount+1],"%d",&variable_25);
                }
                else if(argv[argcount][3]=='6'){
                  sscanf(argv[argcount+1],"%d",&variable_26);
                }
		else if(argv[argcount][3]=='7'){
                  sscanf(argv[argcount+1],"%d",&variable_27);
                }
                else if(argv[argcount][3]=='8'){
                  sscanf(argv[argcount+1],"%d",&variable_28);
                }
		else if(argv[argcount][3]=='9'){
                  sscanf(argv[argcount+1],"%d",&variable_29);
                }
                else{
		  sscanf(argv[argcount+1],"%d",&variable_2);
                }
	      }
              else if(argv[argcount][2]=='3'){
                if(argv[argcount][3]=='0'){
                  sscanf(argv[argcount+1],"%d",&variable_30);
                }
                else if(argv[argcount][3]=='1'){
                  sscanf(argv[argcount+1],"%d",&variable_31);
                }
                else if(argv[argcount][3]=='2'){
                  sscanf(argv[argcount+1],"%d",&variable_32);
                }
		else if(argv[argcount][3]=='3'){
                  sscanf(argv[argcount+1],"%d",&variable_33);
                }
                else if(argv[argcount][3]=='4'){
                  sscanf(argv[argcount+1],"%d",&variable_34);
                }
                else if(argv[argcount][3]=='5'){
                  sscanf(argv[argcount+1],"%d",&variable_35);
                }
                else if(argv[argcount][3]=='6'){
                  sscanf(argv[argcount+1],"%d",&variable_36);
                }
                else if(argv[argcount][3]=='7'){
                  sscanf(argv[argcount+1],"%d",&variable_37);
                }
                else if(argv[argcount][3]=='8'){
                  sscanf(argv[argcount+1],"%d",&variable_38);
                }
                else{
                  sscanf(argv[argcount+1],"%d",&variable_3);
                }
              }
	      else if(argv[argcount][2]=='4')
		sscanf(argv[argcount+1],"%d",&variable_4);
	      else if(argv[argcount][2]=='5')
		sscanf(argv[argcount+1],"%d",&variable_5);
	      else if(argv[argcount][2]=='6')
		sscanf(argv[argcount+1],"%d",&variable_6);
	      else if(argv[argcount][2]=='7')
		sscanf(argv[argcount+1],"%d",&variable_7);
	      else if(argv[argcount][2]=='8')
		sscanf(argv[argcount+1],"%d",&variable_8);
              else if(argv[argcount][2]=='9')
		sscanf(argv[argcount+1],"%d",&variable_9);
	      else exit(0);
	      skip++;
	      break;

 	      case 'n': /* * * *  */
              sscanf(argv[argcount+1],"%d",&SUB_OUTPUT_VARIABLES);
              skip++;
              break;
