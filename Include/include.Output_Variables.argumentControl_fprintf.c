/* Code indetifying which variable is plotted */

              case 'v': /* * * * * * * * */
	      if(argv[argcount][2]=='0')
		fprintf(fp,"-v0  %d  ",variable_0);
	      else if(argv[argcount][2]=='1'){
		if(argv[argcount][3]=='0'){
		  fprintf(fp,"-v10  %d  ",variable_10);
		}
		else if(argv[argcount][3]=='1'){
		  fprintf(fp,"-v11  %d  ",variable_11);
		}
		else if(argv[argcount][3]=='2'){
		  fprintf(fp,"-v12  %d  ",variable_12);
		}
		else if(argv[argcount][3]=='3'){
		  fprintf(fp,"-v13  %d  ",variable_13);
		}
		else if(argv[argcount][3]=='4'){
		  fprintf(fp,"-v14  %d  ",variable_14);
		}
		else if(argv[argcount][3]=='5'){
		  fprintf(fp,"-v15  %d  ",variable_15);
		}
                else if(argv[argcount][3]=='6'){
                  fprintf(fp,"-v16  %d  ",variable_16);
                }
		else if(argv[argcount][3]=='7'){
                  fprintf(fp,"-v17  %d  ",variable_17);
                }
		else if(argv[argcount][3]=='8'){
                  fprintf(fp,"-v18  %d  ",variable_18);
                }
		else if(argv[argcount][3]=='9'){
                  fprintf(fp,"-v19  %d  ",variable_19);
                }
		else{
		  fprintf(fp,"-v1  %d  ", variable_1);
		}
	      }
	      else if(argv[argcount][2]=='2'){
	        if(argv[argcount][3]=='0'){
                  fprintf(fp,"-v20  %d  ",variable_20);
                }
                else if(argv[argcount][3]=='1'){
                  fprintf(fp,"-v21  %d  ",variable_21);
                }
		else if(argv[argcount][3]=='2'){
                  fprintf(fp,"-v22  %d  ",variable_22);
                }
		else if(argv[argcount][3]=='3'){
                  fprintf(fp,"-v23  %d  ",variable_23);
                }
                else if(argv[argcount][3]=='4'){
                  fprintf(fp,"-v24  %d  ",variable_24);
                }
                else if(argv[argcount][3]=='5'){
                  fprintf(fp,"-v25  %d  ",variable_25);
                }
                else if(argv[argcount][3]=='6'){
                  fprintf(fp,"-v26  %d  ",variable_26);
                }
		else if(argv[argcount][3]=='7'){
                  fprintf(fp,"-v27  %d  ",variable_27);
                }
                else if(argv[argcount][3]=='8'){
                  fprintf(fp,"-v28  %d  ",variable_28);
                }
		else if(argv[argcount][3]=='9'){
                  fprintf(fp,"-29  %d  ",variable_29);
                }
                else{
		  fprintf(fp,"-v2  %d  ",variable_2);
                }
	      }
              else if(argv[argcount][2]=='3'){
                if(argv[argcount][3]=='0'){
                  fprintf(fp,"-v30  %d  ",variable_30);
                }
                else if(argv[argcount][3]=='1'){
                  fprintf(fp,"-v31  %d  ",variable_31);
                }
                else if(argv[argcount][3]=='2'){
                  fprintf(fp,"-v32  %d  ",variable_32);
                }
		else if(argv[argcount][3]=='3'){
                  fprintf(fp,"-v33  %d  ",variable_33);
                }
                else if(argv[argcount][3]=='4'){
                  fprintf(fp,"-v34  %d  ",variable_34);
                }
                else if(argv[argcount][3]=='5'){
                  fprintf(fp,"-v35  %d  ",variable_35);
                }
                else if(argv[argcount][3]=='6'){
                  fprintf(fp,"-v36  %d  ",variable_36);
                }
                else if(argv[argcount][3]=='5'){
                  fprintf(fp,"-v37  %d  ",variable_37);
                }
                else if(argv[argcount][3]=='6'){
                  fprintf(fp,"-v38  %d  ",variable_38);
                }
                else{
                  fprintf(fp,"-v3  %d  ",variable_3);
                }
              }
	      else if(argv[argcount][2]=='4')
		fprintf(fp,"-v4  %d  ",variable_4);
	      else if(argv[argcount][2]=='5')
		fprintf(fp,"-v5  %d  ",variable_5);
	      else if(argv[argcount][2]=='6')
		fprintf(fp,"-v6  %d  ",variable_6);
	      else if(argv[argcount][2]=='7')
		fprintf(fp,"-v7  %d  ",variable_7);
	      else if(argv[argcount][2]=='8')
		fprintf(fp,"-v8  %d  ",variable_8);
              else if(argv[argcount][2]=='9')
		fprintf(fp,"-v9  %d  ",variable_9);
	      else exit(0);
	      skip++;
	      break;

 	      case 'n': /* * * *  */
              fprintf(fp,"-n  %d  ", Table->SUB_OUTPUT_VARIABLES);
              skip++;
              break;
