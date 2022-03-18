/* Human Parameters */

       case 'p': /* Model Parameters */
         if(argv[argcount][2]=='K')        sscanf(argv[argcount+1],"%lf", &K_0);     /* K_0  */

	 else if(argv[argcount][2]=='U')   sscanf(argv[argcount+1],"%lf", &Nu_C);    /* Nu_C */

	 else if(argv[argcount][2]=='A')   sscanf(argv[argcount+1],"%lf", &Alpha_C); /* Alpha_C */

	 else if(argv[argcount][2]=='B') {

	   if(argv[argcount][3]=='\0') {
	     printf("Error in Model Parameter Command Line Argument: -pB does not exist!!!\n");
	     printf("Try -pBR [VALUE] or -pBC [VALUE]. The program will exit\n");
	     exit(0);
	   }
	   
	   else if(argv[argcount][3]=='R') sscanf(argv[argcount+1],"%lf", &Beta_R); /* Beta_R */
	   
	   else if(argv[argcount][3]=='C') sscanf(argv[argcount+1],"%lf", &Beta_C); /* Beta_C */ 
	    
	   else {
	     printf(" Error in Model Parameter Command Line Argument: %s does not exist!!!\n",
		    argv[argcount]);
	     printf(" in include.Parameter_Model.argumentControl.c\n");
	     printf(" An example for a command line program call is:\n");
	     printf(" ~$ %s -pK [VALUE] -pBC [VALUE] -pBR [VALUE] -pU [VALUE] -pA [VALUE] -pDC [VALUE] -pDR [VALUE]\n", argv[0]);
	     printf(" The program will exit\n");
	     exit(0);
	   }
	 }

	 else if(argv[argcount][2]=='D') {

	   if(argv[argcount][3]=='\0') {
	     printf(" Error in Model Parameter Command Line Argument: -pD does not exist!!!\n");
	     printf(" Try -pDR [VALUE] or -pDC [VALUE]. The program will exit\n");
	     exit(0);
	   }
	   
	   else if(argv[argcount][3]=='R') sscanf(argv[argcount+1],"%lf", &Delta_R); /* Delta_R */
	   
	   else if(argv[argcount][3]=='C') sscanf(argv[argcount+1],"%lf", &Delta_C); /* Delta_C */
	    
	   else {
	     printf(" Error in Model Parameter Command Line Argument: %s does not exist!!!\n",
		    argv[argcount]);
	     printf(" in include.Parameter_Model.argumentControl.c\n");
	     printf(" An example for a command line program call is:\n");
	     printf(" ~$ %s -pK [VALUE] -pBC [VALUE] -pBR [VALUE] -pU [VALUE] -pA [VALUE] -pDC [VALUE] -pDR [VALUE]\n", argv[0]);
	     printf(" The program will exit\n");
	     exit(0);
	   }
	 }

         else {
	   printf(" Error in Model Parameter Command Line Argument: %s does not exist!!!\n",
		    argv[argcount]);
	   printf(" in include.Parameter_Model.argumentControl.c\n");
	   printf(" An example for a command line program call is:\n");
	   printf(" ~$ %s -pK [VALUE] -pBC [VALUE] -pBR [VALUE] -pU [VALUE] -pA [VALUE] -pDC [VALUE] -pDR [VALUE]\n", argv[0]);
	   printf(" The program will exit\n");
	   exit(0);
	   
	 } 
         skip++;
         break;
