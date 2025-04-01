#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <MODEL.h>

// Function prototypes
void cleanup(Parameter_Table * Table);
void handle_signal(int signal);

/*  To catch segmentation faults and ensure that all memory is properly freed at the end of the program, 
    you can use signal handling in C. The signal function can be used to catch segmentation faults and 
    other signals. When a segmentation fault occurs, you can call a cleanup function to free all allocated 
    memory before exiting the program.

    Here is an example of how you can implement this in your program:

    Include the necessary headers for signal handling.
    Define a cleanup function to free all allocated memory.
    Set up a signal handler for segmentation faults.
*/
/*
    The main function sets up a signal handler for segmentation faults and calls the program logic. 
    At the end of the program, it calls the cleanup function to free all allocated memory.

    The handle_signal function is called when a segmentatioSn fault occurs. It prints a message,
    calls the cleanup function, and exits the program with an error code.

    The cleanup function frees all dynamically allocated memory in the program. You should add
    additional memory freeing code as needed for your specific program.

    This approach ensures that all memory is properly freed when a segmentation fault occurs and 
    when the program exits normally
*/

int main(int argc, char *argv[]) {
    // Set up signal handler for segmentation faults
    signal(SIGSEGV, handle_signal);

    // Your program logic here

    // At the end of the program, free all allocated memory
    cleanup(Table);

    return 0;
}

void handle_signal(int signal) {
    if (signal == SIGSEGV) {
        printf("Segmentation fault caught!\n");
        cleanup(Table);
        exit(1);
    }
}

void cleanup(Parameter_Table * Table) {
    // Free all allocated memory here
    // Example:
    if (Table->Delta_AP) free(Table->Delta_AP);
    if (Table->Beta_AP) free(Table->Beta_AP);
    if (Table->Vector_Model_Variables_Stationarity) free(Table->Vector_Model_Variables_Stationarity);
    // Add other memory freeing code as needed

    // Free other dynamically allocated memory in your program
}