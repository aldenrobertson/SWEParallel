#include <math.h>
#include <iostream>
#include <mpi.h>
#include <netcdf>
#include <vector>
using namespace netCDF;
using namespace netCDF::exceptions;

using namespace std;

/*Lets parallise our script!*/
/*MUst also think about the indexing. Am I doing rows or columns etc*/

int main (int argc, char **argv) {
	int gridsize = 190;
	double dt = 0.0001;
	double dx = 1;
	double dy = 1;
	double g = 9.81;
	double H = 0;
	double f = 0.000113; /*Coriolis parameter at London apparently*/
	double gx, gy;
	gx = g / (2 * dx);
	gy = g / (2 * dy);
	int x, j, i, y;
	double h_grid[gridsize + 2][gridsize + 2]; /* values of h and u which i will find*/
	double u_grid[gridsize + 2][gridsize + 2];
	double v_grid[gridsize + 2][gridsize + 2]; /*added the plus 2 for the nodes at the other end of the time period. This due to the finite difference scheme I used.*/
	double h_grid_old[gridsize + 2][gridsize + 2]; /* values of h and u which i will find*/
	double u_grid_old[gridsize + 2][gridsize + 2];
	double v_grid_old[gridsize + 2][gridsize + 2]; /*added the plus 2 for the nodes at the other end of the time period. This due to the finite difference scheme I used.*/
	double h_grid_older[gridsize + 2][gridsize + 2]; /* values of h and u which i will find*/
	double u_grid_older[gridsize + 2][gridsize + 2];
	double v_grid_older[gridsize + 2][gridsize + 2];


	int size, rank;    /* The size and rank which I will sort out when I start with the MPI.*/

	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int jParallel = (gridsize + 2) / size;     /* This is the vertical strips. Need to think about if we are row or column major etc*/
	double masterbuf[gridsize+2][gridsize+2];
	double masterbufU[gridsize+2][gridsize+2];
	double masterbufV[gridsize + 2][gridsize + 2];
	double masterbufH[gridsize + 2][gridsize + 2];
	double masterbufUOLD[gridsize + 2][gridsize + 2];
	double masterbufVOLD[gridsize + 2][gridsize + 2];
	double masterbufHOLD[gridsize + 2][gridsize + 2];
	double masterbufUOLDER[gridsize + 2][gridsize + 2];
	double masterbufVOLDER[gridsize + 2][gridsize + 2];
	double masterbufHOLDER[gridsize + 2][gridsize + 2];
	double bufU[jParallel][gridsize+2];        /*Sub gridies*/
	double bufV[jParallel][gridsize + 2];        /*Sub gridies*/
	double bufH[jParallel][gridsize + 2];        /*Sub gridies*/
	double bufUOLD[jParallel][gridsize + 2];        /*Sub gridies*/
	double bufVOLD[jParallel][gridsize + 2];        /*Sub gridies*/
	double bufHOLD[jParallel][gridsize + 2];        /*Sub gridies*/
	double bufUOLDER[jParallel][gridsize + 2];        /*Sub gridies*/
	double bufVOLDER[jParallel][gridsize + 2];        /*Sub gridies*/
	double bufHOLDER[jParallel][gridsize + 2];        /*Sub gridies*/
	double topRowU[gridsize+2];
	double bottomRowU[gridsize+2];
	double topRowV[gridsize + 2];
	double bottomRowV[gridsize + 2];
	double topRowH[gridsize + 2];
	double bottomRowH[gridsize + 2];
	double num = 0;
	if (rank == 0){                                   /*SORT THIS OUT*/
		for (i = 0; i <= gridsize + 1; i++){
			for (j = 0; j <= gridsize + 1; j++){
				masterbuf[i][j] = num;
				num = num + 1;
				cout << masterbuf[i][j];
				cout << " ";
			}
			cout << "\n";
		}
	}

	for (j = 0; j <= gridsize + 1; j++) {
		for (i = 0; i <= gridsize + 1;i++) {
			masterbufU[i][j] = 0;
			masterbufV[i][j] = 0;
			masterbufH[i][j] = 1;
			masterbufUOLD[i][j] = 0;
			masterbufVOLD[i][j] = 0;
			masterbufHOLD[i][j] = 1;
			masterbufUOLDER[i][j] = 0;
			masterbufVOLDER[i][j] = 0;
			masterbufHOLDER[i][j] = 1;
		}
	}
	for (i = 0; i <= gridsize + 1; i++)               /*Initial conditions of h,u,v*/
	{
		for (j = 0; j <= gridsize + 1; j++)
		{
			masterbufH[i][j]= 0.8 + 0.5 * exp(-((i - 30)*(i-30) + (j - 30)*(j-30)) / 25.0);
			masterbufHOLD[i][j] = 0.8 + 0.5 * exp(-((i - 30)*(i-30) + (j - 30)*(j-30)) / 25.0);
			masterbufHOLDER[i][j] = 0.8 + 0.5 * exp(-((i - 30)*(i - 30) + (j-30)* (j - 30)) / 25.0);
		}
	}
	double t1, t2; 
	t1 = MPI_Wtime();
        if (rank == 0){                                   /*SORT THIS OUT*/
                for (i = 0; i <= gridsize + 1; i++){
                        for (j = 0; j <= gridsize + 1; j++){
                                cout << masterbufHOLDER[i][j];
                                cout << " ";
                        }
                        cout << "\n";
                }
        }

	cout << "\n";

	if (rank == 0){                                   /*SORT THIS OUT*/
                for (i = 0; i <= gridsize + 1; i++){
                        for (j = 0; j <= gridsize + 1; j++){
                                cout << masterbufUOLD[i][j];
                                cout << " ";
                        }
                        cout << "\n";
                }
        }
	for (i = 0; i <= gridsize + 1; i++){
		for (j = 0; j <= gridsize + 1; j ++){
			H = H + masterbufH[i][j];
		}
	}
	H = H/((gridsize+2)*(gridsize+2));
	cout << H << "\n";

        cout << "\n";
	/*For each timestep I will scatter the (i,j) grid across the processors in columns then each processor will work on a column
	with halo swapping and whatnot, then the master process (rank 0) will gather these to form the (i,j) grid at time n+1 and repeat.*/
        MPI_Scatter(&masterbufU[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &bufU[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&masterbufV[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &bufV[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&masterbufH[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &bufH[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&masterbufUOLD[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &bufUOLD[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&masterbufVOLD[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &bufVOLD[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&masterbufHOLD[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &bufHOLD[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        MPI_Scatter(&masterbufUOLDER[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &bufUOLDER[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&masterbufVOLDER[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &bufVOLDER[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&masterbufHOLDER[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &bufHOLDER[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (x = 1;x <= 3200000;x++) {            /*2nd order euler method*/

		if (rank == 0) {

			for (j = 0; j <= gridsize + 1; j++) {
				bottomRowU[j] = bufUOLD[jParallel - 1][j];
				bottomRowV[j] = bufVOLD[jParallel - 1][j];
				bottomRowH[j] = bufHOLD[jParallel - 1][j];
                                topRowU[j] = bufUOLD[0][j];
                                topRowV[j] = bufVOLD[0][j];
                                topRowH[j] = bufHOLD[0][j];
			}

			MPI_Send(bottomRowU, gridsize + 2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(bottomRowU, gridsize + 2, MPI_DOUBLE, size-1, size-1, MPI_COMM_WORLD, &status);
			MPI_Send(topRowU, gridsize + 2, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD);
			MPI_Recv(topRowU, gridsize + 2, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

			MPI_Send(bottomRowV, gridsize + 2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(bottomRowV, gridsize + 2, MPI_DOUBLE, size-1, size-1, MPI_COMM_WORLD, &status);
			MPI_Send(topRowV, gridsize + 2, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD);
			MPI_Recv(topRowV, gridsize + 2, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);


			MPI_Send(bottomRowH, gridsize + 2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(bottomRowH, gridsize + 2, MPI_DOUBLE, size-1, size-1, MPI_COMM_WORLD, &status);
                        MPI_Send(topRowH, gridsize + 2, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD);
			MPI_Recv(topRowH, gridsize + 2, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);


			double concBufUOLD[jParallel + 2][gridsize + 2];
			double concBufVOLD[jParallel + 2][gridsize + 2];
			double concBufHOLD[jParallel + 2][gridsize + 2];
			for (j = 0; j <= gridsize + 1; j++) {
				concBufUOLD[jParallel+1][j] = topRowU[j];
				concBufVOLD[jParallel+1][j] = topRowV[j];
				concBufHOLD[jParallel+1][j] = topRowH[j];
				concBufUOLD[0][j] = bottomRowU[j];
				concBufVOLD[0][j] = bottomRowV[j];
				concBufHOLD[0][j] = bottomRowH[j];
			}
			for (j = 0; j <= gridsize + 1;j++) {
				for (i = 1; i <= jParallel; i++) {
					concBufUOLD[i][j] = bufUOLD[i-1][j];
					concBufVOLD[i][j] = bufVOLD[i-1][j];
					concBufHOLD[i][j] = bufHOLD[i-1][j];
				}
			}
			/*CONCBUFOLD HOLDS THE VALUES FROM THE PREVIOUS TIME STEP*/

			for (i = 1;i <= jParallel;i++) {
				for (j = 0; j <= gridsize+1; j++) {
					if (j == 0){
						bufU[i-1][j] = (4 * concBufUOLD[i][j] - bufUOLDER[i-1][j] + 2 * dt * (-gx * (concBufHOLD[i + 1][j] - concBufHOLD[i - 1][j]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufUOLD[i][j + 1] - concBufUOLD[i][gridsize+1]) + f * concBufVOLD[i][j])) / 3;
                                        	bufV[i-1][j] = (4 * concBufVOLD[i][j] - bufVOLDER[i-1][j] + 2 * dt * (-gy * (concBufHOLD[i][j + 1] - concBufHOLD[i][gridsize+1]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufVOLD[i + 1][j] - concBufVOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufVOLD[i][j + 1] - concBufVOLD[i][gridsize+1]) - f * concBufUOLD[i][j])) / 3;
                                        	bufH[i-1][j] = (4 * concBufHOLD[i][j] - bufHOLDER[i-1][j] -2 * H * dt * ((concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) / (2 * dx) + (concBufVOLD[i][j + 1] - concBufVOLD[i][gridsize + 1]) / (2 * dy))) / 3;
					}
					else{
						bufU[i-1][j] = (4 * concBufUOLD[i][j] - bufUOLDER[i-1][j] + 2 * dt * (-gx * (concBufHOLD[i + 1][j] - concBufHOLD[i - 1][j]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufUOLD[i][(j + 1)% (gridsize + 2)] - concBufUOLD[i][j - 1]) + f * concBufVOLD[i][j])) / 3;
						bufV[i-1][j] = (4 * concBufVOLD[i][j] - bufVOLDER[i-1][j] + 2 * dt * (-gy * (concBufHOLD[i][(j + 1)% (gridsize + 2)] - concBufHOLD[i][j - 1]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufVOLD[i + 1][j] - concBufVOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufVOLD[i][(j + 1)% (gridsize + 2)] - concBufVOLD[i][j - 1]) - f * concBufUOLD[i][j])) / 3;
						bufH[i-1][j] = (4 * concBufHOLD[i][j] - bufHOLDER[i-1][j] -2 * H * dt * ((concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) / (2 * dx) + (concBufVOLD[i][(j + 1)% (gridsize + 2)] - concBufVOLD[i][j - 1]) / (2 * dy))) / 3;
					}
				}
			}
                for (j = 0;j <= gridsize + 1;j++) {
                        for (i = 0; i <= jParallel - 1; i++) {
                                bufUOLDER[i][j] = bufUOLD[i][j];           /*MUST THINK ABOUT INDEXES*/
                                bufVOLDER[i][j] = bufVOLD[i][j];
                                bufHOLDER[i][j] = bufHOLD[i][j];
                        }
                }

                for (j = 0;j <= gridsize+1; j++) {
                        for (i = 0; i <= jParallel - 1; i++) {
                                bufUOLD[i][j] = bufU[i][j];
                                bufVOLD[i][j] = bufV[i][j];
                                bufHOLD[i][j] = bufH[i][j];                                                                                                                                                          }
                }

		}



		if (rank > 0 && rank < size - 1) {
			for (j = 0; j <= gridsize + 1; j++) {
				bottomRowU[j] = bufUOLD[jParallel - 1][j];
				bottomRowV[j] = bufVOLD[jParallel - 1][j];
				bottomRowH[j] = bufHOLD[jParallel - 1][j];
				topRowU[j] = bufUOLD[0][j];
				topRowV[j] = bufVOLD[0][j];
				topRowH[j] = bufHOLD[0][j];
			}

			/*SEND THE ROWS OF U AND V GRIDS BETWEEN OUR PROCESSES*/

			MPI_Send(bottomRowU, gridsize + 2, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
			MPI_Recv(bottomRowU, gridsize + 2, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
			MPI_Send(topRowU, gridsize + 2, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
			MPI_Recv(topRowU, gridsize + 2, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status);

			MPI_Send(bottomRowV, gridsize + 2, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
			MPI_Recv(bottomRowV, gridsize + 2, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
			MPI_Send(topRowV, gridsize + 2, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
			MPI_Recv(topRowV, gridsize + 2, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status);

			MPI_Send(bottomRowH, gridsize + 2, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
			MPI_Recv(bottomRowH, gridsize + 2, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
			MPI_Send(topRowH, gridsize + 2, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
			MPI_Recv(topRowH, gridsize + 2, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status);
			
			/*cout << "send and recieved on middle ranks \n";*/
			double concBufUOLD[jParallel + 2][gridsize + 2];
			double concBufVOLD[jParallel + 2][gridsize + 2];
			double concBufHOLD[jParallel + 2][gridsize + 2];

                        for (j = 0; j <= gridsize + 1; j++) {
                                concBufUOLD[jParallel+1][j] = topRowU[j];
                                concBufVOLD[jParallel+1][j] = topRowV[j];
                                concBufHOLD[jParallel+1][j] = topRowH[j];
                                concBufUOLD[0][j] = bottomRowU[j];
                                concBufVOLD[0][j] = bottomRowV[j];
                                concBufHOLD[0][j] = bottomRowH[j];
                        }

			for (j = 0; j <= gridsize + 1;j++) {
				for (i = 1; i <= jParallel; i++) {
					concBufUOLD[i][j] = bufUOLD[i - 1][j];
					concBufVOLD[i][j] = bufVOLD[i - 1][j];
					concBufHOLD[i][j] = bufHOLD[i - 1][j];
				}
			}

			for (i = 1;i <= jParallel;i++) {
				for (j = 0; j <= gridsize + 1; j++) {
					if (j == 0){
						bufU[i-1][j] = (4 * concBufUOLD[i][j] - bufUOLDER[i-1][j] + 2 * dt * (-gx * (concBufHOLD[i + 1][j] - concBufHOLD[i - 1][j]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufUOLD[i][j + 1] - concBufUOLD[i][gridsize+1]) + f * concBufVOLD[i][j])) / 3;
						bufV[i-1][j] = (4 * concBufVOLD[i][j] - bufVOLDER[i-1][j] + 2 * dt * (-gy * (concBufHOLD[i][j + 1] - concBufHOLD[i][gridsize+1]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufVOLD[i + 1][j] - concBufVOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufVOLD[i][j + 1] - concBufVOLD[i][gridsize+1]) - f * concBufUOLD[i][j])) / 3;
						bufH[i-1][j] = (4 * concBufHOLD[i][j] - bufHOLDER[i-1][j] - 2 * H * dt * ((concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) / (2 * dx) + (concBufVOLD[i][j + 1] - concBufVOLD[i][gridsize+1]) / (2 * dy))) / 3;
					}
					else{
                                                bufU[i-1][j] = (4 * concBufUOLD[i][j] - bufUOLDER[i-1][j] + 2 * dt * (-gx * (concBufHOLD[i + 1][j] - concBufHOLD[i - 1][j]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufUOLD[i][(j + 1)% (gridsize + 2)] - concBufUOLD[i][j-1]) + f * concBufVOLD[i][j])) / 3;                                                        
						bufV[i-1][j] = (4 * concBufVOLD[i][j] - bufVOLDER[i-1][j] + 2 * dt * (-gy * (concBufHOLD[i][(j + 1)% (gridsize + 2)] - concBufHOLD[i][j-1]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufVOLD[i + 1][j] - concBufVOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufVOLD[i][(j + 1)% (gridsize + 2)] - concBufVOLD[i][j-1]) - f * concBufUOLD[i][j])) / 3;                                                          
						bufH[i-1][j] = (4 * concBufHOLD[i][j] - bufHOLDER[i-1][j] - 2 * H * dt * ((concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) / (2 * dx) + (concBufVOLD[i][(j + 1)% (gridsize + 2)] - concBufVOLD[i][j-1]) / (2 * dy))) / 3;   
					}
				}
			}
                for (j = 0;j <= gridsize + 1;j++) {
                        for (i = 0; i <= jParallel - 1; i++) {
                                bufUOLDER[i][j] = bufUOLD[i][j];           /*MUST THINK ABOUT INDEXES*/
                                bufVOLDER[i][j] = bufVOLD[i][j];
                                bufHOLDER[i][j] = bufHOLD[i][j];
                        }
                }

                for (j = 0;j <= gridsize+1; j++) {
                        for (i = 0; i <= jParallel - 1; i++) {
                                bufUOLD[i][j] = bufU[i][j];
                                bufVOLD[i][j] = bufV[i][j];
                                bufHOLD[i][j] = bufH[i][j];                                                                                                                                                          }
                }


		}

		if (rank == size - 1) {
			/*cout << "Running on final process" << "\n";*/
			for (j = 0; j <= gridsize + 1; j++) {
				topRowU[j] = bufUOLD[0][j];
				topRowV[j] = bufVOLD[0][j];
				topRowH[j] = bufHOLD[0][j];
                                bottomRowU[j] = bufUOLD[jParallel - 1][j];
                                bottomRowV[j] = bufVOLD[jParallel - 1][j];
                                bottomRowH[j] = bufHOLD[jParallel - 1][j];
			}

			MPI_Send(topRowU, gridsize + 2, MPI_DOUBLE, size - 2, rank, MPI_COMM_WORLD);
			MPI_Send(bottomRowU, gridsize + 2, MPI_DOUBLE, 0, size-1, MPI_COMM_WORLD);
			MPI_Recv(topRowU, gridsize + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);         			
			MPI_Recv(bottomRowU, gridsize + 2, MPI_DOUBLE, size - 2, rank - 1, MPI_COMM_WORLD, &status);

			MPI_Send(topRowV, gridsize + 2, MPI_DOUBLE, size - 2, rank, MPI_COMM_WORLD);
                        MPI_Send(bottomRowV, gridsize + 2, MPI_DOUBLE, 0, size-1, MPI_COMM_WORLD);
                        MPI_Recv(topRowV, gridsize + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(bottomRowV, gridsize + 2, MPI_DOUBLE, size - 2, rank - 1, MPI_COMM_WORLD, &status);

			MPI_Send(topRowH, gridsize + 2, MPI_DOUBLE, size - 2, rank, MPI_COMM_WORLD);
                        MPI_Send(bottomRowH, gridsize + 2, MPI_DOUBLE, 0, size-1, MPI_COMM_WORLD);
                        MPI_Recv(topRowH, gridsize + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); 
			MPI_Recv(bottomRowH, gridsize + 2, MPI_DOUBLE, size - 2, rank - 1, MPI_COMM_WORLD, &status);
 
			double concBufUOLD[jParallel + 2][gridsize + 2];
			double concBufVOLD[jParallel + 2][gridsize + 2];
			double concBufHOLD[jParallel + 2][gridsize + 2];


			for (j = 0; j <= gridsize + 1; j++) {
				concBufUOLD[0][j] = bottomRowU[j];
                                concBufUOLD[jParallel + 1][j] = topRowU[j];
				concBufVOLD[0][j] = bottomRowV[j];
                                concBufVOLD[jParallel + 1][j] = topRowV[j];
				concBufHOLD[0][j] = bottomRowH[j];
                                concBufHOLD[jParallel + 1][j] = topRowH[j];
			}
			for (j = 0; j <= gridsize + 1;j++) {
				for (i = 1; i <= jParallel; i++) {
					concBufUOLD[i][j] = bufUOLD[i - 1][j];
					concBufVOLD[i][j] = bufVOLD[i - 1][j];
					concBufHOLD[i][j] = bufHOLD[i - 1][j];
				}
			}

			/*for (i = 0; i <= jParallel; i++) {
				for (j = 0; j <= gridsize + 1; j++) {
					cout << " ";
				}
				cout << "\n";
			}*/
			/*SORT BELOW*/

			for (i = 1;i <= jParallel;i++) {
				for (j = 0; j <= gridsize+1; j++) {
					if (j == 0){
						bufU[i-1][j] = (4 * concBufUOLD[i][j] - bufUOLDER[i-1][j] + 2 * dt * (-gx * (concBufHOLD[i + 1][j] - concBufHOLD[i - 1][j]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufUOLD[i][j + 1] - concBufUOLD[i][gridsize + 1]) + f * concBufVOLD[i][j])) / 3;
						bufV[i-1][j] = (4 * concBufVOLD[i][j] - bufVOLDER[i-1][j] + 2 * dt * (-gy * (concBufHOLD[i][j + 1] - concBufHOLD[i][gridsize + 1]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufVOLD[i + 1][j] - concBufVOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufVOLD[i][j + 1] - concBufVOLD[i][gridsize + 1]) - f * concBufUOLD[i][j])) / 3;
						bufH[i-1][j] = (4 * concBufHOLD[i][j] - bufHOLDER[i-1][j] -2 * H * dt * ((concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) / (2 * dx) + (concBufVOLD[i][j + 1] - concBufVOLD[i][gridsize + 1]) / (2 * dy))) / 3;
					}
				else {
                                                bufU[i-1][j] = (4 * concBufUOLD[i][j] - bufUOLDER[i-1][j] + 2 * dt * (-gx * (concBufHOLD[i + 1][j] - concBufHOLD[i - 1][j]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufUOLD[i][(j + 1) % (gridsize + 2)] - concBufUOLD[i][j - 1]) + f * concBufVOLD[i][j])) / 3;
                                                bufV[i-1][j] = (4 * concBufVOLD[i][j] - bufVOLDER[i-1][j] + 2 * dt * (-gy * (concBufHOLD[i][(j + 1)% (gridsize + 2)] - concBufHOLD[i][j - 1]) - ((concBufUOLD[i][j]) / (2 * dx)) * (concBufVOLD[i + 1][j] - concBufVOLD[i - 1][j]) - ((concBufVOLD[i][j]) / (2 * dy)) * (concBufVOLD[i][(j + 1)% (gridsize + 2)] - concBufVOLD[i][j - 1]) - f * concBufUOLD[i][j])) / 3;
                                                bufH[i-1][j] = (4 * concBufHOLD[i][j] - bufHOLDER[i-1][j] -2 * H * dt * ((concBufUOLD[i + 1][j] - concBufUOLD[i - 1][j]) / (2 * dx) + (concBufVOLD[i][(j + 1)% (gridsize + 2)] - concBufVOLD[i][j - 1]) / (2 * dy))) / 3;
				}
			}
			}
                for (j = 0;j <= gridsize + 1;j++) {
                        for (i = 0; i <= jParallel - 1; i++) {
                                bufUOLDER[i][j] = bufUOLD[i][j];           /*MUST THINK ABOUT INDEXES*/
                                bufVOLDER[i][j] = bufVOLD[i][j];
                                bufHOLDER[i][j] = bufHOLD[i][j];
                        }
                }

                for (j = 0;j <= gridsize+1; j++) {
                        for (i = 0; i <= jParallel - 1; i++) {
                                bufUOLD[i][j] = bufU[i][j];
                                bufVOLD[i][j] = bufV[i][j];
                                bufHOLD[i][j] = bufH[i][j];
                        }
                }

		}
		MPI_Barrier(MPI_COMM_WORLD);

	}

        MPI_Gather(&bufU[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, &masterbufU[0][0], jParallel * (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&bufV[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, &masterbufV[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&bufH[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, &masterbufH[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&bufUOLD[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, &masterbufUOLD[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);   
        MPI_Gather(&bufVOLD[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, &masterbufVOLD[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&bufHOLD[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, &masterbufHOLD[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        MPI_Gather(&bufUOLDER[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, &masterbufUOLDER[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        MPI_Gather(&bufVOLDER[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, &masterbufVOLDER[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        MPI_Gather(&bufHOLDER[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, &masterbufHOLDER[0][0], jParallel* (gridsize + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        cout << "GATHERED\n";     

	if (rank == 0) {
		for (i = 0; i <= gridsize + 1; i++) {
			for (j = 0; j <= gridsize + 1; j++) {
				cout << masterbufH[i][j];
				cout << " ";
			}
			cout << "\n";
		}
	}
	t2 = MPI_Wtime();
	cout << "The size is "<< size << "\n"; 
	cout << "Total time elapsed is " << t2-t1 << "\n\n";
	MPI_Finalize();
	if (rank == 0) {
        cout << "\n";
	cout << "DAta erad time \n";
      // Create the file. The Replace parameter tells netCDF to overwrite
      // this file, if it already exists.
       NcFile dataFile("pog25sec.nc", NcFile::replace);
	cout << "\n Read the data successfully";
      // Create netCDF dimensions
       NcDim xDim = dataFile.addDim("x", gridsize+2);
       NcDim yDim = dataFile.addDim("y", gridsize+2);

      // Define the variable. The type of the variable in this case is
      // ncInt (32-bit integer).
       vector<NcDim> dims;
       dims.push_back(xDim);
       dims.push_back(yDim);
       NcVar data = dataFile.addVar("data", ncDouble, dims);

      // Write the data to the file. Although netCDF supports
      // reading and writing subsets of data, in this case we write all
      // the data in one operation.
       data.putVar(masterbufH);
      cout << "\n data read ok! \n";

      // The file will be automatically close when the NcFile object goes
      // out of scope. This frees up any internal netCDF resources
      // associated with the file, and flushes any buffers.

      //cout << "*** SUCCESS writing example file simple_xy.nc!" << endl;
      }
	return 0;
}
