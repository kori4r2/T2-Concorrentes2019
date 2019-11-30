#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <omp.h>

#define SLAVE_PROGRAM "slave_cidade"
#define NUM_THREADS 4

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

int main (int argc, char *argv[]){
    int R, C, A, seed, NTA, w_rank, w_size;
    double wtime;
    int *notas;
    int *menorc, *maiorc, *menorr, *maiorr, menorb, maiorb;
    double *medianac, *mediac, *dpc, *medianar, *mediar, *dpr, medianab, mediab, dpb;
    int mr;       //num da melhor regiao
    int mcr, mcc; //num da regiao da melhor cidade e da melhor cidade respectivamente
	MPI_Comm commWorld = MPI_COMM_WORLD;

//--------------- Código executado no inicio de cada processo, avalia os argumentos de entrada e inicializa a matriz
    //entrada parametros
	if(argc != 5){
		printf("Erro de entrada, o comando de execução deve ser mpirun -np \'nº de processos\' Trabalho2 \'nº de regiões\' \'nº de Cidades por região \'nº de alunos por cidade\' \'seed para aleatorizar\'");
		return 1;
	}
	R = atoi(argv[1]);
	C = atoi(argv[2]);
	A = atoi(argv[3]);
	seed = atoi(argv[4]);

    srand(seed);

    NTA = R * C * A;
    notas = (int *)malloc(NTA * sizeof(int));
    //geracao das notas
    for (int i = 0; i < NTA; ++i)
        notas[i] = rand() % 101;

	omp_set_num_threads(NUM_THREADS);
	omp_set_nested(1);

//--------------- A partir daqui começam as comunicações usando o MPI
    wtime = omp_get_wtime();
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &w_size);

	// Se o processo for desnecessário, ele é encerrado
	if(w_size > R){
		MPI_Group activeProcesses;
		int color = 0;
		if(w_rank >= R){
			color = 1;
		}else{
			color = 0;
		}
		MPI_Comm_split(MPI_COMM_WORLD, color, w_rank, &commWorld);
		MPI_Comm_size(commWorld, &w_size);
		MPI_Comm_rank(commWorld, &w_rank);
		if(w_rank >= R){
			MPI_Finalize();
			free(notas);
			return 0;
		}
	}

	// Calcula a carga de dados de cada processo
	// Distribui a carga igualmente, mas tenta deixar o processo 0 mais livre quando possível
	int nRegioes = R / w_size;
	if((w_size - w_rank) <= (R % w_size))
		nRegioes++;
	int nCidades = nRegioes * C;
	int nAlunos = nCidades * A;

	int startIndex = nRegioes * w_rank;
	if((R % w_size) - (w_size - w_rank) > 0)
		startIndex += ((R % w_size) - (w_size - w_rank));

//--------------- Alocação de memoria necessaria
	// Alocação de memoria, o processo 0 precisa ser capaz de armazenar todas as informações para realizar um gatherv
	if(w_rank == 0){
		// menorc armazena os menores valores de cada cidade
		menorc = (int *)malloc(R * C * sizeof(int));
		// maiorc armazena os maiores valores de cada cidade
		maiorc = (int *)malloc(R * C * sizeof(int));
		// medianac armazena as medianas de cada cidade
		medianac = (double *)malloc(R * C * sizeof(double));
		// mediac armazena as medias de cada cidade
		mediac = (double *)malloc(R * C * sizeof(double));
		// dpc armazena os desvios padrao de cada cidade
		dpc = (double *)malloc(R * C * sizeof(double));

		// menorr armazena os menores valores de cada regiao
		menorr = (int *)malloc(R * sizeof(int));
		// maiorr armazena os maiores valores de cada regiao
		maiorr = (int *)malloc(R * sizeof(int));
		// medianarr armazena as medianas de cada regiao
		medianar = (double *)malloc(R * sizeof(double));
		// mediar armazena as medias de cada regiao
		mediar = (double *)malloc(R * sizeof(double));
		// dpr armazena os desvios padrao de cada regiao
		dpr = (double *)malloc(R * sizeof(double));

		// O processo 0 precisa de arrays de contagem de elementos para realizar os gatherv
	}
	// Todos os processos precisam receber ter memoria suficiente para armazenar os dados das suas regiões -----------
	// menorc_aqui armazena os menores valores de cada cidade dentro das suas regioes
	int *menorc_aqui = (int *)malloc(nCidades * sizeof(int));
	// maiorc_aqui armazena os maiores valores de cada cidade dentro das suas regioes
	int *maiorc_aqui = (int *)malloc(nCidades * sizeof(int));
	// medianac_aqui armazena as medianas de cada cidade dentro das suas regioes
	double *medianac_aqui = (double *)malloc(nCidades * sizeof(double));
	// mediac_aqui armazena as medias de cada cidade dentro das suas regioes
	double *mediac_aqui = (double *)malloc(nCidades * sizeof(double));
	// dpc_aqui armazena os desvios padrao de cada cidade dentro das suas regioes
	double *dpc_aqui = (double *)malloc(nCidades * sizeof(double));

	// menorr_aqui armazena os menores valores de cada regiao sua
	int *menorr_aqui = (int *)malloc(nRegioes * sizeof(int));
	// maiorr_aqui armazena os maiores valores de cada regiao sua
	int *maiorr_aqui = (int *)malloc(nRegioes * sizeof(int));
	// medianarr_aqui armazena as medianas de cada regiao sua
	double *medianar_aqui = (double *)malloc(nRegioes * sizeof(double));
	// mediar_aqui armazena as medias de cada regiao sua
	double *mediar_aqui = (double *)malloc(nRegioes * sizeof(double));
	// dpr_aqui armazena os desvios padrao de cada regiao sua
	double *dpr_aqui = (double *)malloc(nRegioes * sizeof(double));

//--------------- Calcula os valores para todas as regioes processadas
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			// calcula_menor
		}
		#pragma omp section
		{
			// calcula_maior
		}
		#pragma omp section
		{
			// calcula_mediana (OBS.: alterado)
		}
		#pragma omp section
		{
			// calcula_media
		}
		#pragma omp section
		{
			// calcula_dp
		}
		#pragma omp section
		{
//--------------- Ordenação de todas as notas para calculo de mediana do brasil (processo 0)
			if(w_rank == 0){
				int *notas2 (int*)malloc(sizeof(int) * NTA);
				memcpy(notas2, notas, sizeof(int) * NTA);
				qsort(notas2, NTA, sizeof(int), cmpfunc);
				if(NTA % 2){
					medianab = notas[NTA / 2];
				}else{
					medianab = (notas[NTA / 2] + notas[NTA / 2 - 1]) / 2.0;
				}
			}
		}
	}
//--------------- Gather no processo 0 para calculo dos valores do brasil

//--------------- Liberação de memoria de todos os processos
	// Libera a memória alocada
	free(menorc_aqui);
	free(maiorc_aqui);
	free(medianac_aqui);
	free(mediac_aqui);
	free(dpc_aqui);
	free(menorr_aqui);
	free(maiorr_aqui);
	free(medianar_aqui);
	free(mediar_aqui);
	free(dpr_aqui);


//--------------- O processo zero faz a impressão dos resultados
	if(w_rank == 0){
		wtime = omp_get_wtime() - wtime;

		//printando o resultado das cidades
		for (int r = 0; r < R; r++)
		{
			for (int c = 0; c < C; c++)
			{
				printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n", r, c, menorc[r * C + c], maiorc[r * C + c], medianac[r * C + c], mediac[r * C + c], dpc[r * C + c]);
			}
			printf("\n");
		}

		//printando o resultado das regiões
		for (int r = 0; r < R; r++)
		{
			printf("Reg %d menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n", r, menorr[r], maiorr[r], medianar[r], mediar[r], dpr[r]);
		}
		printf("\n");

		//printando o resultado do Brasil
		printf("Brasil: menor: %d, maior: %d, mediana: %.2f, média: %.2f e DP: %.2f\n\n", menorb, maiorb, medianab, mediab, dpb);

		//printando melhores
		printf("Melhor região: Região %d\n", mr);
		printf("Melhor cidade: Região %d, Cidade %d\n", mcr, mcc);
		printf("\n");

		printf("Tempo de resposta sem considerar E/S, em segundos: %.3lfs\n", wtime);

		// Libera a memoria alocada exclusivamente no processo 0
		free(menorc);
		free(maiorc);
		free(medianac);
		free(mediac);
		free(dpc);
		free(menorr);
		free(maiorr);
		free(medianar);
		free(mediar);
		free(dpr);
		// Libera a memoria das arrays de contagem de elementos
	}


	// Encerra as comunicações do MPI
	MPI_Finalize();

	// Libera a memoria alocada antes da comunicação por MPI começar
    free(notas);

	return 0;
}
