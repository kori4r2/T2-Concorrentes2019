#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <limits.h>
#include <omp.h>

#define NUM_THREADS 4

void calcula_menor(int *notas, int *menorc, int *menorr, int *menorb, int startReg, int nRegioes, int C, int A)
{
    int mb = INT_MAX;
	int endReg = startReg + nRegioes;
    #pragma omp parallel for reduction(min: mb) schedule(guided)
    for (int r = startReg; r < endReg; ++r)
    {
        menorr[r - startReg] = INT_MAX;
        int menorr_r = INT_MAX;
        #pragma omp parallel for reduction(min: menorr_r) schedule(guided)
        for (int c = 0; c < C; ++c)
        {
            menorc[(r - startReg) * C + c] = INT_MAX;
            int menorc_rCx = INT_MAX;
            #pragma omp parallel for reduction(min: menorc_rCx) schedule(guided)
            for (int a = 0; a < A; ++a)
            {
                if(notas[r * C * A + c * A + a] < menorc_rCx)
                    menorc_rCx = notas[r * C * A + c * A + a];
            }
            menorc[(r - startReg) * C + c] = menorc_rCx;
            if(menorc[(r - startReg) * C + c] < menorr_r)
                menorr_r = menorc[(r - startReg) * C + c];
        }
        menorr[r - startReg] = menorr_r;
        if(menorr[r - startReg] < mb)
            mb = menorr[r - startReg];
    }
    *menorb = mb;
}

void calcula_maior(int *notas, int *maiorc, int *maiorr, int *maiorb, int startReg, int nRegioes, int C, int A)
{
    //valores default
    memset(maiorr, 0, nRegioes * sizeof(int));
    memset(maiorc, 0, nRegioes * C * sizeof(int));
    
    int mb = 0;
	int endReg = startReg + nRegioes;
    #pragma omp parallel for reduction(max: mb) schedule(guided)
    for (int r = startReg; r < endReg; ++r)
    {
        int maiorr_r = maiorr[r - startReg];
        #pragma omp parallel for reduction(max: maiorr_r) schedule(guided)
        for (int c = 0; c < C; ++c)
        {
            int maiorc_rCc = maiorc[(r - startReg) * C + c];
            #pragma omp parallel for reduction(max: maiorc_rCc) schedule(guided)
            for (int a = 0; a < A; ++a)
            {
                if(notas[r * C * A + c * A + a] > maiorc_rCc)
                    maiorc_rCc = notas[r * C * A + c * A + a];
            }
            maiorc[(r - startReg) * C + c] = maiorc_rCc;
            if(maiorc[(r - startReg) * C + c] > maiorr_r)
                maiorr_r = maiorc[(r - startReg) * C + c];
        }
        maiorr[r - startReg] = maiorr_r;
        if(maiorr[r - startReg] > mb)
            mb = maiorr[r - startReg];
    }
    
    *maiorb = mb;
}

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

int mergeSortedBlocks(int *block, int**sortedBlocks, int firstSortedBlock, int blockSize, int size, int blockStart){
    if(size == blockSize){
        memcpy(block, sortedBlocks[firstSortedBlock + blockStart/blockSize], blockSize * sizeof(int));
        return 1;
    }
    
    int nBlocks = size / blockSize;
    int sizeBlockA = (nBlocks/2) * blockSize;
    int sizeBlockB = (nBlocks/2 + ((nBlocks % 2 == 0)? 0 : 1)) * blockSize;
    
    // Ordena a primeira metade do bloco
    int* blockA = (int*)malloc(sizeBlockA * sizeof(int));
    if(sizeBlockA == blockSize){
        memcpy(blockA, sortedBlocks[firstSortedBlock + blockStart/blockSize], blockSize * sizeof(int));
    }else{
        mergeSortedBlocks(blockA, sortedBlocks, firstSortedBlock, blockSize, sizeBlockA, blockStart);
    }
    
    // Ordena a segunda metade do bloco
    int* blockB = (int*)malloc(sizeBlockB * sizeof(int));
    if(sizeBlockB == blockSize){
        memcpy(blockB, sortedBlocks[firstSortedBlock + (blockStart + sizeBlockA)/blockSize], blockSize * sizeof(int));
    }else{
        mergeSortedBlocks(blockB, sortedBlocks, firstSortedBlock, blockSize, sizeBlockB, blockStart + sizeBlockA);
    }
    
    // Faz um merge dos dois blocos no bloco no atual
    int i, j, k, min;
    for(i = 0, j = 0, k = 0; (i < sizeBlockA || j < sizeBlockB) && k < size;k++){
        // Calcula qual o valor a ser colocado na proxima posição do bloco atual
        if(i >= sizeBlockA){
            min = blockB[j++];
        }else if(j >= sizeBlockB){
            min = blockA[i++];
        }else{
            if(blockB[j] <= blockA[i]){
                min = blockB[j++];
            }else{
                min = blockA[i++];
            }
        }
        
        block[k] = min;
    }
    
    return 0;
}

void calcula_mediana(int *notas, double *medianac, double *medianar, int startReg, int nRegioes, int C, int A)
{
    // Copia as notas de cada cidade separadamente e ordena usando quicksort
    int **notasOrdenadasCidades = (int**)malloc(nRegioes * C * sizeof(int*));
    #pragma omp parallel for schedule(guided)
    for(int i = 0; i < nRegioes * C; i++){
        notasOrdenadasCidades[i] = (int*)malloc(A * sizeof(int));
        memcpy(notasOrdenadasCidades[i], &notas[startReg * C * A + i * A], A * sizeof(int));
        qsort(notasOrdenadasCidades[i], A, sizeof(int), cmpfunc);
        if (A % 2)
        {
            medianac[i] = notasOrdenadasCidades[i][A / 2];
        }
        else
        {
            medianac[i] = (notasOrdenadasCidades[i][A / 2] + notasOrdenadasCidades[i][A / 2 - 1]) / 2.0;
        }
    }
    
    // Pega as notas de cidades ordenadas, e faz merge em blocos para cada região
    int **notasOrdenadasRegioes = (int**)malloc(nRegioes * sizeof(int*));
    #pragma omp parallel for schedule(guided)
    for(int i = 0; i < nRegioes; i++){
        notasOrdenadasRegioes[i] = (int*)malloc(C * A * sizeof(int));
        mergeSortedBlocks(notasOrdenadasRegioes[i], notasOrdenadasCidades, i * C, A, C * A, 0);
        if ((C * A) % 2)
        {
            medianar[i] = notasOrdenadasRegioes[i][(C * A) / 2];
        }
        else
        {
            medianar[i] = (notasOrdenadasRegioes[i][(C * A) / 2] + notasOrdenadasRegioes[i][(C * A) / 2 - 1]) / 2.0;
        }
    }
    
    // Não precisa mais da matriz de notas ordenadas das cidades
    #pragma omp parallel for schedule(guided)
    for(int i = 0; i < nRegioes * C; i++){
        free(notasOrdenadasCidades[i]);
    }
    free(notasOrdenadasCidades);

    // Não precisa mais da matriz de notas ordenadas das regioes
    #pragma omp parallel for schedule(guided)
    for(int i = 0; i < nRegioes; i++){
        free(notasOrdenadasRegioes[i]);
    }
    free(notasOrdenadasRegioes);
    
}

void calcula_media(int *notas, double *mediac, double *mediar, double *mediab, int *mr, int *mc, int startReg, int nRegioes, int C, int A)
{
    //valores default
    #pragma omp parallel for schedule(guided)
    for (int r = 0; r < nRegioes; ++r)
    {
        mediar[r] = 0;
        #pragma omp parallel for schedule(guided)
        for (int c = 0; c < C; ++c)
        {
            mediac[r * C + c] = 0;
        }
    }
    *mediab = 0;
	double mb = 0;

	int endReg = startReg + nRegioes;
    #pragma omp parallel for schedule(guided) reduction(+: mb)
    for (int r = startReg; r < endReg; ++r)
    {
        double mediar_r = mediar[r - startReg];
        #pragma omp parallel for schedule(guided) reduction(+: mediar_r)
        for (int c = 0; c < C; ++c)
        {
            double mediac_rCc = mediac[(r - startReg) * C + c];
            #pragma omp parallel for schedule(guided) reduction(+: mediac_rCc)
            for (int a = 0; a < A; ++a)
            {
                mediac_rCc += notas[r * C * A + c * A + a];
            }
            mediac[(r - startReg) * C + c] = mediac_rCc/A;
            mediar_r += mediac[(r - startReg) * C + c];

			#pragma omp critical
			{
				if (mediac[(r - startReg) * C + c] > mc[2] || ((mediac[(r - startReg) * C + c] == mc[2]) && (c < mc[1])))
				{
					mc[0] = r;
					mc[1] = c;
					mc[2] = mediac[(r - startReg) * C + c];
				}
			}
        }
        mediar[r - startReg] = mediar_r/C;
        mb += mediar[r - startReg];

		#pragma omp critical
		{
			if (mediar[r - startReg] > mr[1] || (mediar[r - startReg] == mr[1] && r < mr[0]))
			{
				mr[0] = r;
				mr[1] = mediar[r - startReg];
			}
		}
    }
	mb /= nRegioes;
    *mediab = mb;
}

void media_dp_aux(int *notas, double *mediac, double *mediar, double *mediab, int startReg, int nRegioes, int C, int A)
{
    //valores default
    #pragma omp parallel for schedule(guided)
    for (int r = 0; r < nRegioes; ++r)
    {
        mediar[r] = 0;
        #pragma omp parallel for schedule(guided)
        for (int c = 0; c < C; ++c)
        {
            mediac[r * C + c] = 0;
        }
    }
    *mediab = 0;
	double mb = 0;

	int endReg = startReg + nRegioes;
    #pragma omp parallel for schedule(guided) reduction(+: mb)
    for (int r = startReg; r < endReg; ++r)
    {
        double mediar_r = mediar[r - startReg];
        #pragma omp parallel for schedule(guided) reduction(+: mediar_r)
        for (int c = 0; c < C; ++c)
        {
            double mediac_rCc = mediac[(r - startReg) * C + c];
            #pragma omp parallel for schedule(guided) reduction(+: mediac_rCc)
            for (int a = 0; a < A; ++a)
            {
                mediac_rCc += notas[r * C * A + c * A + a];
            }
            mediac[(r - startReg) * C + c] = mediac_rCc/A;
            mediar_r += mediac[(r - startReg) * C + c];
        }
        mediar[r - startReg] = mediar_r/C;
        mb += mediar[r -startReg];
    }
	mb /= nRegioes;
    *mediab = mb;
}

void calcula_dp(int *notas, double *dpc, double *dpr, double *dpb, int startReg, int nRegioes, int C, int A)
{
    double *mediac, *mediar, mediab;
    mediac = (double *) malloc(nRegioes * C * sizeof(double));
    mediar = (double *) malloc(nRegioes * sizeof(double));
    
    //calcula a media
    media_dp_aux(notas, mediac, mediar, &mediab, startReg, nRegioes, C, A);
    
    //valores default
    double x = 0, y = 0, z = 0;
	int endReg = startReg + nRegioes;
    #pragma omp parallel for schedule(guided) reduction(+: z) private(y)
    for (int r = startReg; r < endReg; ++r)
    {
		y = 0;
        #pragma omp parallel for schedule(guided) reduction(+: y, z) private(x)
        for (int c = 0; c < C; ++c)
        {
			x = 0;
            #pragma omp parallel for reduction(+: x, y, z) schedule(guided)
            for (int a = 0; a < A; ++a)
            {
                x += pow(notas[r * C * A + c * A + a] - mediac[(r - startReg) * C + c], 2);
                y += pow(notas[r * C * A + c * A + a] - mediar[r - startReg], 2);
                z += pow(notas[r * C * A + c * A + a] - mediab, 2);
            }
            dpc[(r - startReg) * C + c] = sqrt(x / A);
        }
        dpr[r - startReg] = sqrt(y / (C * A));
    }
    *dpb = sqrt(z / (nRegioes * C * A));
    
    free(mediac);
    free(mediar);
}

int main (int argc, char *argv[]){
    int R, C, A, seed, NTA, w_rank, w_size, recvBuffer;
    double wtime;
    int *notas;
    int *menorc, *maiorc, *menorr, *maiorr, menorb, maiorb;
    double *medianac, *mediac, *dpc, *medianar, *mediar, *dpr, medianab, mediab, dpb;
    int mr[2] = {-1, -1}; // pos 0 = numero da melhor regiao; pos 1 = valor da media na melhor regiao
    int mc[3] = {-1, -1, -1}; // pos 0 = numero da regiao da melhor cidade; pos 1 = numero da melhor cidade; pos 2 = valor da media na melhor cidade
	MPI_Comm commWorld = MPI_COMM_WORLD;

//--------------- Código executado no inicio de cada processo, avalia os argumentos de entrada e inicializa a matriz (Replicação de 
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
    for (int i = 0; i < NTA; ++i){
        notas[i] = rand() % 101;
	}

	omp_set_num_threads(NUM_THREADS);
	omp_set_nested(1);

//--------------- A partir daqui começam as comunicações usando o MPI
    wtime = omp_get_wtime();
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &w_size);

	// Se o processo for desnecessário, ele é encerrado
	if(w_size > R){
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

	int primeiraRegiao = (R / w_size) * w_rank;
	if((R % w_size) - (w_size - w_rank) > 0)
		primeiraRegiao += ((R % w_size) - (w_size - w_rank));

//--------------- Alocação de memoria necessaria
	int *distribuicaoCidades = NULL;
	int *offsetCidades = NULL;
	int *distribuicaoRegioes = NULL;
	int *offsetRegioes = NULL;
	int *melhorRegiao = NULL;
	int *melhorCidade = NULL;
	double *mediabProcessos = NULL;
	double *dpbProcessos = NULL;
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
		distribuicaoCidades = (int*)malloc(sizeof(int) * w_size);
		offsetCidades = (int*)malloc(sizeof(int) * w_size);
		distribuicaoRegioes = (int*)malloc(sizeof(int) * w_size);
		offsetRegioes = (int*)malloc(sizeof(int) * w_size);

		for(int i = 0; i < w_size; i++){
			// Calcula a quantos elementos o processo i tem
			distribuicaoRegioes[i] = R / w_size;
			if((w_size - i) <= (R % w_size))
				distribuicaoRegioes[i]++;
			distribuicaoCidades[i] = distribuicaoRegioes[i] * C;

			// Calcula o offset do processo i nos vetores de resultados para cidades e regioes
			if(i > 0){
				offsetCidades[i] = offsetCidades[i-1] + distribuicaoCidades[i-1];
				offsetRegioes[i] = offsetRegioes[i-1] + distribuicaoRegioes[i-1];
			}else{
				offsetCidades[i] = 0;
				offsetRegioes[i] = 0;
			}
		}

		// Ele tambem precisa de arrays auxiliares para armazenar as medias e dp brasileiras calculadas em cada processo, assim como index de melhor regiao e cidade
		mediabProcessos = (double*)malloc(sizeof(double) * w_size);
		dpbProcessos = (double*)malloc(sizeof(double) * w_size);
		melhorCidade = (int*)malloc(sizeof(int) * 3 * w_size);
		melhorRegiao = (int*)malloc(sizeof(int) * 2 * w_size);
	}else{
		menorc = NULL;
		maiorc = NULL;
		medianac = NULL;
		mediac = NULL;
		dpc = NULL;
		menorr = NULL;
		maiorr = NULL;
		medianar = NULL;
		mediar = NULL;
		dpr = NULL;
	}

	// Cada Recebe o numero de cidades e regioes em cada processo para definir os pesos dos receives
	// Pode ser mais eficiente replicar o calculo
//	MPI_Gather(&nCidades, 1, MPI_INT, distribuicaoCidades, 1, MPI_INT, 0, commWorld);
//	MPI_Gather(&nRegioes, 1, MPI_INT, distribuicaoRegioes, 1, MPI_INT, 0, commWorld);

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
			calcula_menor(notas, menorc_aqui, menorr_aqui, &menorb, primeiraRegiao, nRegioes, C, A);
		}
		#pragma omp section
		{
			// calcula_maior
			calcula_maior(notas, maiorc_aqui, maiorr_aqui, &maiorb, primeiraRegiao, nRegioes, C, A);
		}
		#pragma omp section
		{
			// calcula_mediana 
            calcula_mediana(notas, medianac_aqui, medianar_aqui, primeiraRegiao, nRegioes, C, A);
		}
		#pragma omp section
		{
			// calcula_media
            calcula_media(notas, mediac_aqui, mediar_aqui, &mediab, mr, mc, primeiraRegiao, nRegioes, C, A);
		}
		#pragma omp section
		{
			// calcula_dp
			calcula_dp(notas, dpc_aqui, dpr_aqui, &dpb, primeiraRegiao, nRegioes, C, A);
		}
		#pragma omp section
		{
//--------------- Ordenação de todas as notas para calculo de mediana do brasil (processo 0)
			if(w_rank == 0){
				int *notas2 = (int*)malloc(sizeof(int) * NTA);
				memcpy(notas2, notas, sizeof(int) * NTA);
				qsort(notas2, NTA, sizeof(int), cmpfunc);
				if(NTA % 2){
					medianab = notas2[NTA / 2];
				}else{
					medianab = (notas2[NTA / 2] + notas2[NTA / 2 - 1]) / 2.0;
				}
				free(notas2);
			}
		}
	}
//--------------- Gatherv no processo 0 para calculo dos valores do brasil

	// Recebe os resultados dos calculos para cada cidade
	MPI_Gatherv(menorc_aqui, nCidades, MPI_INT, menorc, distribuicaoCidades, offsetCidades, MPI_INT, 0, commWorld);
	MPI_Gatherv(maiorc_aqui, nCidades, MPI_INT, maiorc, distribuicaoCidades, offsetCidades, MPI_INT, 0, commWorld);
	MPI_Gatherv(medianac_aqui, nCidades, MPI_DOUBLE, medianac, distribuicaoCidades, offsetCidades, MPI_DOUBLE, 0, commWorld);
	MPI_Gatherv(mediac_aqui, nCidades, MPI_DOUBLE, mediac, distribuicaoCidades, offsetCidades, MPI_DOUBLE, 0, commWorld);
	MPI_Gatherv(dpc_aqui, nCidades, MPI_DOUBLE, dpc, distribuicaoCidades, offsetCidades, MPI_DOUBLE, 0, commWorld);

	// Recebe os resultados dos calculos para cada regiao
	MPI_Gatherv(menorr_aqui, nRegioes, MPI_INT, menorr, distribuicaoRegioes, offsetRegioes, MPI_INT, 0, commWorld);
	MPI_Gatherv(maiorr_aqui, nRegioes, MPI_INT, maiorr, distribuicaoRegioes, offsetRegioes, MPI_INT, 0, commWorld);
	MPI_Gatherv(medianar_aqui, nRegioes, MPI_DOUBLE, medianar, distribuicaoRegioes, offsetRegioes, MPI_DOUBLE, 0, commWorld);
	MPI_Gatherv(mediar_aqui, nRegioes, MPI_DOUBLE, mediar, distribuicaoRegioes, offsetRegioes, MPI_DOUBLE, 0, commWorld);
	MPI_Gatherv(dpr_aqui, nRegioes, MPI_DOUBLE, dpr, distribuicaoRegioes, offsetRegioes, MPI_DOUBLE, 0, commWorld);

	// Usa um reduce para calcular o menor e o maior valor do brasil
	MPI_Reduce(&menorb, &recvBuffer, 1, MPI_INT, MPI_MIN, 0, commWorld);
	menorb = recvBuffer;
	MPI_Reduce(&maiorb, &recvBuffer, 1, MPI_INT, MPI_MAX, 0, commWorld);
	maiorb = recvBuffer;

	// Recebe as medias e desvios padrao calculados por cada processo
	MPI_Gather(&mediab, 1, MPI_DOUBLE, mediabProcessos, 1, MPI_DOUBLE, 0, commWorld);
	MPI_Gather(&dpb, 1, MPI_DOUBLE, dpbProcessos, 1, MPI_DOUBLE, 0, commWorld);

	// Recebe os resultados de melhor regiao e cidade de cada processo
	MPI_Gather(mr, 2, MPI_INT, melhorRegiao, 2, MPI_INT, 0, commWorld);
	MPI_Gather(mc, 3, MPI_INT, melhorCidade, 3, MPI_INT, 0, commWorld);
	
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

		// Calcula a media e o desvio padrao global entre os processos
		mediab = 0;
		for(int i = 0; i < w_size; i++){
			mediab += distribuicaoRegioes[i] * mediabProcessos[i];
		}
		mediab /= R;

		dpb = 0;
		for(int i = 0; i < w_size; i++){
			double desvioGlobal = mediabProcessos[i] - mediab;
			dpb += distribuicaoRegioes[i] * ((dpbProcessos[i] * dpbProcessos[i]) + (desvioGlobal * desvioGlobal));
		}
		dpb /= R;
		dpb = sqrt(dpb);

		// Avalia qual a melhor regiao e qual a melhor cidade
		mr[0] = -1;
		mr[1] = -1;
		mc[0] = -1;
		mc[1] = -1;
		mc[2] = -1;
		for(int i = 0; i < w_size; i++){
			if(melhorRegiao[i*2 + 1] > mr[1]){
				mr[1] = melhorRegiao[i*2 + 1];
				mr[0] = melhorRegiao[i*2];
			}

			if(melhorCidade[i*3 + 2] > mc[2]){
				mc[2] = melhorCidade[i*3 + 2];
				mc[1] = melhorCidade[i*3 + 1];
				mc[0] = melhorCidade[i*3];
			}
		}

		// Calcula quanto tempo demorou
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
		printf("Melhor região: Região %d\n", mr[0]);
		printf("Melhor cidade: Região %d, Cidade %d\n", mc[0], mc[1]);
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
		free(mediabProcessos);
		free(dpbProcessos);
		free(melhorCidade);
		free(melhorRegiao);
		// Libera a memoria das arrays de contagem de elementos
		free(distribuicaoCidades);
		free(distribuicaoRegioes);
	}


	// Encerra as comunicações do MPI
	MPI_Finalize();

	// Libera a memoria alocada antes da comunicação por MPI começar
    free(notas);

	return 0;
}
