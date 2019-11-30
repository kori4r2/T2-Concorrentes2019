#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <omp.h>

void calcula_menor(int *notas, int *menorc, int *menorr, int *menorb, int R, int C, int A)
{
    //valores default
    for (int r = 0; r < R; ++r)
    {
        menorr[r] = INT_MAX;
        for (int c = 0; c < C; ++c)
        {
            menorc[r * C + c] = INT_MAX;
        }
    }
    *menorb = notas[0];

    for (int r = 0; r < R; ++r)
    {
        for (int c = 0; c < C; ++c)
        {
            for (int a = 0; a < A; ++a)
            {
                if (notas[r * C * A + c * A + a] < menorc[r * C + c])
                {
                    menorc[r * C + c] = notas[r * C * A + c * A + a];
                    if (notas[r * C * A + c * A + a] < menorr[r])
                    {
                        menorr[r] = notas[r * C * A + c * A + a];
                        if (notas[r * C * A + c * A + a] < *menorb)
                            *menorb = notas[r * C * A + c * A + a];
                    }
                }
            }
        }
    }
}

void calcula_maior(int *notas, int *maiorc, int *maiorr, int *maiorb, int R, int C, int A)
{
    memset(maiorr, 0, R);
    memset(maiorc, 0, C);
    
    *maiorb = notas[0];

    for (int r = 0; r < R; ++r)
    {
        for (int c = 0; c < C; ++c)
        {
            for (int a = 0; a < A; ++a)
            {
                if (notas[r * C * A + c * A + a] > maiorc[r * C + c])
                {
                    maiorc[r * C + c] = notas[r * C * A + c * A + a];
                    if (notas[r * C * A + c * A + a] > maiorr[r])
                    {
                        maiorr[r] = notas[r * C * A + c * A + a];
                        if (notas[r * C * A + c * A + a] > *maiorb)
                            *maiorb = notas[r * C * A + c * A + a];
                    }
                }
            }
        }
    }
}


int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

void calcula_mediana(int *notas, double *medianac, double *medianar, double *medianab, int R, int C, int A)
{
    //mediana nas cidades
    for (int r = 0; r < R; ++r)
    {
        for (int c = 0; c < C; ++c)
        {
            //mediana nas cidades
            qsort(&notas[r * C * A + c * A], A, sizeof(int), cmpfunc);
            if (A % 2)
            {
                medianac[r * C + c] = notas[r * C * A + c * A + A / 2 ];
            }
            else
            {
                medianac[r * C + c] = (notas[r * C * A + c * A + A / 2] + notas[r * C * A + c * A + A / 2 - 1]) / 2.0;
            }
        }
        //mediana nas regioes
        qsort(&notas[r * C * A], C * A, sizeof(int), cmpfunc);
        if (C * A % 2)
        {
            medianar[r] = notas[r * C * A + C * A / 2];
        }
        else
        {
            medianar[r] = (notas[r * C * A + C * A / 2] + notas[r * C * A + C * A / 2 - 1]) / 2.0;
        }
    }
    //mediana no brasil
    qsort(notas, R * C * A, sizeof(int), cmpfunc);
    if (R * C * A % 2)
    {
        *medianab = notas[R * C * A / 2];
    }
    else
    {
        *medianab = (notas[R * C * A / 2] + notas[R * C * A / 2 - 1]) / 2.0;
    }
}

void calcula_media(int *notas, double *mediac, double *mediar, double *mediab, int *mr, int *mcr, int *mcc, int R, int C, int A)
{
    //valores default
    for (int r = 0; r < R; ++r)
    {
        mediar[r] = 0;
        for (int c = 0; c < C; ++c)
        {
            mediac[r * C + c] = 0;
        }
    }
    *mediab = 0;
    *mr = 0;
    *mcr = 0;
    *mcc = 0;

    for (int r = 0; r < R; ++r)
    {
        for (int c = 0; c < C; ++c)
        {
            for (int a = 0; a < A; ++a)
            {
                mediac[r * C + c] += notas[r * C * A + c * A + a];
            }
            mediac[r * C + c] /= A;
            mediar[r] += mediac[r * C + c];

            if (mediac[r * C + c] > mediac[(*mcr) * C + (*mcc)])
            {
                *mcr = r;
                *mcc = c;
            }
        }
        mediar[r] /= C;
        *mediab += mediar[r];

        if (mediar[r] > mediar[(*mr)])
        {
            *mr = r;
        }
    }
    *mediab /= R;
}

void calcula_dp(int *notas, double *dpc, double *dpr, double *dpb, double *mediac, double *mediar, double mediab, int R, int C, int A)
{
    //valores default
    double x = 0, y = 0, z = 0;

    for (int r = 0; r < R; ++r)
    {
        for (int c = 0; c < C; ++c)
        {
            for (int a = 0; a < A; ++a)
            {
                x += pow(notas[r * C * A + c * A + a] - mediac[r * C + c], 2);
                y += pow(notas[r * C * A + c * A + a] - mediar[r], 2);
                z += pow(notas[r * C * A + c * A + a] - mediab, 2);
            }
            dpc[r * C + c] = sqrt(x / A);
            x = 0;
        }
        dpr[r] = sqrt(y / (C * A));
        y = 0;
    }
    *dpb = sqrt(z / (R * C * A));
}

int main(int argc, char *argv[])
{
    int R, C, A, seed, size;
    double wtime;
    int *notas;
    int *menorc, *maiorc, *menorr, *maiorr, menorb, maiorb;
    double *medianac, *mediac, *dpc, *medianar, *mediar, *dpr, medianab, mediab, dpb;
    int mr;       //num da melhor regiao
    int mcr, mcc; //num da regiao da melhor cidade e da melhor cidade respectivamente

    //entrada parametros
	if(argc != 5){
		printf("Erro de entrada, o comando de execução deve ser ./Trabalho2 \'nº de regiões\' \'nº de Cidades por região \'nº de alunos por cidade\' \'seed para aleatorizar\'");
		return 1;
	}
	R = atoi(argv[1]);
	C = atoi(argv[2]);
	A = atoi(argv[3]);
	seed = atoi(argv[4]);

    srand(seed);

    size = R * C * A;
    notas = (int *)malloc(size * sizeof(int));

    menorc = (int *)malloc(R * C * sizeof(int));
    maiorc = (int *)malloc(R * C * sizeof(int));
    medianac = (double *)malloc(R * C * sizeof(double));
    mediac = (double *)malloc(R * C * sizeof(double));
    dpc = (double *)malloc(R * C * sizeof(double));

    menorr = (int *)malloc(R * sizeof(int));
    maiorr = (int *)malloc(R * sizeof(int));
    medianar = (double *)malloc(R * sizeof(double));
    mediar = (double *)malloc(R * sizeof(double));
    dpr = (double *)malloc(R * sizeof(double));

    //geracao das notas
    for (int i = 0; i < size; ++i)
        notas[i] = rand() % 101;
        
    // for (int i = 0; i < R; ++i) {
    //     for (int j = 0; j < C; ++j) {
    //         for (int k = 0; k < A; ++k) {
    //             printf("%d ", notas[i*C*A + j*A + k]);
    //         }
    //         printf("\n");
    //     }
    // }

    wtime = omp_get_wtime();

    calcula_menor(notas, menorc, menorr, &menorb, R, C, A);
    calcula_maior(notas, maiorc, maiorr, &maiorb, R, C, A);
    calcula_media(notas, mediac, mediar, &mediab, &mr, &mcr, &mcc, R, C, A);
    calcula_dp(notas, dpc, dpr, &dpb, mediac, mediar, mediab, R, C, A);
    calcula_mediana(notas, medianac, medianar, &medianab, R, C, A); //caga todas as notas tem que ser por ultimo

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

    free(notas);
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
    

    return 0;
}
