//0-th qubit represent the LSB

#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

typedef unsigned int _uint;

void print_bin(_uint a, _uint n)
{
	int marker = 0;
	for(int i = n - 1; i >= 4; i--)
	{
		marker++;
		if(a & (1 << i)) printf("1");
		else printf("0");
		if(marker % 8 == 0) printf(" ");
	}
}

void add_key(Qureg *qubits, int key)
{
	int base = 4;
	while(key)
	{
		if(key & 1) pauliX((*qubits), base);
		key >>= 1;
		base++;
	}
}

void p_box(Qureg *qubits)
{
	swapGate((*qubits), 5, 8);
	swapGate((*qubits), 6, 10);
	swapGate((*qubits), 7, 9);
}

void s_box(Qureg *qubits)
{
	const int m = 4;
	int S_Box_anf[][12] = { { 0xE, 0xC, 0x8, 0x7, 0x6, 0x4, 0x1, 0x0, 0, 0, 0, 0 }, 
							{ 0xB, 0xA, 0x8, 0x5, 0x4, 0x3, 0x0, 0, 0, 0, 0, 0 },
							{ 0xE, 0xD, 0xC, 0xA, 0x9, 0x6, 0x5, 0x3, 0x2, 0x1, 0x0, 0 },
							{ 0xB, 0x9, 0x8, 0x5, 0x2, 0, 0, 0, 0, 0, 0, 0 } };
	int S_Box_size[4] = { 8, 7, 11, 5 };
	int S_Box_Inverse_anf[][9] = {  { 0x0, 0x1, 0x2, 0x4, 0x7, 0x8, 0, 0, 0 }, 
									{ 0x0, 0x1, 0x4, 0x9, 0xA, 0xB, 0xE, 0, 0 },
									{ 0x0, 0x2, 0x4, 0x9, 0xC, 0xD, 0, 0, 0 },
									{ 0x1, 0x3, 0x4, 0x6, 0x7, 0x8, 0x9, 0xC, 0xD } };
	int S_Box_Inverse_size[4] = { 6, 7, 6, 9 };

	ComplexMatrix2 ux = {
		.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
	};

	for(int byte = 1; byte <= 2; byte++)
	{
		int base = 3;
		for(int fi = 0; fi < m; fi++)
		{
			for(int i = 0; i < S_Box_size[fi]; i++)
			{
				if(S_Box_anf[fi][i] == 0x0)
				{
					pauliX((*qubits), base);
					continue;
				}
				int *ctrls = (int *)malloc(sizeof(int) * m);
				int ctrl_size = 0;
				int term = S_Box_anf[fi][i];
				int qb = byte * m;
				while(term)
				{
					if(term & 1) ctrls[ctrl_size++] = qb;
					term >>= 1;
					qb++;
				}
				multiControlledUnitary((*qubits), ctrls, ctrl_size, base, ux);
				free(ctrls);
			}
			base--;
		}

		base = byte * m + 3;
		for(int fi = 0; fi < m; fi++)
		{
			for(int i = 0; i < S_Box_Inverse_size[fi]; i++)
			{
				if(S_Box_Inverse_anf[fi][i] == 0x0)
				{
					pauliX((*qubits), base);
					continue;
				}
				int *ctrls = (int *)malloc(sizeof(int) * m);
				int ctrl_size = 0;
				int term = S_Box_Inverse_anf[fi][i];
				int qb = 0;
				while(term)
				{
					if(term & 1) ctrls[ctrl_size++] = qb;
					term >>= 1;
					qb++;
				}
				multiControlledUnitary((*qubits), ctrls, ctrl_size, base, ux);
				free(ctrls);
			}
			base--;
		}
		base = byte * m;
		for(int i = 0; i < m; i++)
			controlledNot((*qubits), i, i + base);
		for(int i = base; i < base + m; i++)
			controlledNot((*qubits), i, i - base);
	}
}

int main (int narg, char *varg[])
{
	int n = 8;
	int req_qubits = 2 * n + 4;

	QuESTEnv env = createQuESTEnv();

	Qureg qubits = createQureg(req_qubits, env);
	initZeroState(qubits);
	//reportQuregParams(qubits);
	//reportQuESTEnv(env);

	//Superposition for all plaintext
	for(int i = req_qubits - n; i < req_qubits; i++)
		hadamard(qubits, i);
    for(int i = req_qubits - n; i < req_qubits; i++)
        controlledNot(qubits, i, i - n);

	//START - Encryption
	//SubKey 0: E3 (1110 0011)
	add_key(&qubits, 0xE3);
	s_box(&qubits);
	p_box(&qubits);

	//SubKey 1: 2B (0010 1011)
	add_key(&qubits, 0x2B);
	s_box(&qubits);
	p_box(&qubits);

	//SubKey 2: F5 (1111 0101)
	add_key(&qubits, 0xF5);
	s_box(&qubits);
	p_box(&qubits);

	//SubKey 3: 9A (1001 1010)
	add_key(&qubits, 0x9A);
	s_box(&qubits);

	//SubKey 4: 75 (0111 0101)
	add_key(&qubits, 0x75);
	//ENDS - Encryption

	qreal prob = 0.0;
	for(int i = 0; i < (int)pow(2, req_qubits); i++)
	{
		prob = getProbAmp(qubits, i);
		if((int)(prob * 1000000) == 0) continue;
		printf("Prob. of ");
		print_bin(i, req_qubits);
		printf(" state is %f with amp %f\n", prob, getAmp(qubits, i).real);
	}

	destroyQureg(qubits, env);
	destroyQuESTEnv(env);

	return 0;
}