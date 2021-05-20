//0-th qubit represent the LSB

#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

void add_key(Qureg *qubits, int key, int base_idx)
{
	while(key)
	{
		if(key & 1) pauliX((*qubits), base_idx);
		key >>= 1;
		base_idx++;
	}
}

void p_box(Qureg *qubits, int base_idx)
{
	swapGate((*qubits), 1 + base_idx, 4 + base_idx);
	swapGate((*qubits), 2 + base_idx, 6 + base_idx);
	swapGate((*qubits), 3 + base_idx, 5 + base_idx);
}

void s_box(Qureg *qubits, int base_idx)
{
	int s_box_width = 4;

	//S-Box and its inverse in ANF representation
	int S_Box_anf[4][12] = { { 0xE, 0xC, 0x8, 0x7, 0x6, 0x4, 0x1, 0x0, 0, 0, 0, 0 }, 
							{ 0xB, 0xA, 0x8, 0x5, 0x4, 0x3, 0x0, 0, 0, 0, 0, 0 },
							{ 0xE, 0xD, 0xC, 0xA, 0x9, 0x6, 0x5, 0x3, 0x2, 0x1, 0x0, 0 },
							{ 0xB, 0x9, 0x8, 0x5, 0x2, 0, 0, 0, 0, 0, 0, 0 } };
	int S_Box_size[4] = { 8, 7, 11, 5 };
	int S_Box_Inverse_anf[4][12] = {  { 0x0, 0x1, 0x2, 0x4, 0x7, 0x8, 0, 0, 0, 0, 0, 0 }, 
									{ 0x0, 0x1, 0x4, 0x9, 0xA, 0xB, 0xE, 0, 0, 0, 0, 0 },
									{ 0x0, 0x2, 0x4, 0x9, 0xC, 0xD, 0, 0, 0, 0, 0, 0 },
									{ 0x1, 0x3, 0x4, 0x6, 0x7, 0x8, 0x9, 0xC, 0xD, 0, 0, 0 } };
	int S_Box_Inverse_size[4] = { 6, 7, 6, 9 };

	ComplexMatrix2 ux = {
		.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
	};

	for(int sbox_num = 0; sbox_num <= 1; sbox_num++)
	{
		int base = 3;
		for(int fi = 0; fi < s_box_width; fi++)
		{
			for(int i = 0; i < S_Box_size[fi]; i++)
			{
				if(S_Box_anf[fi][i] == 0x0)
				{
					pauliX((*qubits), base);
					continue;
				}
				int *ctrls = (int *)malloc(sizeof(int) * s_box_width);
				int ctrl_size = 0;
				int term = S_Box_anf[fi][i];
				int qb = base_idx + (sbox_num * s_box_width);
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

		base = base_idx + (sbox_num * s_box_width) + s_box_width - 1;
		for(int fi = 0; fi < s_box_width; fi++)
		{
			for(int i = 0; i < S_Box_Inverse_size[fi]; i++)
			{
				if(S_Box_Inverse_anf[fi][i] == 0x0)
				{
					pauliX((*qubits), base);
					continue;
				}
				int *ctrls = (int *)malloc(sizeof(int) * s_box_width);
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
		
		base = base_idx + (sbox_num * s_box_width);
		for(int i = 0; i < s_box_width; i++)
			controlledNot((*qubits), i, i + base);
		for(int i = base; i < base + s_box_width; i++)
			controlledNot((*qubits), i, i - base);
	}
}

void s_box_inverse(Qureg *qubits, int base_idx)
{
	int s_box_width = 4;

	//S-Box and its inverse in ANF representation
	int S_Box_anf[4][12] = { { 0x0, 0x1, 0x2, 0x4, 0x7, 0x8, 0, 0, 0, 0, 0, 0 }, 
							 { 0x0, 0x1, 0x4, 0x9, 0xA, 0xB, 0xE, 0, 0, 0, 0, 0 },
							 { 0x0, 0x2, 0x4, 0x9, 0xC, 0xD, 0, 0, 0, 0, 0, 0 },
							 { 0x1, 0x3, 0x4, 0x6, 0x7, 0x8, 0x9, 0xC, 0xD, 0, 0, 0 } };
	int S_Box_size[4] = { 6, 7, 6, 9 };
	int S_Box_Inverse_anf[4][12] = { 	{ 0xE, 0xC, 0x8, 0x7, 0x6, 0x4, 0x1, 0x0, 0, 0, 0, 0 }, 
										{ 0xB, 0xA, 0x8, 0x5, 0x4, 0x3, 0x0, 0, 0, 0, 0, 0 },
										{ 0xE, 0xD, 0xC, 0xA, 0x9, 0x6, 0x5, 0x3, 0x2, 0x1, 0x0, 0 },
										{ 0xB, 0x9, 0x8, 0x5, 0x2, 0, 0, 0, 0, 0, 0, 0 } };
	int S_Box_Inverse_size[4] = { 8, 7, 11, 5 };

	ComplexMatrix2 ux = {
		.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
	};

	for(int sbox_num = 0; sbox_num <= 1; sbox_num++)
	{
		int base = 3;
		for(int fi = 0; fi < s_box_width; fi++)
		{
			for(int i = 0; i < S_Box_size[fi]; i++)
			{
				if(S_Box_anf[fi][i] == 0x0)
				{
					pauliX((*qubits), base);
					continue;
				}
				int *ctrls = (int *)malloc(sizeof(int) * s_box_width);
				int ctrl_size = 0;
				int term = S_Box_anf[fi][i];
				int qb = base_idx + (sbox_num * s_box_width);
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

		base = base_idx + (sbox_num * s_box_width) + s_box_width - 1;
		for(int fi = 0; fi < s_box_width; fi++)
		{
			for(int i = 0; i < S_Box_Inverse_size[fi]; i++)
			{
				if(S_Box_Inverse_anf[fi][i] == 0x0)
				{
					pauliX((*qubits), base);
					continue;
				}
				int *ctrls = (int *)malloc(sizeof(int) * s_box_width);
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
		
		base = base_idx + (sbox_num * s_box_width);
		for(int i = 0; i < s_box_width; i++)
			controlledNot((*qubits), i, i + base);
		for(int i = base; i < base + s_box_width; i++)
			controlledNot((*qubits), i, i - base);
	}
}

void oracle_enc(Qureg *qubits, int n, int base_idx)
{
	//START - Encryption
	//Create entanglement of all plaintext qubits to ciphertext storage qubits
    for(int i = base_idx + 8; i < base_idx + 16; i++)
        controlledNot((*qubits), i, i - n);

	//SubKey 0: E3 (1110 0011)
	add_key(&(*qubits), 0xE3, base_idx);
	s_box(&(*qubits), base_idx);
	p_box(&(*qubits), base_idx);

	//SubKey 1: 2B (0010 1011)
	add_key(&(*qubits), 0x2B, base_idx);
	s_box(&(*qubits), base_idx);
	p_box(&(*qubits), base_idx);

	//SubKey 2: F5 (1111 0101)
	add_key(&(*qubits), 0xF5, base_idx);
	s_box(&(*qubits), base_idx);
	p_box(&(*qubits), base_idx);

	//SubKey 3: 9A (1001 1010)
	add_key(&(*qubits), 0x9A, base_idx);
	s_box(&(*qubits), base_idx);

	//SubKey 4: 75 (0111 0101)
	add_key(&(*qubits), 0x75, base_idx);
	//ENDS - Encryption
}

void oracle_dec(Qureg *qubits, int n, int base_idx)
{
	//START - Decryption
	//SubKey 4: 75 (0111 0101)
	add_key(&(*qubits), 0x75, base_idx);
	
	//SubKey 3: 9A (1001 1010)
	s_box_inverse(&(*qubits), base_idx);
	add_key(&(*qubits), 0x9A, base_idx);
		
	//SubKey 2: F5 (1111 0101)
	p_box(&(*qubits), base_idx);
	s_box_inverse(&(*qubits), base_idx);
	add_key(&(*qubits), 0xF5, base_idx);
	
	//SubKey 1: 2B (0010 1011)
	p_box(&(*qubits), base_idx);
	s_box_inverse(&(*qubits), base_idx);
	add_key(&(*qubits), 0x2B, base_idx);
		
	//SubKey 0: E3 (1110 0011)
	p_box(&(*qubits), base_idx);
	s_box_inverse(&(*qubits), base_idx);
	add_key(&(*qubits), 0xE3, base_idx);
	
	//Remove entanglement of all plaintext qubits to ciphertext storage qubits
    for(int i = base_idx + 8; i < base_idx + 16; i++)
        controlledNot((*qubits), i, i - n);
	//ENDS - Decryption
}

int main (int narg, char *varg[])
{	
	int subkey = 0x05; // Change to 0x70 to find suitable second linear characteristic

	int n = 8;

	//First n qubits for plaintext
	//Next n qubits for ciphertext
	//Last 4 qubits (LSB) ancilla qubits
	int req_qubits = 2 * n + 4;

	QuESTEnv env = createQuESTEnv();

	Qureg qubits = createQureg(req_qubits, env);
	initZeroState(qubits);
	//reportQuregParams(qubits);
	//reportQuESTEnv(env);

	ComplexMatrix2 ux = {
		.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
	};

	int base_idx = 4;
	//Superposition for all plaintext in qubits indexed from (2 * n - 1) to (2 * n - 8)
	for(int i = base_idx + 8; i < base_idx + 16; i++)
		hadamard(qubits, i);
	
	//Call to encryption oracle
	oracle_enc(&qubits, n, base_idx);

	//Last-round inversion
	add_key(&qubits, subkey, base_idx);
	s_box_inverse(&qubits, base_idx);

	for(int i = 4; i < req_qubits; i++)
		hadamard(qubits, i);

	//Printing superposition instead of measuring
	qreal prob1;
	for(int i = 0; i < (int)pow(2, req_qubits); i++)
	{
		prob1 = getProbAmp(qubits, i);
		if((int)(prob1 * 10000) == 0) continue;
		printf("Prob. of %02X %02X ", (i & 0xFF000) >> 12, (i & 0xFF0) >> 4);
		printf("linear characteristic is %f\n", prob1);
	}
	printf("\nSelect a suitable linear characteristic with sufficient high prob.");
	
	destroyQureg(qubits, env);
	destroyQuESTEnv(env);

	return 0;
}