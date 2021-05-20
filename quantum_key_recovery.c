#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

typedef unsigned __int8 block;

//S-Box of the 8-bit custom SPN cipher
int s_table[16] = { 0xE,0x4,0xD,0x1,
					0x2,0xF,0xB,0x8,
					0x3,0xA,0x6,0xC,
					0x5,0x9,0x0,0x7 };

//Performs inverse S-Box
block s_box_inv(block a)
{
	block temp = 0;
	for (int i = 0; i < 2; i++)
	{
		int idx;
		block d = (a & (0x0F << (i * 4))) >> (i * 4);
		for (int i = 0; i < 16; i++)
			if (s_table[i] == d) idx = i;
		temp = temp | (idx) << (i * 4);
	}
	return temp;
}

int bit_product(block a, block mask)
{
	a = a & mask;
	mask = 0;
	while (a > 0)
	{
		mask = mask ^ (a & 0x1);
		a >>= 1;
	}
	return mask;
}

int oracle_e[256] = { 	0x3C,0xBD,0x22,0x5,0x15,0xEE,0xF0,0x68,0x64,0x72,0x88,0xBE,0xCC,0xD2,0xA0,0x5C,
						0x31,0xB9,0x3E,0x5E,0x9D,0xE3,0xFB,0x11,0xE1,0x9E,0x7B,0xB3,0xC1,0xFD,0xAB,0x51,
						0xF,0x82,0x37,0x57,0x9A,0xC5,0x26,0xAF,0x8F,0x97,0xC6,0xA2,0x46,0xFA,0x4F,0xD6,
						0x49,0x9,0x7C,0xB8,0x2C,0x81,0xA1,0xE8,0xDC,0x0,0xB4,0xCB,0x23,0x10,0xD3,0x43,
						0x4,0xEC,0x5A,0xF7,0xA7,0x70,0x28,0x66,0x6F,0xCA,0x86,0x1C,0x48,0x3A,0x44,0xD8,
						0x17,0xDA,0x4E,0xAC,0xDE,0x69,0x63,0xC0,0xD,0x2E,0x8C,0x83,0x7,0x4D,0x2A,0x7A,
						0x3D,0xB1,0x2F,0x6,0x16,0xEB,0xFE,0x65,0x62,0x7F,0x85,0xBB,0xCD,0xDF,0xAE,0x5D,
						0x94,0xF4,0x47,0xA9,0xD7,0xE2,0x12,0xC3,0xA,0x27,0x89,0x75,0x58,0x4A,0x38,0x98,
						0x1B,0xD1,0xC2,0xB6,0x52,0x6D,0x6E,0xE6,0x32,0xF5,0xBF,0x8E,0xB,0xA5,0x21,0x71,
						0x1E,0xDD,0xC4,0xB5,0x54,0x6C,0x60,0xE5,0x34,0xF8,0xB2,0x80,0xE,0xA8,0x2D,0x7D,
						0x2,0x84,0x33,0x53,0x99,0xC8,0x25,0x1A,0xEA,0x93,0x77,0xA4,0x45,0xF9,0x42,0xD5,
						0x13,0xD9,0x4B,0xAD,0xDB,0x61,0x6B,0xCE,0x1,0x2B,0x8D,0x8B,0x3,0x41,0x29,0x79,
						0x92,0xF2,0x96,0xAA,0x36,0xEF,0x1F,0xC7,0xFF,0x56,0x8A,0x76,0x55,0x9F,0x35,0x95,
						0xF1,0xE9,0x24,0x8,0x18,0x73,0x5B,0xB0,0xBC,0x74,0xE0,0x19,0x9B,0xD4,0x91,0x3B,
						0x4C,0xC,0xCF,0x67,0x5F,0xE4,0x14,0x87,0x3F,0xF6,0x6A,0x78,0x20,0xA6,0xD0,0x40,
						0x39,0xBA,0x30,0x50,0x9C,0xE7,0xF3,0x1D,0xED,0x90,0x7E,0xB7,0xC9,0xFC,0xA3,0x59 };

double est_eps(int j_prime, int a_0, int b_3, int shift)
{
	QuESTEnv env = createQuESTEnv();

    //Value of Pi
    const double pi = 3.14;

    const int p = 8, P = (int)pow(2, p); //Precision qubits for theta of Grover's iteration
    const int n = 8, N = 256;
    Qureg qubits = createQureg((p + n), env);
    initZeroState(qubits);
    //reportQuregParams(qubits);
    //reportQuESTEnv(env);

	ComplexMatrixN e = createComplexMatrixN(n);
    for(int i = 0; i < 256; i++)
    {
    	block c_i_j = oracle_e[i] ^ (j_prime << (shift * 4));
		c_i_j = s_box_inv(c_i_j);
		if ((bit_product(i, a_0) ^ bit_product(c_i_j, b_3)) == 0) e.real[i][i] = -1;
    	else e.real[i][i] = 1;
    }

    for(int i = 0; i < (p + n); i++)
        hadamard(qubits, i);
    
    int targs[n];
    for(int i = 0; i < n; i++)
    	targs[i] = p + i;

    //Applying controlled grover's operator
    int iterations = 1;
    for(int gi = p - 1; gi >= 0; gi--)
    {
    	int ctrls[] = { gi };
    	for(int j = 0; j < iterations; j++)
    	{
    		//Marking
	    	multiControlledMultiQubitUnitary(qubits, ctrls, 1, targs, n, e);
	    	
	        //Diffusion
	        for(int i = p; i < (n + p); i++)
	            hadamard(qubits, i);
	        for(int i = p; i < (n + p); i++)
	            pauliX(qubits, i);
	        multiControlledPhaseFlip(qubits, targs, n);
	        for(int i = p; i < (n + p); i++)
	            pauliX(qubits, i);
	        for(int i = p; i < (n + p); i++)
	            hadamard(qubits, i);
    	}
    	iterations *= 2;
    }
	
    //Inverse QFT on the first p-qubits
    for(int i = 0; i < p; i++)
    {
    	for(int j = 0; j < i; j++)
    	{
    		qreal angle = -2.0 * pi / (pow(2, i - j + 1));
    		controlledPhaseShift(qubits, j, i, angle);
    	}
    	hadamard(qubits, i);
    }
    
    //Measures the qubits to get the value of first p-qubits
    qreal prob = 0, max_prob = 0;
    int p_val = 0;
    for(int i = 0; i < P; i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(prob > max_prob)
    	{
    		max_prob = prob;
    		p_val = i;
    	}
    }
    if(p_val > (P / 2)) p_val = P - p_val;
    //printf("\nThe measured value of (f_tilde) is %d", p_val);
    //printf("\nThe approximate count of the correct answers is: %f\n", (N * pow(sin(p_val * pi / P), 2)));

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

	double eps = (((N * pow(sin(p_val * pi / P), 2)) - 128.0) / 256.0);
	if(eps < 0.0) eps = (-1) * eps;

    return eps;
}

int main (int narg, char *varg[])
{
	int lin_ch[2][2] = { { 0x10, 0x0E }, { 0x04, 0xA0 } };
	int subkey[2] = { 0x0, 0x0 } ;

	//Value of Pi
    const double pi = 3.14;

	int m = 2;
	int c = 4;
	int s_n = 4;

	double eps[2][16] = { { 0.012662, 0.036367, 0.024919, 0.121043,
							0.036367, 0.214036, 0.011866, 0.073722,
							0.073722, 0.072934, 0.060779, 0.048588,
							0.024124, 0.061570, 0.109107, 0.073722 },
						  { 0.049380, 0.049380, 0.097887, 0.049380,
							0.072934, 0.024124, 0.024919, 0.132907,
							0.037161, 0.012662, 0.085830, 0.012662,
							0.097106, 0.000398, 0.012662, 0.109107 } };
	

	for(int repeat = 0; repeat < 2; repeat++)
	{
		printf("\n-----------------------------------------------");
		printf("\nRunning for linear characteristic (0x%02X, 0x%02X):", lin_ch[repeat][0], lin_ch[repeat][1]);
		printf("\n-----------------------------------------------");

		int j_prime = 0x1; // A Random threshold key to begin with
		double est_eps_prime = est_eps(j_prime, lin_ch[repeat][0], lin_ch[repeat][1], repeat);

		printf("\nStart with threshold subkey 0x%X with estimated epsilon %f\n", j_prime, est_eps_prime);
		
		for(int times = 1; times <= 2 * m * c; times++)
		{
			QuESTEnv env = createQuESTEnv();
			Qureg qubits = createQureg(s_n, env);
			initZeroState(qubits);
			//reportQuregParams(qubits);
			//reportQuESTEnv(env);

			int ti = 0;
			ComplexMatrixN e = createComplexMatrixN(s_n);
			for(int i = 0; i < 16; i++)
			{
				if(eps[repeat][i] > est_eps_prime) {e.real[i][i] = -1; ti++;}
				else e.real[i][i] = 1;
			}

			for(int i = 0; i < s_n; i++)
				hadamard(qubits, i);
			int targs[] = { 0,1,2,3 };
			ti = (pi / 4.0) * sqrt(16.0 / ti);
			for(int gi = 0; gi < ti; gi++)
			{
				//Marking
				multiQubitUnitary(qubits, targs, s_n, e);
				
				//Diffusion
				for(int i = 0; i < s_n; i++)
					hadamard(qubits, i);
				for(int i = 0; i < s_n; i++)
					pauliX(qubits, i);
				multiControlledPhaseFlip(qubits, targs, s_n);
				for(int i = 0; i < s_n; i++)
					pauliX(qubits, i);
				for(int i = 0; i < s_n; i++)
					hadamard(qubits, i);
			}
			qreal prob;
			int j_m = 0;
			for(int i = 0; i < s_n; i++)
			{
				int outcome = measureWithStats(qubits, i, &prob);
				if(outcome) j_m ^= (outcome << i);
			}

			double temp = est_eps(j_m, lin_ch[repeat][0], lin_ch[repeat][1], repeat);
			printf("New measured subkey is 0x%X with estimated epsilon %f\n", j_m, temp);
			if(temp > est_eps_prime)
			{
				j_prime = j_m;
				est_eps_prime = temp;
				printf("Updated the subkey to 0x%X with bias %f\n\n", j_prime, est_eps_prime);
			}
			else
			{
				printf("Subkey remains same 0x%X with bias %f\n\n", j_prime, est_eps_prime);
			}

			destroyComplexMatrixN(e);
			destroyQureg(qubits, env);
			destroyQuESTEnv(env);
		}
		subkey[repeat] = j_prime;
		printf("The partial subkey with largest bias %f is 0x%X\n", est_eps_prime, subkey[repeat]);
	}

	printf("\nThe complete recovered last round key is 0x%X%X\n", subkey[1], subkey[0]);

    return 0;
}