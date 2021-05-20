/*
8-bit block SPN Cipher based on Howard Heys' Design
Linear aprroximation relation: (0x10,0x0E), (0x04, 0xA0)
Linear cryptanalysis produces last-round's subkey's rightmost 4-bit
correctly with highest bias of 0.210938 and leftmost 4-bit
correctly with highest bias of -0.1328

Designed and coded by: ASHWINI KUMAR MALVIYA (ashwinixar@gmail.com)
*/

#include<iostream>
#include<stdlib.h>
#include<random>

using namespace std;

random_device rd;
mt19937 gen(rd());

typedef unsigned __int8 block;

//STARTS - Printing in Hex

void printHex(std::ostream& x, int num)
{
	ios::fmtflags f(x.flags());
	x << endl << std::hex << std::uppercase << num;
	x.flags(f);
}

//ENDS - Printing in Hex

//STARTS - Encryption Oracle Implementation

block k[5] = { 0xE3, 0x2B, 0xF5, 0x9A, 0x75 };

int s_table[16] = { 0xE,0x4,0xD,0x1,
					0x2,0xF,0xB,0x8,
					0x3,0xA,0x6,0xC,
					0x5,0x9,0x0,0x7 };

int p_table[8] = { 0,4,6,5,1,3,2,7 };

void s_box(block &a)
{
	block temp = 0;
	for (int i = 0; i < 2; i++)
	{
		int idx = (a & (0x0F << (i * 4))) >> (i * 4);
		temp = temp | (s_table[idx]) << (i * 4);
	}
	a = temp;
}

void s_box_inv(block &a)
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
	a = temp;
}

void p_box(block &a)
{
	block temp = 0;
	for (int i = 0; i < 8; i++)
	{
		if (a & (0x1 << i))
			temp = temp | (0x1 << p_table[i]);
	}
	a = temp;
}

block oracle_enc(block p)
{
	block c = p;
	for (int r = 0; r < 4; r++)
	{
		c = c ^ k[r];
		s_box(c);
		if (r != 3) p_box(c);
		else c = c ^ k[r + 1];
	}
	return c;
}

block oracle_dec(block c)
{
	block p = c;
	for (int r = 3; r >= 0; r--)
	{
		if (r != 3) p_box(p); //Self Invertible
		else p = p ^ k[r + 1];
		s_box_inv(p);
		p = p ^ k[r];
	}
	return p;
}

//ENDS - Encryption Oracle Implementation


//STARTS - Linear Cryptanalysis

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

int exists(block p, block a[], int n)
{
	for (int i = 0; i < n; i++)
		if (a[i] == p) return 1;
	return 0;
}

void lin_last_rounds_attack()
{
	//Data (Plaintext-Ciphertext pairs) generation and storage
	const int db = 256;
	block P[db], C[db];
	//block p;
	for (int i = 0; i < db; i++)
	{
		//p = gen() % 256;
		//while (exists(p, P, db)) p = gen() % 256;
		P[i] = i;// p;
		C[i] = oracle_enc(P[i]);
	}

	//Running for all target partial subkeys
	//Last 4-bits belong to the target partial subkey
	block subkeys[16];
	for (int i = 0; i < 16; i++)
		subkeys[i] = 0;
	for (int i = 0; i < 16; i++) //Loop for each subkey
	{
		for (int j = 0; j < db; j++) //Loop for each data (ciphertext) last round backward operation
		{
			block c = C[j] ^ i;
			s_box_inv(c);
			if ((bit_product(P[j], 0x10) ^ bit_product(c, 0x0E)) == 0)
				subkeys[i]++;
		}
	}

	cout << "\n\nLast 4-bits of the subkey:";
	cout << "\n--------------------------";
	for (int i = 0; i < 16; i++)
	{
		printHex(cout, (unsigned int)i);
		cout << ": " << (((db / 2.0) - subkeys[i]) / (float)db) << ", " << (int)subkeys[i] << ", " << (int)(128 - subkeys[i]);
	}

	//First 4-bits belong to the target partial subkey
	for (int i = 0; i < 16; i++)
		subkeys[i] = 0;
	for (int i = 0; i < 16; i++) //Loop for each subkey
	{
		for (int j = 0; j < db; j++) //Loop for each data (ciphertext) last round backward operation
		{
			block c = C[j] ^ (i << 4);
			s_box_inv(c);
			if ((bit_product(P[j], 0x04) ^ bit_product(c, 0xA0)) == 0)
				subkeys[i]++;
		}
	}

	cout << "\n\nFirst 4-bits of the subkey:";
	cout << "\n---------------------------";
	for (int i = 0; i < 16; i++)
	{
		printHex(cout, (unsigned int)i);
		cout << ": " << (((db / 2.0) - subkeys[i]) / (float)db) << ", " << (int)subkeys[i] << ", " << (int)(128 - subkeys[i]);
	}
}

//ENDS - Linear Cryptanalysis

int main()
{
	cout << "\nLinear Cryptanalysis:";
	lin_last_rounds_attack();
	cout << endl;

	return 0;
}