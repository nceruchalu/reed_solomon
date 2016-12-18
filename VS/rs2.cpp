/************************************************************************
* @file		rs2.cpp
* @date		2016-12-18
* @author	yuan_yuanxiang@163.com
* @brief	RS2���������ԣ����У�2���ر�ʾһ������
* 00		0
* 01		1
* 10		2
* 11		3
************************************************************************/

#include "../reedSolomon.h"
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#define NUM 8

int main()
{
	int temp[3], result[3 * NUM] = { 0 };
	int array[NUM] = {3, 1, 2, 1, 2, 3, 1, 2};
	reedSolomon rs(2, 1);

	// ���ԭʼ����
	printf("Origin:\n");
	for (int i = 0; i < NUM; ++i)
		printf("%d\t", array[i]);
	printf("\n\n");

	// ��ԭ���ݽ��б���
	for (int j = 0; j < NUM; ++j)
	{
		memset(temp, 0, 3 * sizeof(int));
		temp[0] = array[j];
		rs.rs_encode(temp);
		memcpy(result + 3 * j, temp, 3 * sizeof(int));
	}
	// ����������������
	printf("Eecode:\n");
	for (int i = 0; i < 3 * NUM; ++i)
		printf("%d, ", result[i]);
	printf("\n\n");

	// ��Ӵ��󵽱�������
	printf("Error:\n");
	memset(result, 0, 16 * sizeof(int));
	for (int i = 0; i < 3 * NUM; ++i)
		printf("%d, ", result[i]);
	printf("\n\n");

	// �Դ�������ݽ��о���
	memset(array, 0, NUM * sizeof(int));
	for (int j = 0; j < NUM; ++j)
	{
		memcpy(temp, result + 3 * j, 3 * sizeof(int));
		rs.rs_decode(temp);
		if(false == rs.compare())
		{
			printf(" * Error occurs at %d.\n\n", j);
		}
		array[j] = temp[2];
	}

	// ���ݾ������Ľ��
	printf("Decode:\n");
	for (int i = 0; i < NUM; ++i)
		printf("%d\t", array[i]);
	printf("\n\n");

	system("pause");
	return 0;
}