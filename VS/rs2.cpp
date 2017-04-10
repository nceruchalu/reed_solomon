/************************************************************************
* @file		rs2.cpp
* @date		2016-12-18
* @author	yuan_yuanxiang@163.com
* @brief	RS2纠错编码测试，其中，2比特表示一个符号
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

	// 输出原始数据
	printf("Origin:\n");
	for (int i = 0; i < NUM; ++i)
		printf("%d\t", array[i]);
	printf("\n\n");

	// 对原数据进行编码
	for (int j = 0; j < NUM; ++j)
	{
		memset(temp, 0, 3 * sizeof(int));
		temp[0] = array[j];
		rs.rs_encode(temp);
		memcpy(result + 3 * j, temp, 3 * sizeof(int));
	}
	// 输出纠错编码后的数据
	printf("Eecode:\n");
	for (int i = 0; i < 3 * NUM; ++i)
		printf("%d, ", result[i]);
	printf("\n\n");

	// 添加错误到编码数据
	printf("Error:\n");
	memset(result, 0, 16 * sizeof(int));
	for (int i = 0; i < 3 * NUM; ++i)
		printf("%d, ", result[i]);
	printf("\n\n");

	// 对错误的数据进行纠错
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

	// 数据纠错编码的结果
	printf("Decode:\n");
	for (int i = 0; i < NUM; ++i)
		printf("%d\t", array[i]);
	printf("\n\n");

	system("pause");
	return 0;
}