package main

import (
	"fmt"
	"lpp"
	"math"
	"sort"
)

func COORD_SORT(array [][2]int, status int) [][2]int {
	cache_hash := make(map[int][][2]int)
	result := [][2]int{}
	cache_data := []int{}
	for _, data := range array {
		num := data[0]
		_, ok := cache_hash[num]
		if !ok {
			cache_hash[num] = [][2]int{}

		}
		cache_hash[num] = append(cache_hash[num], data)

	}
	for key, _ := range cache_hash {
		cache_data = append(cache_data, key)
	}
	if status == 0 {

		sort.Ints(cache_data)
	} else {

		sort.Sort(sort.Reverse(sort.IntSlice(cache_data)))
	}
	for _, element := range cache_data {
		for _, value := range cache_hash[element] {
			result = append(result, value)
		}
	}

	return result
}
func COORD_MERGE(array [][2]int, status int) int {
	array = COORD_SORT(array, status)

	length := 0
	result_list := [][2]int{}
	for i, data := range array {
		if i == 0 {
			result_list = append(result_list, data)
		} else {
			length := len(result_list) - 1
			if data[0] > data[1] {
				if data[0] >= result_list[length][1] {
					if data[1] < result_list[length][1] {
						result_list[length][1] = data[1]
					}
				} else {
					result_list = append(result_list, data)
				}
			} else {
				if data[0] <= result_list[length][1] {
					if data[1] >= result_list[length][1] {
						result_list[length][1] = data[1]
					}
				} else {
					result_list = append(result_list, data)
				}
			}
		}
	}

	for _, data := range result_list {
		length += int(math.Abs(float64(data[0]-data[1])) + 1)
	}
	return length
}
func COORD_CHAIN(array [][2]int, raw_array [][2]int, status int) ([][2]int, [][2]int, int) {
	//status 0是递增，status 1是递减

	i, j, k, max := 0, 0, 0, 0
	var result, result2 [][2]int
	length := len(array)

	//变长数组参数，C99新特性，用于记录当前各元素作为最大元素的最长递增序列长度
	liss := make([]int, length)

	//前驱元素数组，记录当前以该元素作为最大元素的递增序列中该元素的前驱节点，用于打印序列用
	pre := make([]int, length)

	for i = 0; i < length; i++ {

		liss[i] = 0
		pre[i] = i
	}

	max = 0
	k = 0
	for i = 1; i < length; i++ {
		//找到以array[i]为最末元素的最长递增子序列
		for j = 0; j < i; j++ {
			if status == 0 {
				if array[j][1] <= array[i][1] && array[i][0] < array[i][1] {
					cache := liss[j] + array[i][1] - array[i][0]

					if cache > liss[i] {
						liss[i] = cache
						if cache > max {
							max = cache
							k = i
						}
						if array[j][0] < array[j][1] {

							pre[i] = j
						}
					}
				}
			} else {
				if array[j][0] > array[i][0] && array[j][1] > array[i][1] && array[i][0] > array[i][1] {
					cache := liss[j] + array[i][0] - array[i][1]

					if cache > liss[i] {
						liss[i] = cache
						if cache > max {
							max = cache
							k = i
						}
						if array[j][0] > array[j][1] {
							pre[i] = j
						}
					}
				}

			}
			//如果要求非递减子序列只需将array[j] < array[i]改成<=，
			//如果要求递减子序列只需改为>

		}
	}
	//	fmt.Println(pre, k)
	//输出序列
	//	= max - 1

	result = [][2]int{array[k]}
	result2 = [][2]int{raw_array[k]}
	for {
		if pre[k] == k {
			break
		}

		k = pre[k]
		result = append([][2]int{array[k]}, result...)
		result2 = append([][2]int{raw_array[k]}, result2...)
	}
	//	fmt.Println(result)
	total_length := COORD_MERGE(result, status)
	return result, result2, total_length
}

func main() {
	data := [][2]int{
		{21, 7},
		{2, 1},
		{1, 2},
		{2, 141},
		{5, 9},
		{8, 9},
		{5, 11},
		{8, 7},
		{10, 9},
		{13, 21}, {14, 6}, {142, 227}, {139, 232}, {7, 5}, {5, 1}}

	fmt.Println(data)
	fmt.Println(lpp.COORD_CHAIN(data, data, 0))
	fmt.Println(COORD_CHAIN(data, data, 1))
}
